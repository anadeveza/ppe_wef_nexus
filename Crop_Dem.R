library(gsubfn)
library(xts)
library(dplyr)
library(hydroTSM)
library(lubridate)
library(data.table)
library(MASS)
library(ggpubr)
library(prodlim)
library(Hmisc)

# install.packages(c("gsubfn", "xts", "dplyr", "hydroTSM", "lubridate",
#                    "data.table", "v", "hydroTSM", "ggpubr", "prodlim", "Hmisc"))

basePath <- c("D:\\Dropbox (PSR)\\Mestrado\\Modelo\\Raiz\\")
resultsPath <- c("D:\\Dropbox (PSR)\\Mestrado\\Modelo\\Raiz\\Resultados\\")

# basePath <- c("C:\\Users\\User\\Dropbox\\Mestrado\\Modelo\\Raiz\\")
# resultsPath <- c("C:\\Users\\User\\Dropbox\\Mestrado\\Modelo\\Raiz\\Resultados\\")

#Leitura dos arquivos base
postos_temp_min <- as.data.frame(fread(paste(basePath, "temp_min.csv", sep="") , sep = ",")) 
postos_temp_max <- as.data.frame(fread(paste(basePath, "temp_max.csv", sep="") , sep = ",")) 
# postos_temp_med <- as.data.frame(fread(paste(basePath, "temp_med.csv", sep="") , sep = ",")) 
postos_ET0 <- as.data.frame(fread(paste(basePath, "temp_max.csv", sep="") , sep = ","))  #Para aproveitar a estrutura do arquivo
postos_ET0[,-(1)] <-0
postos_prec <- as.data.frame(fread(paste(basePath, "precipitacao.csv", sep="") , sep = ",")) 
######## Garantir que a ordem dos postos nesses 5 primeiros arquivos seja igual
postos_Pef <-as.data.frame(fread(paste(basePath, "precipitacao.csv", sep="") , sep = ",")) #Para aproveitar o formato
postos_Pef[,-(1)] <-0
localidades <- as.data.frame(fread(paste(basePath, "localidades.csv", sep="") , sep = ",")) 
areas <-as.data.frame(fread(paste(basePath, "areas_calendario.csv", sep="") , sep = ",")) 
culturas <- as.data.frame(fread(paste(basePath, "culturas.csv", sep="") , sep = ",")) 
culturas_headers <- colnames(culturas)
lotes <- as.data.frame(fread(paste(basePath, "lotes.csv", sep="") , sep = ",")) # sao as areas de cultivo
lotes_headers <- colnames(lotes)
# eficiencias <- as.data.frame(fread(paste(basePath, "eficiencias.csv", sep="") , sep = ",")) #Ainda nao esta sendo usada
latitudes <- read.csv("estacoes_latitude.csv", header = TRUE, sep = ",", quote = "\"",
                      dec = ".", stringsAsFactors = FALSE)
# coef_referencia <- as.data.frame(fread(paste(basePath, "coef_referencia.csv", sep="") , sep = ","))
coef_retirada <- as.data.frame(fread(paste(basePath, "areas_calendario.csv", sep="") , sep = ",")) #Para aproveitar a estrutura
coef_retirada[,-(1)] <-0
Q_demandada <-as.data.frame(fread(paste(basePath, "areas_calendario.csv", sep="") , sep = ","))  #Para aproveitar a estrutura
Q_demandada[,-(1)] <-0
ETRc <- as.data.frame(fread(paste(basePath, "areas_calendario.csv", sep="") , sep = ",")) # sao as areas de cultivo
ETRc[,-(1)] <-0
IRA <- as.data.frame(fread(paste(basePath, "areas_calendario.csv", sep="") , sep = ",")) # sao as areas de cultivo
IRA[,-(1)] <-0


Evap_HS <- function(tmin , tmax , Lat){
  N = length(tmin)
  tmean = (tmax+tmin)/2
  #Conversao de latitude de grau pra radianos
  radLat = Lat*pi/180
  #Number of the day in the year
  nDay = rep(c(15,45,76,106,137,167,197,228,258,289,319,349),N/12)
  #Inverse relative distance Earth-Sun
  dr = 1+0.033*cos(2*pi*nDay/365)
  #Solar declination
  sigma = 0.409*sin(2*pi*nDay/365-1.39)
  #Sunset hour angle
  w = acos(-tan(radLat)*tan(sigma))
  gsc = 0.082
  #Extraterrestrial radiation (mm/day)
  Ra = 24*60/pi*gsc*dr*(w*sin(radLat)*sin(sigma)+sin(w)*cos(radLat)*cos(sigma))*0.408
  #Evapotranspiracao anual
  ET = 0.0023*(tmean+17.8)*(tmax-tmin)^(0.5)*Ra
  #Evapotranspiracao mensal
  ET = ET*365/12
  return(ET)
}

#Calculo da ET0
nPostos <- length(colnames(postos_ET0))-1

  for(k in 1:nPostos){
    station <- colnames(postos_ET0)[k+1]
    Lat <- latitudes[latitudes$Código_INMET==station,2]
    tmax <- as.numeric(postos_temp_max[,k+1])
    tmin <- as.numeric(postos_temp_min[,k+1])
    postos_ET0[,k+1] <- Evap_HS(tmin , tmax , Lat)
  }

#Calculo da precipitacao efetiva (Pef)
nPostos <- length(colnames(postos_prec))-1

  for(j in 1:nPostos){
    
    a <- as.numeric(postos_prec[,j+1])
    b <- as.numeric(postos_Pef[,j+1]) #Para deixar b na estrutura
      
    b <- ifelse(a >250, b <- a - 0.2/125*a^2, b <- 125 + 0.1*a)
     
    postos_Pef[,j+1] <- b
  
  }

####Conferir se a conta do posto certo esta sendo colocada na colunas certas dos arquivos finais  
nLotes <- length(colnames(ETRc))-1
  for(lote in 1:nLotes){  
    #Le o codigo do lote
    codigo <- colnames(ETRc)[1+lote]
    #Le a cultura
    cultura_AC <- lotes[2,match(codigo,lotes_headers)] #Indica a linha das culturas e a coluna no arquivo de lotos referente ao codigo em questao
    #Le o Ks  
    Ks <- as.numeric(lotes[5,match(codigo,lotes_headers)])
    #Le o Kc
    Kc_cultura <- culturas[10,match(cultura_AC,colnames(culturas))] #Calcular o ponderado e chamar a linha 10!
    #Le a localidade
    localidade_AC <- lotes[1,match(codigo,lotes_headers)] #Indica a linha das localidades e a coluna no arquivo de lotos referente ao codigo em questao
    #Le o posto de ET0 da localidade
    posto_AC <-localidades[4,match(toupper(localidade_AC),toupper(colnames(localidades)))]# toupper muda para lowercase
    #Le a serie do posto
    seriesET0_AC <- as.numeric(postos_ET0[,match(posto_AC,colnames(postos_ET0))])
    ETRc[1+lote] <-  seriesET0_AC
  
    #Identifica posto de Pef da localidade
    series_Pef_AC <- postos_Pef[,match(posto_AC,colnames(postos_Pef))]
    #Ler a Ea da AC
    Ea_AC <- as.numeric(lotes[4,match(codigo,lotes_headers)])
    #Ler a Tu da AC
    Tu_AC <- as.numeric(lotes[6,match(codigo,lotes_headers)])
    #Calcula IRA
    series_IRA_AC <- (seriesET0_AC*Ks*Kc_cultura-series_Pef_AC)*Tu_AC/Ea_AC
    #Criar df com a serie de IRA
    IRA[1+lote] <- series_IRA_AC
    
    #Atencao - se quem limita as datas das demandas eh o tamanho do arq de area ou de postos
    filtro_IRA <- IRA[match(areas[,1], IRA$data, nomatch=0),]
    datas <- filtro_IRA[,1]
    nDias <- monthDays(as.Date(datas, format="%m/%d/%Y"))
    IRA_mes <- filtro_IRA[match(codigo,colnames(filtro_IRA))]
    filtro_area <- areas[match(datas, areas$data, nomatch=0),]
    area_mes <- filtro_area[match(codigo,colnames(filtro_area))]
    coef_retirada[,1+lote] <- IRA_mes[,1]*10000/(nDias*24*3600)
    colnames(coef_retirada)[1+lote] <- codigo
    Q_demandada[,1+lote] <- IRA_mes[,1]*area_mes[,1]*10000/(nDias*24*3600*1000)
    colnames(Q_demandada)[1+lote] <- codigo
  }    


#Para garantir que a ordem das colunas seja igual
# b <- coef_referencia
# Cod_ref <- length(colnames(coef_referencia))-1
#   for(Cod in 1:Cod_ref){  
#     codigo <- colnames(coef_referencia)[1+Cod]
#     localidade_2 <- coef_retirada[match(codigo,colnames(coef_retirada))]
#     b[,1+Cod] <- coef_retirada[match(codigo,colnames(coef_retirada))]
#   } 
# 
# referencia <- as.matrix(coef_referencia[(-1)])
# calculado <- as.matrix(b[(-1)])
# correlacao <- cor(referencia,calculado)
# diagonal <- diag(correlacao)
# 
# #Plota coeficiente de referencia x coeficiente calculado
# png('correlacoes.png')
# matplot(referencia,calculado,xlim=c(0,2),ylim=c(0,2))
# # plot(diagonal)
# dev.off()

#Avaliar os clusters da mesma cultura e da mesma localidade
# #Diminuir o numero de linhas do bloco abaixo
# x <- colnames(lotes)[(-1)]
# y <- t(lotes[1,(-1)])
# rownames(y) <- NULL
# df <- data.frame(cbind(x,y),stringsAsFactors=FALSE)
# colnames(df) <- c("lote","localidade")
# lotes_labels <- data.frame(rownames(as.data.frame(diagonal)))
# colnames(lotes_labels) <- c("lote")
# lote_localidade <- left_join(lotes_labels, df, by = c("lote")) 

# write.csv(lote_localidade, file="lote_localidade.csv",row.names=FALSE)  
# write.csv(diagonal,  paste(resultsPath,"correlacoes.csv", sep=""), sep=",",row.names=TRUE)  
write.csv(postos_ET0,  paste(resultsPath,"ET0_postos.csv", sep=""), sep=",",row.names=TRUE)  
write.csv(postos_Pef,  paste(resultsPath,"Pef.csv", sep=""), sep=",",row.names=TRUE)  
write.csv(ETRc,  paste(resultsPath,"ETRc.csv", sep=""), sep=",",row.names=TRUE)  
write.csv(IRA,  paste(resultsPath,"IRA.csv", sep=""), sep=",",row.names=TRUE)  
write.csv(coef_retirada,  paste(resultsPath,"coef_retirada.csv", sep=""), sep=",",row.names=TRUE)  
write.csv(Q_demandada,  paste(resultsPath,"Q_demandada.csv", sep=""), sep=",",row.names=FALSE)  
  