##S�rie: taxa mensal IGP-M jan 2000 - julho 2018

setwd('C:\\Users\\godinho\\Desktop\\projeto ipea')
dados <- read.csv2("igpm.csv", header = TRUE,sep = ',', dec = '.')


# Transformando em s�rie temporal e selecionando o per�odo

library(tseries)
igpm <- ts(dados[,2][127:349], freq = 12, start = c(2000,1), end = c(2018, 7))


# N�mero �ndice (100 =  jan 2016)

data <- data.frame(igpm)

indice <- data[1:223,]*100/data[193,]

indice_ts <- ts(indice, freq = 12, start = c(2000,1), end = c(2018, 7))



ts.plot(indice_ts, col = c("darkred"),
        gpars = list(ylab="�ndice", ylim = c(-100,400), xlab="ano",lwd=1, bty='l',lty=1 ))
legend("topright", cex=.7 ,legend=c("�ndice Jan 2016 = 100"), fill = c("darkred"),bty="n")
axis(1, 2018)
grid(col='lightgray', lwd =.05, lty = 2)
dev.off()





# Estat�siticas descritivas
library(moments)

summary(igpm)
sd(igpm)
coef.var=sd(igpm)/mean(igpm)*100
skewness(igpm)
kurtosis(igpm)
# Gr�ficos dos dados brutos


#Pacotepara integra��o do R com o Latex
require('tikzDevice')
#tikz(file="igpm_bruto.tex", width = 6, height = 3.5)
ts.plot(igpm, col = c("darkred"),
        gpars = list(ylab="\\%", xlab="ano",lwd=1, bty='l',lty=1 ))
legend("topright", cex=.7 ,legend=c("IGP-M"), fill = c("darkred"),bty="n")
axis(1, 2018)
grid(col='lightgray', lwd =.05, lty = 2)
dev.off()



# Teste de nomalidade:
## Jarque-Bera

jarque.bera.test(igpm)
### A s�rie rejeita H0, ou seja, a s�rie n�o segue a normal.

## Confirma��o: Histograma

higpm<-hist(igpm, breaks=10, col="gray", xlab="", ylab="") 
igpmfit<-seq(min(igpm),max(igpm),length=40) 
yigpmfit<-dnorm(igpmfit,mean=mean(igpm),sd=sd(igpm)) 
yigpmfit <- yigpmfit*diff(higpm$mids[1:2])*length(igpm) 
lines(igpmfit, yigpmfit, col="black", lwd=2)

dev.off()




# Decomposi��o gr�fica


#tikz(file="decompose.tex", width = 6, height = 3.5)
plot(decompose(igpm))
dev.off()

#Tratamento dos dados

#Verificar sazonalidade
##M�todo 1: An�lise gr�fica
###A princ�pio, uma an�lise gr�fica � poss�vel detectar que n�o existe sazonalidade
##M�todo 2: Gerar dummies mensais e regredir na s�rie
### O teste mostra que n�o h� sazonalidade determin�stica
#Gerar as dummies e regredir o igpm 
M <- ordered(cycle(igpm))

igpm.reg <- lm(igpm~M)

summary(igpm.reg)


#Criando dummies mensais, sazonalidade determin�stica
mensal <- factor(cycle(igpm))
mensal <- model.matrix(~mensal)[,-1]
# Regress�o das s�ries contra as dummies
igpm.saz <- lm(igpm~mensal)
summary(igpm.saz)

# Nenhuma das s�ries precisa de desazonaliza��o
# Confirma��o gr�fica: boxplot
boxplot(igpm ~ cycle (igpm), las = 1, names = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10", "M11", "M12"), horizontal = T)

## n�o h� sazonalidade


##Res�duos da regress�o: s�rie dessazonalizada(igpm.des), e o componente sazonal(igpm.sazo)
igpm.des <- ts(resid(igpm.reg), start = 2000, freq=12)
igpm.sazo <- ts(fitted(igpm.reg), start=2000, freq=12)
par(mfrow=c(1,2))
plot(igpm.des)
plot(igpm.sazo)

#normalizar a serie, para retirar os valores negativos

igpm.desn <- igpm.des + mean(fitted(igpm.reg))

#gr�fico igpm vs igpm dessazonalizado

#tikz(file="sazonalidade.tex", width = 6, height = 3.5)
par(mfrow=c(1,1))
plot(igpm,
     main='',
     xlab='Ano', ylab='',
     col='blue',
     bty='l')
par(new=TRUE)
plot(igpm.desn,
     axes=F, ann=F,
     col='red',
     lty=2)
legend('topright',
       c('IGP-M', 'IGP-M dessazonalizado'),
       col=c('blue', 'red'), lty=1:2,
       bty='n')
axis(1,2018)
grid(col='darkgrey')

dev.off()

#Verificar tend�ncia
##Analisando os gr�ficos das s�ries � razo�vel supor que existe tend�ncia.

##M�todo 1: Retirar a tend�ncia: regredindo cada s�rie contra o tempo e recolhendo os res�duos
tempo <- time(igpm)
igpm.td <- lm(igpm ~ tempo)
linigpm <- residuals(igpm.td)

##Convertendo res�duos para s�ries temporais sem tend�ncia
linigpm <- ts(linigpm, frequency = 12, start = c(2000,1))


##M�todo 2: Filtro hp
filtrohp <- function(y, lambda){
  id <- diag(length(y))
  d <- diff(id, d=2)
  tendhp <- solve(id+lambda*crossprod(d), y)
  tendhp
}

lambda <- 14400
tendencia <- filtrohp(igpm, lambda)
hpigpm <- igpm - tendencia



# apresentando graficamente os dois m�todos e a s�rie em n�vel  
#tikz(file="tendencia.tex", width = 6, height = 3.5)
par(mfrow=c(1,1))
plot(linigpm,
     main='',
     xlab='Ano', ylab='',
     col='blue',
     bty='l')
par(new=TRUE)
plot(hpigpm,
     axes=F, ann=F,
     col='red',
     lty=2)
par(new=TRUE)
plot(igpm,
     axes=F, ann=F,
     col='green',
     lty=3)

legend('topright',
       c('Linear', 'Filtro HP','IGP-M'),
       col=c('blue', 'red', 'green'), lty=1:3,
       bty='n')
axis(1, 2018)
grid(col='darkgrey')

dev.off()


## Supondo uma tend�ncia determin�stica para a s�rie, � escolhido o m�todo linear,
##mais simples que o filtro HP, por�m eficiente nesse tipo de caso.

##Detectando e tratando Outliers
## A partir de uma an�lise gr�fica � poss�vel verificar outliers
## A op��o de tratamento foi pegar a s�rie a partir de 2005
ig <- window(igpm, start = c(2005,01))
#tikz(file="igpm.tex", width = 6, height = 3.5)
ts.plot(ig, col = c("darkred"),
        gpars = list(ylab="\\%", xlab="ano",lwd=1, bty='l',lty=1 ))
legend("topright", cex=.7 ,legend=c("IGP-M"), fill = c("darkred"),bty="n")
grid(col='lightgray', lwd =.05, lty = 2)
dev.off()




#Estima��o do modelo de previs�o ARIMA (p,d,q)
##Metodologia Box-Jekins

##1� Passo:
### Determinar a defasagem (d)
### Testes de ra�z unit�ria
####Dickey Fuller
summary(lm(diff(ig)~lag(ig,-1)[-length(ig)] - 1))
##Rejeita a hip�tese nula
####Dickey Fuller Aumentado Ho = possui ra�z unit�ria
library(urca)
##Sem tend�ncia
adf <- ur.df(ig, type = "none", selectlags = "AIC")
summary(adf)
## Com tend�ncia
adf1 <- ur.df(ig, type = "trend", selectlags = "AIC")
summary(adf1)
## Com constante
adf2 <- ur.df(ig, type='drift', selectlags = "AIC")
summary(adf2)
acf(adf2@res, ci.type='ma',
    main='Res�duos de teste ADF ',
    xlab='Defasagem')

## O teste rejeita a hip�tese nula de ra�z unit�ria, ou seja a s�rie � estacion�ria.

####Portanto, ser� trabalhado a s�rie em n�vel (d=0)

##2� Passo 
###Escolher p e q
### Fun��o de autocorrel��o e autocorrela��o parcial
par(mfrow=c(1,1))
acf(ig, ci.type='ma', main='IGP-M (varia��o mensal)')
pacf(ig, main='IGP-M (varia��o mensal)')

### Como � um m�todo mais visual, muitas vezes pode ser confuso
### Opto por partir de modelos mais gerais e testar a signific�ncia da estima��o
library(lmtest)



###Analisando o melhor modelo
###Diagn�stico dos res�duos
#1.Os res�duos da regress�o s�o n�o autocorrelacionados?
#2. Os res�duos s�o normalmente distribu�dos?
#3. Devemos incluir/exluir regressores
#Para responder a primeira quest�o � utilizado o teste Ljung-Box.
#Hip�tese nula a independ�ncia de uma dada s�rie de tempo, isto �, assume que as observa��es que temos
#(neste caso os res�duos de uma regress�o) s�o conjuntamente n�o correlacionados ao longo do tempo.

###Modelos escolhidos como poss�veis
###reg1.i = ARIMA(1,0,0)
###reg2.i = ARIMA(0,0,4)
###reg3.i = ARIMA(1,0,1)
###reg4.i = ARIMA(1,0,4)
#Para responder a primeira quest�o � utilizado o teste Ljung-Box.
#Hip�tese nula: h� independ�ncia de uma dada s�rie de tempo, isto �, assume que as observa��es que temos
#(neste caso os res�duos de uma regress�o) s�o conjuntamente n�o correlacionados ao longo do tempo.

reg1.i <- arima(ig, order = c(1,0,0))
Box.test(resid(reg1.i), lag=6, type='Ljung-Box', fitdf=1)

reg2.i <- arima(ig, order = c(0,0,4))
Box.test(resid(reg2.i), lag=6, type='Ljung-Box', fitdf=1)


reg3.i <- arima(ig, order = c(1,0,1))
Box.test(resid(reg3.i), lag=6, type='Ljung-Box', fitdf=1)

reg4.i <- arima(ig, order = c(1,0,4))
Box.test(resid(reg4.i), lag=6, type='Ljung-Box', fitdf=1)


###Testar se as s�ries s�o normalmente distribu�da
##Shapiro test, Hip�tese nula: a s�rie � normalmente distribu�da
shapiro.test(resid(reg1.i))

shapiro.test(resid(reg2.i))

shapiro.test(resid(reg3.i))         

shapiro.test(resid(reg4.i))

### Aceita a hip�tese nula nos  modelos, ou seja, s�o normalmente distribu�dos

###Escolhendo o melhor modelo pelo crit�rio Akaike
#Menor o crit�rio de informa��o melhor � seu poder de informa��o
AIC(reg1.i)
AIC(reg2.i)
AIC(reg3.i)
AIC(reg4.i)

BIC(reg1.i)
BIC(reg2.i)
BIC(reg3.i)
BIC(reg4.i)
#Pelo AIC o melhor modelo � o ARIMA(1,0,0)

###Outro m�todo de escolha do modelo � realizar uma previs�o dentro da amostra
###Analisar a soma do erro ao quadrado m�dio
###Previs�o dentro da amostra, horizonte igual a 12



#Avaliando os erros de previs�o


erro1 <- matrix(NA, nrow=length(window(ig, start=2005)), ncol=1)
for (i in 1:length(erro1)){erro1[i] <- ig[(length(window(ig, end=c(2017, 07)))+i)] -
  predict(arima(ig[1:(length(window(ig, end=c(2017, 07)))+i-1)], order=c(1,0,0)),
          n.ahead=1)$pred
erro1 <- ts(erro1, start=c(2017,08) , freq=12)
}

erro2 <- matrix(NA, nrow=length(window(ig, start=2005)), ncol=1)
for (i in 1:length(erro2)){erro2[i] <- ig[(length(window(ig, end=c(2017, 07)))+i)] -
  predict(arima(ig[1:(length(window(ig, end=c(2017, 07)))+i-1)], order=c(0,0,4)),
          n.ahead=1)$pred
erro2 <- ts(erro2, start=c(2017, 08), freq=12)
}

erro3 <- matrix(NA, nrow=length(window(ig, start=2005)), ncol=1)
for (i in 1:length(erro3)){erro3[i] <- ig[(length(window(ig, end=c(2017, 07)))+i)] -
  predict(arima(ig[1:(length(window(ig, end=c(2017, 07)))+i-1)], order=c(1,0,1)),
          n.ahead=1)$pred
erro3 <- ts(erro3, start=c(2017, 08), freq=12)
}

erro4 <- matrix(NA, nrow=length(window(ig, start=2005)), ncol=1)
for (i in 1:length(erro4)){erro4[i] <- ig[(length(window(ig, end=c(2017, 07)))+i)] -
  predict(arima(ig[1:(length(window(ig, end=c(2017, 07)))+i-1)], order=c(1,0,4)),
          n.ahead=1)$pred
erro4 <- ts(erro4, start=c(2017, 08), freq=12)
}
#Ajustando a janela da previs�o
erro1 <- window(erro1, start=c(2017, 08), end = c(2018, 07))
erro2 <- window(erro2, start=c(2017, 08), end = c(2018, 07))
erro3 <- window(erro3, start=c(2017, 08), end = c(2018, 07))
erro4 <- window(erro4, start=c(2017, 08), end = c(2018, 07))
#Calculando o menor erro ao quadrado m�dio

sqe1 <- (sum(erro1^2)/length(erro1))
sqe2 <- (sum(erro2^2)/length(erro2))
sqe3 <- (sum(erro3^2)/length(erro3))
sqe4 <- (sum(erro4^2)/length(erro4))
sqe1
sqe2
sqe3
sqe4
# A partir da an�lise do erro ao quadrado m�dio 
#a melhor previs�o tamb�m se encaixa no modelo ARIMA(1,0,1)
#previs�o para o m�s de julho
library(forecast)
ig.b <- window(ig, end=c(2018, 6))
reg.ig.b <- arima(ig.b, order=c(1,0,1))

prev1 <- forecast(reg.ig.b, h=1, level=c(25,50,90))
summary(prev1)


plot.forecast(prev1, include=24, 
              main='Previs�o do IGP-M com um modelo ARMA(1,1)')

#Previs�o para agosto
reg.ig.a <- arima(ig, order=c(0,0,4))
prev2 <- forecast(reg.ig.a, h=1,level=c(30,60))
summary(prev2)

plot.forecast(prev2, include=24, 
              main='Previs�o do IGP-M para dezembro com um modelo MA(4)')
