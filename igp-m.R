##Série: taxa mensal IGP-M jan 2000 - julho 2018

setwd('C:\\Users\\godinho\\Desktop\\projeto ipea')
dados <- read.csv2("igpm.csv", header = TRUE,sep = ',', dec = '.')


# Transformando em série temporal e selecionando o período

library(tseries)
igpm <- ts(dados[,2][127:349], freq = 12, start = c(2000,1), end = c(2018, 7))


# Número índice (100 =  jan 2016)

data <- data.frame(igpm)

indice <- data[1:223,]*100/data[193,]

indice_ts <- ts(indice, freq = 12, start = c(2000,1), end = c(2018, 7))



ts.plot(indice_ts, col = c("darkred"),
        gpars = list(ylab="Índice", ylim = c(-100,400), xlab="ano",lwd=1, bty='l',lty=1 ))
legend("topright", cex=.7 ,legend=c("Índice Jan 2016 = 100"), fill = c("darkred"),bty="n")
axis(1, 2018)
grid(col='lightgray', lwd =.05, lty = 2)
dev.off()





# Estatísiticas descritivas
library(moments)

summary(igpm)
sd(igpm)
coef.var=sd(igpm)/mean(igpm)*100
skewness(igpm)
kurtosis(igpm)
# Gráficos dos dados brutos


#Pacotepara integração do R com o Latex
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
### A série rejeita H0, ou seja, a série não segue a normal.

## Confirmação: Histograma

higpm<-hist(igpm, breaks=10, col="gray", xlab="", ylab="") 
igpmfit<-seq(min(igpm),max(igpm),length=40) 
yigpmfit<-dnorm(igpmfit,mean=mean(igpm),sd=sd(igpm)) 
yigpmfit <- yigpmfit*diff(higpm$mids[1:2])*length(igpm) 
lines(igpmfit, yigpmfit, col="black", lwd=2)

dev.off()




# Decomposição gráfica


#tikz(file="decompose.tex", width = 6, height = 3.5)
plot(decompose(igpm))
dev.off()

#Tratamento dos dados

#Verificar sazonalidade
##Método 1: Análise gráfica
###A princípio, uma análise gráfica é possível detectar que não existe sazonalidade
##Método 2: Gerar dummies mensais e regredir na série
### O teste mostra que não há sazonalidade determinística
#Gerar as dummies e regredir o igpm 
M <- ordered(cycle(igpm))

igpm.reg <- lm(igpm~M)

summary(igpm.reg)


#Criando dummies mensais, sazonalidade determinística
mensal <- factor(cycle(igpm))
mensal <- model.matrix(~mensal)[,-1]
# Regressão das séries contra as dummies
igpm.saz <- lm(igpm~mensal)
summary(igpm.saz)

# Nenhuma das séries precisa de desazonalização
# Confirmação gráfica: boxplot
boxplot(igpm ~ cycle (igpm), las = 1, names = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10", "M11", "M12"), horizontal = T)

## não há sazonalidade


##Resíduos da regressão: série dessazonalizada(igpm.des), e o componente sazonal(igpm.sazo)
igpm.des <- ts(resid(igpm.reg), start = 2000, freq=12)
igpm.sazo <- ts(fitted(igpm.reg), start=2000, freq=12)
par(mfrow=c(1,2))
plot(igpm.des)
plot(igpm.sazo)

#normalizar a serie, para retirar os valores negativos

igpm.desn <- igpm.des + mean(fitted(igpm.reg))

#gráfico igpm vs igpm dessazonalizado

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

#Verificar tendência
##Analisando os gráficos das séries é razoável supor que existe tendência.

##Método 1: Retirar a tendência: regredindo cada série contra o tempo e recolhendo os resíduos
tempo <- time(igpm)
igpm.td <- lm(igpm ~ tempo)
linigpm <- residuals(igpm.td)

##Convertendo resíduos para séries temporais sem tendência
linigpm <- ts(linigpm, frequency = 12, start = c(2000,1))


##Método 2: Filtro hp
filtrohp <- function(y, lambda){
  id <- diag(length(y))
  d <- diff(id, d=2)
  tendhp <- solve(id+lambda*crossprod(d), y)
  tendhp
}

lambda <- 14400
tendencia <- filtrohp(igpm, lambda)
hpigpm <- igpm - tendencia



# apresentando graficamente os dois métodos e a série em nível  
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


## Supondo uma tendência determinística para a série, é escolhido o método linear,
##mais simples que o filtro HP, porém eficiente nesse tipo de caso.

##Detectando e tratando Outliers
## A partir de uma análise gráfica é possível verificar outliers
## A opção de tratamento foi pegar a série a partir de 2005
ig <- window(igpm, start = c(2005,01))
#tikz(file="igpm.tex", width = 6, height = 3.5)
ts.plot(ig, col = c("darkred"),
        gpars = list(ylab="\\%", xlab="ano",lwd=1, bty='l',lty=1 ))
legend("topright", cex=.7 ,legend=c("IGP-M"), fill = c("darkred"),bty="n")
grid(col='lightgray', lwd =.05, lty = 2)
dev.off()




#Estimação do modelo de previsão ARIMA (p,d,q)
##Metodologia Box-Jekins

##1º Passo:
### Determinar a defasagem (d)
### Testes de raíz unitária
####Dickey Fuller
summary(lm(diff(ig)~lag(ig,-1)[-length(ig)] - 1))
##Rejeita a hipótese nula
####Dickey Fuller Aumentado Ho = possui raíz unitária
library(urca)
##Sem tendência
adf <- ur.df(ig, type = "none", selectlags = "AIC")
summary(adf)
## Com tendência
adf1 <- ur.df(ig, type = "trend", selectlags = "AIC")
summary(adf1)
## Com constante
adf2 <- ur.df(ig, type='drift', selectlags = "AIC")
summary(adf2)
acf(adf2@res, ci.type='ma',
    main='Resíduos de teste ADF ',
    xlab='Defasagem')

## O teste rejeita a hipótese nula de raíz unitária, ou seja a série é estacionária.

####Portanto, será trabalhado a série em nível (d=0)

##2º Passo 
###Escolher p e q
### Função de autocorrelção e autocorrelação parcial
par(mfrow=c(1,1))
acf(ig, ci.type='ma', main='IGP-M (variação mensal)')
pacf(ig, main='IGP-M (variação mensal)')

### Como é um método mais visual, muitas vezes pode ser confuso
### Opto por partir de modelos mais gerais e testar a significância da estimação
library(lmtest)



###Analisando o melhor modelo
###Diagnóstico dos resíduos
#1.Os resíduos da regressão são não autocorrelacionados?
#2. Os resíduos são normalmente distribuídos?
#3. Devemos incluir/exluir regressores
#Para responder a primeira questão é utilizado o teste Ljung-Box.
#Hipótese nula a independência de uma dada série de tempo, isto é, assume que as observações que temos
#(neste caso os resíduos de uma regressão) são conjuntamente não correlacionados ao longo do tempo.

###Modelos escolhidos como possíveis
###reg1.i = ARIMA(1,0,0)
###reg2.i = ARIMA(0,0,4)
###reg3.i = ARIMA(1,0,1)
###reg4.i = ARIMA(1,0,4)
#Para responder a primeira questão é utilizado o teste Ljung-Box.
#Hipótese nula: há independência de uma dada série de tempo, isto é, assume que as observações que temos
#(neste caso os resíduos de uma regressão) são conjuntamente não correlacionados ao longo do tempo.

reg1.i <- arima(ig, order = c(1,0,0))
Box.test(resid(reg1.i), lag=6, type='Ljung-Box', fitdf=1)

reg2.i <- arima(ig, order = c(0,0,4))
Box.test(resid(reg2.i), lag=6, type='Ljung-Box', fitdf=1)


reg3.i <- arima(ig, order = c(1,0,1))
Box.test(resid(reg3.i), lag=6, type='Ljung-Box', fitdf=1)

reg4.i <- arima(ig, order = c(1,0,4))
Box.test(resid(reg4.i), lag=6, type='Ljung-Box', fitdf=1)


###Testar se as séries são normalmente distribuída
##Shapiro test, Hipótese nula: a série é normalmente distribuída
shapiro.test(resid(reg1.i))

shapiro.test(resid(reg2.i))

shapiro.test(resid(reg3.i))         

shapiro.test(resid(reg4.i))

### Aceita a hipótese nula nos  modelos, ou seja, são normalmente distribuídos

###Escolhendo o melhor modelo pelo critério Akaike
#Menor o critério de informação melhor é seu poder de informação
AIC(reg1.i)
AIC(reg2.i)
AIC(reg3.i)
AIC(reg4.i)

BIC(reg1.i)
BIC(reg2.i)
BIC(reg3.i)
BIC(reg4.i)
#Pelo AIC o melhor modelo é o ARIMA(1,0,0)

###Outro método de escolha do modelo é realizar uma previsão dentro da amostra
###Analisar a soma do erro ao quadrado médio
###Previsão dentro da amostra, horizonte igual a 12



#Avaliando os erros de previsão


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
#Ajustando a janela da previsão
erro1 <- window(erro1, start=c(2017, 08), end = c(2018, 07))
erro2 <- window(erro2, start=c(2017, 08), end = c(2018, 07))
erro3 <- window(erro3, start=c(2017, 08), end = c(2018, 07))
erro4 <- window(erro4, start=c(2017, 08), end = c(2018, 07))
#Calculando o menor erro ao quadrado médio

sqe1 <- (sum(erro1^2)/length(erro1))
sqe2 <- (sum(erro2^2)/length(erro2))
sqe3 <- (sum(erro3^2)/length(erro3))
sqe4 <- (sum(erro4^2)/length(erro4))
sqe1
sqe2
sqe3
sqe4
# A partir da análise do erro ao quadrado médio 
#a melhor previsão também se encaixa no modelo ARIMA(1,0,1)
#previsão para o mês de julho
library(forecast)
ig.b <- window(ig, end=c(2018, 6))
reg.ig.b <- arima(ig.b, order=c(1,0,1))

prev1 <- forecast(reg.ig.b, h=1, level=c(25,50,90))
summary(prev1)


plot.forecast(prev1, include=24, 
              main='Previsão do IGP-M com um modelo ARMA(1,1)')

#Previsão para agosto
reg.ig.a <- arima(ig, order=c(0,0,4))
prev2 <- forecast(reg.ig.a, h=1,level=c(30,60))
summary(prev2)

plot.forecast(prev2, include=24, 
              main='Previsão do IGP-M para dezembro com um modelo MA(4)')
