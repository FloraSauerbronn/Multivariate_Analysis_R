## PPGOceano - ANALISE MULTIVARIADA DE DADOS OCEANOGR?FICOS 
# Profa. Carla Bonetti

### PRATICA: DISTRIBUI??ES E TESTES DE NORMALIDADE  ###


setwd("C:/Users/flora/OneDrive/Documentos/MESTRADO_UFSC/rotinas/R")


# Ler dados CSV (MS DOS) 
df <-read.csv2("Tab_TestesNormalidade.csv")  

# Explore os dados     
summary (df)

### HISTOGRAMAS ###

hist(df$N_1, freq = TRUE,                     #se freq=FALSE retorna um histograma de densidade
     breaks = 7,
     main = "Histograma", 
     xlab = "classes (N_1)", ylab = "Freq. Absoluta", 
     col = c("seagreen"), 
     border = c("gray"),
     right = TRUE,
     axes = TRUE) 

# opcional (para histogramas de densidade)
lines(density(df$N_1, na.rm = TRUE ))        # apenas para histogramas de densidade (freq = FALSE)

## Histograma mais customizado

# numero de classes (k) e intervalo de classe (h)
# k = numero de classes (opcional  >>  usar formula de Sturges k =  1 + 3.3 * log10(n))
# h =  amplitude/k  sendo 

k <- 1+3.3*log10(50)        # arredondar para numero inteiro

min <- min(df$N_1)
max <- max(df$N_1)

h <- (max-min)/k            # pode ser arredondado
h <- 0.5

limites <- seq(min,max,h)   # a fun??o seq(), informa o local das quebras de classe e gera sequ?ncias regulares
#OU
limites <- seq(3.5,7,0.5)   # valores arredondados

Xm <- (limites+(h/2))       # Ponto m?dio

hist(df$N_1,  
     main = "Distribuição de Frequencias", 
     xlab = "ponto médio das classes (N_1)", ylab = "Freq. Absoluta", 
     col = c("seagreen"), 
     border = c("gray"),
     breaks = limites,
     right = TRUE,
     xlim = c(3,7), ylim = c(0,15),
     axes = FALSE)
axis(1,Xm,Xm)   # O primeiro argumento (1) referencia o eixo das abscissas (x); Xm, Xm indica os locais do ticks e o vetor dos valores (pontos m?dios)
axis(2, at=seq(0,15,5))   # O primeiro argumento (2) referencia o eixo das ordenadas (y); seq(min, max e intervalo)

# OPCIONAL (substituir express?o abaixo na linha de defini??o da escala do eixo X)
axis(1, at=seq(3.5,7,0.5))


# GR?FICO QQ-Plot (quantil-quantil plot)
qqnorm(df$N_1, pch = 20, frame = FALSE, main = "Q-Q plot N_1")
qqline(df$N_1, col = "steelblue", lwd = 2)


### TESTES DE NORMALIDADE UNIVARIADA

#install.packages("nortest")    # usado para os testes de Normalidade
library(nortest)    

# Teste de Normalidade de Shapiro-Wilk (W)
shapiro.test(df$N_1)

# Teste de Normalidade de Lilliefors (Kolmogorov-Smirnov) 
lillie.test(df$N_1)

# Teste de Normalidade de Anderson-Darling 
ad.test(df$N_1)


### TESTE DE NORMALIDADE MULTIVARIADA 

#install.packages("MVN")     
library(MVN)                 # testes de normalidade uni e multivariados - plot de graficos 

# testes univariados:  Shapiro-Wilk ("SW"), Lilliefors ("Lillie"), Anderson-Darling ("AD"), outros
# testes multivariados: "mardia" for Mardia's test, "hz" for Henze-Zirkler's test, "royston" for Royston's test, outros

mvn(data = df[,8:9], mvnTest = "royston", univariateTest = "SW", univariatePlot = "histogram")     # OU univariatePlot = "qqplot" , "box"

# em casos particulares
mvn(data = df[,2:9], univariatePlot = "box")


### TRANSFORMACOES ###

# logaritmica -  log2(VAR+1)  opcoes log - retorna logaritmo natural (2,718); log2; log10
# raiz quadrada -  sqrt(VAR+1)
# exponencial - (VAR^2) OU exp(VAR)  - calcula e^x (sendo e = 2,718)
# inversa - (1/VAR)

# Transforma a variavel original e cria uma nova coluna no dataframe

df$ASpos_t <- log10(df$ASpos+1)


##   Transformacao BOX COX 

#install.packages("forecast")   #BoxCox LAMBDA
library(forecast)

# Calcular valor de lambda
Lambda_ASneg <- BoxCox.lambda(df$ASneg,lower=-5,upper=5)
Lambda_ASpos <- BoxCox.lambda(df$ASpos,lower=-5,upper=5)
Lambda_Lept <- BoxCox.lambda(df$Lept,lower=-5,upper=5)
Lambda_Bimodal <- BoxCox.lambda(df$Bimodal,lower=-5,upper=5)
###Essa aqui eu alterei porque a função da professora não esta rodando ###
#boxcox(df$ASpos ~ 1, lambda = seq(-5,5,0.1))

# Transforma a variavel original e cria uma nova coluna no dataframe

df$ASneg_t <- (df$ASneg^(2.5) - 1)/2.5     # (var^lambda - 1/lambda)
df$ASpos_t <- (df$ASpos^(-0.27) - 1)/-0.27
df$Lept_t <- (df$Lept^(1.76) - 1)/1.76
df$Bimodal_t <- (df$Bimodal^(0.84) - 1)/0.84

mvn(data = df[,10:13], mvnTest = "royston", univariateTest = "SW", univariatePlot = "histogram")     # OU univariatePlot = "qqplot" , "box"

# Opcional
#install.packages("MASS")       # Boxcox GRAFICO
library(MASS)
Box = boxcox(df$ASpos ~ 1, lambda = seq(-3,3,0.1))


### OUTRAS TRANSFORMACOES (SEM AJUSTE DE NORMALIDADE)  ###

# Centralizar  (nova media = zero)
df$N_1_cent <- (df$N_1 - mean(df$N_1))

summary(df$N_1_cent)


# Re-escalar (nova distribuicao entre 0 e 1)

df$N_1_esc <- ((df$N_1 - min(df$N_1)))/(max(df$N_1)-min(df$N_1))

summary(df$N_1_esc)


# Padronizar (valor de z - adimensional; media zero e desvio padrao 1)

df$ASpos_z <- scale(df$ASpos)

mean(df$ASpos_z)
sd(df$ASpos_z)


### EXTRAS

# numero de matrizes/graficos por tela
par(mfrow = c(1, 1))                #apenas 1
par(mfrow = c(1, 2))                # 2 por linha

# SEM NOTACAO CIENTIFICA
options(scipen=999)
print(resultados, digits = 2)

# salvar todo o ambiente de trabalho
save.image() #ou
save.image(file = "Pratica_Normalidade.RData")


#REFERENCIAS

# Korkmaz et al 2015_MVN a package for assesing multivariate normality
# http://www.sthda.com/english/wiki/normality-test-in-r
# http://www.sthda.com/english/wiki/qq-plots-quantile-quantile-plots-r-base-graphs
# https://cran.r-project.org/web/packages/MVN/vignettes/MVN.html#14_Royston%E2%80%99s_MVN_test
# https://rpubs.com/melinatarituba/356739
# https://rdrr.io/cran/MASS/man/boxcox.html