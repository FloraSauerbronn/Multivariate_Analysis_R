## PPG Oceanografia UFSC 
# Analise Multivariada de Dados OceanogrAficos - 2023
# Profa. Carla Bonetti


#### PRATICA: REGRESSAO LINEAR - MODELOS LINEARES ORDINARIOS ####

#install.packages("nortest")     # testes de normalidade
#install.packages("MVN")         # testes de normalidade multivariada
#install.packages("dplyr")       # Medidas estatistiscas por grupo (estatistica descritiva) - pipe
#install.packages("ggplot2")     # representacoes graficas
#install.packages("ggpubr")      # representacoes graficas
#install.packages("broom")       # gera data frame com os resultados da regressao
#install.packages("corrplot")    # matrizes graficas de correlacao
#install.packages("PerformanceAnalytics") 
#install.packages("MASS")        # modelos
#install.packages("ggfortify")   # gera graficos de residuos

library(MVN)
library(ggplot2)
library(ggpubr)
library(nortest)                 
library(broom)
library(ggfortify)
library(corrplot)
library(MASS)
library(PerformanceAnalytics)

setwd("C:/Users/flora/OneDrive/Documentos/MESTRADO_UFSC/rotinas/R")

df <-read.csv2("Tab_Regressao.csv") 

# TESTAR A NORMALIDADE
library(MVN)    
mvn(data = df[,2:5], univariateTest = "SW", univariatePlot = "scatter", mvnTest = "royston")  # ou histogram

# SE NECESSARIO - para ajustar a janela de graficos
par(mfrow = c(1, 1))

#DIAGRAMAS DE DISPERSAO

ggplot(df, aes(x = df$INDEP1, y = df$DEP)) +
  geom_point(color = "black", fill = "seagreen", shape = 21, size = 3, alpha = .5) +   #ref de shapes: https://www.datanovia.com/en/lessons/ggplot-scatter-plot/
  xlab("INDEP1") +          #opcional: + lims(x = c(16, 40), y = c(3, 7)) +
  ylab("DEP") + 
  geom_smooth(method="lm", se=FALSE, color="red")+
  ggpubr::stat_cor(method = "pearson", cor.coef.name = c("R"), p.accuracy = 0.001) +
  theme_classic()

# CORRELACAO DE PEARSON
cor.test(df$INDEP1, df$DEP, method=c("pearson"))  #OPCIONAL ("kendall", "spearman")

# COEFICIENTE DE DETERMINACAO (R2)
cor(df$INDEP1, df$DEP)^2

# MATRIZ DE CORRELACAO
# matriz com valores de r representados com escala de cores ou circulos
 
library (corrplot)   

corrplot(cor(df[,2:5], method ='pearson'), type = "lower", diag = FALSE, method = c("number"),      # OU method = "circle"   
         tl.col = "black", number.digits = 2, tl.srt = 0, tl.cex = 1)                          # tl.srt = 90  angulo de dire??o do titulo das colunas; tl.cex = tamanho do texto

# matriz grafica com dispersao, histograma e valores de r e significancia p
library(PerformanceAnalytics)

chart.Correlation(df[,2:5], histogram=TRUE, pch="+")   #method = c("pearson", "kendall", "spearman") - pearson ? o default


## REGRESSAO LINEAR 

# variavel dependente (eixo Y) e independente (eixo x)

# simples
modelo <- lm(DEP ~ INDEP1, data = df1) # Var dependente ~ var independente

# multiplo
modelo.rm <- lm(DEP ~ INDEP1 + INDEP3 + INDEP3,
                data = df)


# a funcao 'anova' apresenta o modelo de regressao atraves da distribuicao F
anova(modelo.rm)

# resumo do modelo linear e teste de significancia dos coeficientes de regressao
summary(modelo.rm)  

# Interpretando a saida dos dados:
# Estimate (intercep)= coef a; Estimate (VARIAVEL) = coef b
# Std.Error: o erro padr?o das estimativas dos coeficientes. (Isso representa a precis?o dos coeficientes (qto maior o erro padr?o, menor a confian?a na estimativa)
# t value: resultado do teste t (signific?ncia dos coeficientes da equa??o da reta)
# Pr (> | t |): O valor p correspondente ? estat?stica t. 
# Residual standard error (RSE): erro de predi??o ou res?duo (diferen?a m?dia entre observado e previsto).
# Multiple R-squared: coeficiente de Determina??o
# Adjusted R-squared:  coeficiente de Determina??o ajustado ? multiplas vari?veis 
# F-statistic: raz?o F para verifica??o do ajuste da equa??o da reta e valor de p associado a distribui??o F 

###ANALISE DOS RESIDUOS  

ajuste.modelo <- augment(modelo.rm)   #cria um dataframe com as estatisticas da regressao

head(ajuste.modelo)

#significado das metricas 
#.fitted: valores estimados
#.resid: res?duos
#.std.resid: res?duos padronizados(pode ser usado para detectar outliers)
#	extra -- .hat: valores usados para detectar pontos de alta alavancagem (leverage) (influencia de valores extremos)
#	extra -- .cooksd: dist?ncia de Cook, usada para detectar valores influentes (outlier, leverage)


# TESTAR NORMALIDADE DOS RESIDUOS
shapiro.test(ajuste.modelo$.resid)

#Grafico de distribui??o dos residuos padrao
plot(x = ajuste.modelo$INDEP1, 
     y = ajuste.modelo$.std.resid,
     ylim = c(-3, 3),
     pch = 16,
     col = "black",
     xlab = "valores observados - INDEP1", ylab = "res?duos padronizados")
abline(0, 0, col = "darkgray") 
abline(2, 0, col = "gray")
abline(-2, 0, col = "gray")


# DIAGNOSTICO DOS RESIDUOS

par(mfrow = c(1, 1))        # se necessario - ajustar visualizacao de graficos

autoplot(modelo.rm)         # diagnostico geral


autoplot(modelo.rm, 1)      # plota um grafico de residuos x estimados
                                            # pode ser usado para analisar visualmente a homocedasticidade dos residuos

autoplot(modelo.rm,2)       # Normal QQ (avalia a Normalidade)

autoplot(modelo.rm,3)       # localiza??o-escala (res padronizado x val estimados) 
                            # espera-se uma linha reta horizontal - util para verificar heterocedastia

autoplot(modelo.rm,4)       # Distancia de Cook por casos  - # como default mostra os 3 casos mais extremos, para alterar use plot(model, 4, id.n = 5)

autoplot(modelo.rm,5)       # Residuos x Leverage (alavancagem)



# Excluir caso anomalo
library(dplyr)
df1<- slice(df, -100, -99, -95)


## Predi??o (estimar valor de Y usando a equa??o da reta)   y = a + b*x

predict( modelo.rm, data.frame( INDEP1=c(7.0), INDEP2=c(50.0), INDEP3=c(0.04) ) )



### PARTI??O - Treino e Teste dataset

# Amostragem randomica
samplesize <- 0.70*nrow(df)    # 70% dos casos s?o selecioandos aleatoriamente para o dataset treino
set.seed(100)
index <- sample(seq_len(nrow(df)), size = samplesize)
# Definir os dataset
treino = df[index,]
teste = df[-index,]

## Construcao do modelo a partir do dataset treino

mod_treino <- MASS::polr(DEP ~ INDEP1 + INDEP2 + INDEP3,
                         data = datatrain, Hess = T)

summary(mod_treino)

mod_treino <- lm(DEP ~ INDEP1 + INDEP3 + INDEP3,
                data = treino)


#REF
#http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software   
#https://acervolima.com/visualizacao-de-uma-array-de-correlacao-usando-ggplot2-em-r/
#https://bookdown.org/mpfoley1973/statistics/generalized-linear-models-glm.html
#https://stats.oarc.ucla.edu/r/dae/multinomial-logistic-regression/
#https://cran.r-project.org/web/packages/ordinal/vignettes/clm_article.pdf
#https://cehs-research.github.io/eBook_regression/logistic-regression-ex-bronchopulmonary-dysplasia-in-premature-infants.html#logistic-regresion-fit-the-model-to-the-data
# http://www.sthda.com/english/wiki/correlation-test-between-two-variables-in-r
# http://www.sthda.com/english/articles/40-regression-analysis/165-linear-regression-essentials-in-r
# http://www.sthda.com/english/articles/39-regression-model-diagnostics/161-linear-regression-assumptions-and-diagnostics-in-r-essentials
# Bruce, Peter, and Andrew Bruce. 2017. Practical Statistics for Data Scientists. O'Reilly Media.
# https://rpubs.com/leonardoreffatti/627365





