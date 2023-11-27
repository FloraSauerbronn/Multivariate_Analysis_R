## PPG Oceanografia UFSC 
# Analise Multivariada de Dados OceanogrAficos - 2023

# Profa. Carla Bonetti


# PRATICA: ANOVA & MANOVA  & PERMANOVA

# PACOTES
install.packages("nortest")     # testes de normalidade
install.packages("MVN")         # testes de normalidade multivariada
install.packages("dplyr")       # Medidas estatistiscas por grupo (estatistica descritiva) - pipe
install.packages("ggplot2")     # representacoes graficas
install.packages("ggpubr")      # representacoes graficas
install.packages("car")         # testes de Levene e ANOVA  # se necessario instalar tambem pacote FSA
install.packages("FSA")         # Teste de Dunn (multipla comparacao a posteriori)
install.packages("rstatix")     # Teste Box?s M (homocedastia multivariada)

setwd("C:/Users/flora/OneDrive/Documentos/MESTRADO_UFSC/rotinas/R")

# Ler dados CSV (MS DOS) - criando uma tabela de dados - dataframe (df)
df <-read.csv2("Tab_Manova.csv", stringsAsFactors = TRUE)  # le celulas de texto como "fator"

#Se necessario caso precise corrigir as variaveis
df$FATOR<-as.factor(df$FATOR)
df$DEP<-as.numeric(df$DEP)

#Observando dataframe
summary(df)

#### ANALISE DESCRITIVA-GRAFICA DOS GRUPOS E TESTES DE NORMALIDADE ####

library(MVN)

mvn(data = df[,7:8], mvnTest = "royston", univariateTest = "SW", univariatePlot = "histogram")  

## Se necessario - Transformacao Box Cox (se necessario)

BoxCox.lambda(df$DEP, lower=-5, upper=5)                   # method = c("loglik")
df$DEP_t <- (df$DEP^(1.45) - 1)/1.45          # (var^lambda - 1/lambda)


# DESCRICAO ESTATISTICA E GRAFICA DOS GRUPOS

library(dplyr)                   # Medidas estatistiscas por grupo (estatistica descritiva)
library(ggplot2)                 # representacoes graficas

# obter valores de medida central e dispers?o por grupo (fator)
group_by(df, FATOR) %>%                #(df, var categorica)
  summarise(                          #fun??o do dplyr
    count = n(),
    mean = mean(DEP, na.rm = TRUE),    #na.rm ignora celulas NA
    sd = sd(DEP, na.rm = TRUE),
    median = median(DEP, na.rm = TRUE),
    quantile1 = quantile(DEP, probs = 0.25,na.rm = TRUE),
    quantile3 = quantile(DEP, probs = 0.75,na.rm = TRUE),
    IQR = IQR(DEP, na.rm = TRUE))

# Histogramas por grupos (densidade de probabilidade)
ggplot(df, aes(x = DEP, fill = FATOR)) + geom_density(alpha = 0.4)+
  scale_fill_manual(values = c("lightblue","darkgray","cyan4"))+
  labs(x = "FATOR",
       y = "Densidade - DEP") +
  theme_classic()

# Diagrama de caixa por grupo (1o. e 3o. quartis; 1,5*IQR)
ggplot(df, aes(x=df$FATOR, y=df$DEP)) +
  geom_boxplot(fill='seagreen', color="black")+
  labs(x = "FATOR",
       y = "DEP")+
  stat_summary(fun=mean, geom="point", shape=20, size=4) +
  theme_classic()


## TESTES DE NORMALIDADE UNIVARIADOS POR GRUPOS  

library(nortest)   
# Qdo se tem muitos niveis (ANALISE FATORIAL), esse procedimento pode ser substituido 
# pelo teste da normalidade dos residuos da ANOVA (indicado ao final da ANOVA)
# Teste de Normalidade a aplicada a cada grupo de amostras - mais importante para ANOVA
with(df, shapiro.test(DEP1[FATOR1 == "FA"]))
with(df, shapiro.test(VAR[FATOR == "FB"]))
with(df, shapiro.test(VAR[FATOR == "FC"]))


# TESTAR HOMOCEDASTICIDADE DAS VARIANCIAS
library(car)

leveneTest(DEP ~ FATOR, data = df, center = mean)    

#OU
bartlett.test(DEP ~ FATOR, data = df)


####   ANOVA   #### 

# ANALISE DE VARIANCIA FATOR UNICO - ONE-WAY (UM FATOR)
res.aov <- aov(DEP ~ FATOR, data = df)

summary(res.aov)

# Se nao houver homocedasticidade, pode-se usar "Welch one-way" test" como alternativa a ANOVA 
oneway.test(DEP ~ FATOR, data = df, var.equal=FALSE)

## POST HOC - Teste de TUKEY 
TukeyHSD(res.aov)

# ou 
pairwise.t.test(df$DEP, df$FATOR, p.adj = "bonf")           # Bonferroni



### RESULTADOS GRAFICOS DA ANOVA - DIAGRAMAS DE CAIXA

library(ggpubr)

# diagrama de caixa com resultado do teste de Tukey
tukey.test <- aov(df$DEP ~ FATOR, data = df) %>%
  tukey_hsd()

ggboxplot(df, x = "FATOR", y = "DEP",                                              
          color = "black", fill="seagreen",                                   # order = c("NIVEL1", "NIVEL2", "NIVEL3"),
          ylab = "DEP (unidade)", xlab = "FATOR") + 
          stat_summary(fun=mean, geom="point", shape=20, size=3)+ 
          stat_compare_means(label.y = 0,                                     # altura no eixo Y onde sera apresentado o resultado da ANOVA
                               hide.ns = FALSE,                               # apresenta o simbolo ns (n?o significativo)
                               method = "anova")   +                          # t.test, wilcox.test, anova, kruskal.test
          stat_pvalue_manual(tukey.test, label = "p.adj.signif", 
                               y.position = c(20, 22, 24))                    # altura no eixo Y onde serao apresentados os valores de p; o numero de valores deve ser igual ao numero de pares de compara??o

# ns: p > 0.05 *   p <= 0.05 **   p <= 0.01  ***   p <= 0.001  ****  p <= 0.0001


## ANOVA UNIFATORIAL NAO PARAMETRICA (TESTE de KRUSKAL WALLIS) ##

kruskal.test(DEP ~ FATOR, data = df)

# POST HOC test: Dunn test for multiple comparisons of groups
library(FSA)

dunn <- dunnTest(DEP ~ FATOR, data=df, method="bonferroni")

# ouTeste de Comparacoes Multiplas 

pairwise.wilcox.test(df$DEP, df$FATOR, paired = FALSE,                   # PAIRED = FALSE calcula Mann-Whitney 
                     p.adjust.method = "BH")                             #correcao de Bonferroni


# diagrama de caixa

stat.test <- compare_means(DEP ~ FATOR, data = df,
  method = "wilcox.test")                                               # default wilcox.test = n?o pareado (Mann-Whitney)

ggboxplot(df, x = "FATOR", y = "DEP",                                               
          color = "black", fill="seagreen",
  ylab = "DEP (unidade)", xlab = "FATOR") + 
  stat_summary(fun=mean, geom="point", shape=20, size=3)+ 
  stat_compare_means(label.y = 0,                                               # altura no eixo Y onde sera apresentado o resultado da ANOVA
                     hide.ns = FALSE,                                           # apresenta o simbolo ns (n?o significativo)
                     method = "kruskal.test")   +                                         # t.test, wilcox.test, anova, kruskal.test
  stat_pvalue_manual(stat.test, label = "p.signif", 
                     y.position = c(20, 22, 24))                   # altura no eixo Y onde serao apresentados os valores de p; o numero de valores deve ser igual ao numero de pares de compara??o


## ANALISE DE VARIANCIA FATORIAL ##

# Modelo aditivo (assume-se que n?o h? intera??o entre os dois fatores - ou seja, s?o fatores independentes)
# Se h? intera??o entre os fatores (efeito sinerg?tico), subsituir + por * (DEP ~ fator 1 * fator2)

res.anova2 <- aov(DEP ~ FATOR1 * FATOR2, data = df)
summary(res.anova2)

# Teste de TUKEY para comparacoes multiplas
TukeyHSD(res.anova2, which = "FATOR1")
TukeyHSD(res.anova2, which = "FATOR2")

#ANOVA FATORIAL (DESIGN NAO BALANCEADO) 
#SE O NUMERO DE CASOS NAO ? IGUAL ENTRE OS FATORES, USE A FORMULACAO ABAIXO (UNBALANCED DESIGN)
res.anova_unbal <- aov(DEP ~ FATOR1 + FATOR2, data = df)
Anova(res.anova_unbal, type = "III")


#####  ANALISES MULTIVARIADAS   #####

### MANOVA ###

# verificar normalidade e homocedasticidade das covariancias (MVN e Box's M-test)
library(rstatix)

mvn(data = df[,c("DEP1", "DEP2")], mvnTest = "royston", univariateTest = "SW", univariatePlot = "histogram")  
box_m(df[, c("DEP1", "DEP2")], df$FATOR)


# Executar a MANOVA

res.man <- manova(cbind(df$DEP1,df$DEP2) ~ FATOR, data = df)

summary(res.man,test="Wilks")    # se o teste de Wilks nao for especificado, o default  ? o teste de Pillai


### PERMANOVA ###

#install.packages('vegan')       # PERMANOVA (ordination, diversity & dissimilarities)
library(vegan)

# padronizar as variaveis 
df$DEP1_z <- scale(df$VDEP1)
df$DEP2_z <- scale(df$VDEP2)

tab <- df[, c("DEP1", "DEP2")]   ## selecionar as var dependentes

permanova <- adonis2(tab ~ FATOR, data=df, permutations = 999, method="euclidean")  #roda a Permanova

permanova

# representacao grafica por NMDS  (# https://analises-ecologicas.com/cap9#permanova)


#Teste Post Hoc (comparacoes multiplas)

#install.packages("pairwiseAdonis")
library(pairwiseAdonis)

pairwise.adonis(tab, df$FATOR, sim.method = "euclidean",
                p.adjust.m = "bonferroni")



### EXTRA 
### EM CASOS ESPECIAIS: ANALISE DOS RESIDUOS  DA ANOVA

# Teste de normalidade dos residuos
anova1_residuals <- residuals(object = res.aov)    # extrai os residuos
shapiro.test(x = anova1_residuals )  

# Teste de homocedasticidade dos residuos
leveneTest(anova1_residuals ~ SAL, data = df, center = mean)    # ou bartlett.test(anova1_residuals~df$AREA)

# Teste de normalidade dos residuos - Anova Fatorial
anova2_residuals <- residuals(object = res.anova2)
shapiro.test(x = anova2_residuals )

# Teste de homocedasticidade dos residuos - Anova Fatorial
leveneTest(anova2_residuals ~ AREA, data = df, center = mean)    # ou bartlett.test(anova1_residuals~df$AREA)

# TESTE DE NORMALIDADE DOS RESIDUOS - MANOVA
manova_residuals <- residuals(object = res.man)    # extrai os residuos
shapiro.test(x = manova_residuals ) 


#### COMPLEMENTOS  ####

# Remover notacao cientifica
options(scipen=999, digits = 4)

#Exportando dados em CSV
library(readr)
write_delim(df, "teste123.csv", delim = ";", na = "-")


# salvar todo o ambiente de trabalho
save.image() #ou
save.image(file = "Anova_Manova.RData")

#Carregando arquivos salvos
load("teste.RData")

#REFERENCIAS
#https://www.statology.org/multivariate-normality-test-r/
#http://www.sthda.com/english/wiki/one-way-anova-test-in-r
#http://www.sthda.com/english/wiki/two-way-anova-test-in-r
#https://www.datanovia.com/en/lessons/one-way-manova-in-r/
#http://www.sthda.com/english/wiki/wiki.php?title=manova-test-in-r-multivariate-analysis-of-variance
#https://rpubs.com/collnell/manova
#https://f-santos.gitlab.io/2020-05-07-npmanova.html
#https://rpubs.com/vitorks/665678