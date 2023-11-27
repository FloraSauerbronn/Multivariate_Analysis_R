###############################################################################
#SCRIPT REFERENTE A TAREFA FINAL DA DISCIPLINA DE POSGRADUACAO EM OCEANOGRAFIA#
#OCE410006 - Analise Multivariada de dados Oceanograficos                     #
#Professora Dra Carla Bonetti                                                 #
#Aluna Flora Medeiros Sauerbronn                                              # 
#Matricula:202302384                                                          #
###############################################################################

##### Carregando os Dados #####

##Lendo o Arquivo csv como data frame

setwd("C:/Users/flora/OneDrive/Documentos/MESTRADO_UFSC/rotinas/R/Dados")
df <-read.csv2("RS-4_2019-02-19_2021-09-08_process.csv", stringsAsFactors = TRUE) 

#Notamos que os dados provenientes do ADCP estão com muitos NaNs, prejudicando 
#as análises, vamos removelos

#Retirada de campos ADCP
df <- df[, !colnames(df) %in% "Avg_W_Tmp2"]
df <- df[, !colnames(df) %in% "Avg_Cel1_Mag"]
df <- df[, !colnames(df) %in% "Avg_Cel1_Dir"]
df <- df[, !colnames(df) %in% "Avg_Cel1_Dir_N"]

#Retirar os indices do python 
df <- df[, !colnames(df) %in% "X"]

#Retirar declividade magnéica pois é negativa e não serve para o box-cox
df <- df[, !colnames(df) %in% "M_Decl"]


#Acertando os campos de datas como Date



## Explorando os dados 
summary (df)

########################### VERIFICAÇÃO NORMALIDADE  ###########################



## Q-Q plot

## Altura de Onda
qqnorm(df$Hsig, pch = 20, frame = FALSE, main = "Q-Q plot Altura Significativa de Onda")
qqline(df$Hsig, col = "steelblue", lwd = 2)

# Temperatura Média, utilizar a 1 de CTD pois a 2 de ADCP tem muitos dados nulos
qqnorm(df$Avg_W_Tmp1, pch = 20, frame = FALSE, main = "Q-Q plot Temperatura Média CTD (Cº)")
qqline(df$Avg_W_Tmp1, col = "steelblue", lwd = 2)


## Salinidade Média
qqnorm(df$Avg_Sal, pch = 20, frame = FALSE, main = "Q-Q plot Salinidade Média (psu)")
qqline(df$Avg_Sal, col = "steelblue", lwd = 2)

## Clorofila Média
qqnorm(df$Avg_Chl, pch = 20, frame = FALSE, main = "Q-Q plot Clorofila Média (ug/L)")
qqline(df$Avg_Chl, col = "steelblue", lwd = 2)

## Clorofila Média
qqnorm(df$Avg_DO, pch = 20, frame = FALSE, main = "Q-Q plot Oxigênio Dissolvido Média (ml/L)")
qqline(df$Avg_DO, col = "steelblue", lwd = 2)

##(ANOTAÇÕES)
##Após a série de Plots observamos que importantes variáveis não apresentavam uma distribuição normal
##Vamos verificar utilizando um teste diferente

##Shapiro-Wilk
#Vamos testar com ele por ser um teste mais versátil para diferentes tipos de 
#tamanho de amostras e comparar com o Q-Q plot

##Instalando o pacote necessário
# install.packages("MVN")

# Carregue o pacote MVN
library(MVN)

#Aplicando o teste Shapiro
#teste de normalidade multivariado - royston
#teste de normalidade univariada - Shapiro - Wilk
#Pegamos apenas as colunas de 5 à 15 que são nossas variáveis ambientais
mvn(data = df[,5:13], mvnTest = "royston", univariateTest = "SW", univariatePlot = "histogram")  

##(ANOTAÇÕES)
##Realmente os dados não tem uma dstribuição normal
##Vamos Realizar uma padronização


###########################  TRANSFORMACAO  ####################################

#vamos realizara  transformação para normalizar os dados e assim aplicarmos o PCA
#e o K-means.
#Escolhemos esse tipo de transformação pois
##   Transformacao BOX COX 

#install.packages("forecast")   #BoxCox LAMBDA
#install.packages("MASS")
library(forecast)
library(MASS)


#Trnsformando as variáveis e inserindo na matriz


#uma por uma infelizmente ...
# Calcular valor de lambda 1  ALTURA ONDA
Lambda_var <- BoxCox.lambda(df$Hsig,lower=-5,upper=5)
#Transforma o dado
df$Hsig_t <- (df$Hsig^(Lambda_var) - 1)/Lambda_var   # (var^lambda - 1/lambda)

# Calcular valor de lambda  2  TEMPERATURA
Lambda_var <- BoxCox.lambda(df$Avg_W_Tmp1,lower=-5,upper=5)
#Transforma o dado
df$Avg_W_Tmp1_t <- (df$Avg_W_Tmp1^(Lambda_var) - 1)/Lambda_var   # (var^lambda - 1/lambda)

# Calcular valor de lambda  3  SALINIDADE
Lambda_var <- BoxCox.lambda(df$Avg_Sal,lower=-5,upper=5)
#Transforma o dado
df$Avg_Sal_t <- (df$Avg_Sal^(Lambda_var) - 1)/Lambda_var   # (var^lambda - 1/lambda)

# Calcular valor de lambda  4 WAVE DIRECTION
Lambda_var <- BoxCox.lambda(df$Avg_Wv_Dir,lower=-5,upper=5)
#Transforma o dado
df$Avg_Wv_Dir_t <- (df$Avg_Wv_Dir^(Lambda_var) - 1)/Lambda_var   # (var^lambda - 1/lambda)

# Calcular valor de lambda  5 WAVE DIRECTION NORTH MAGNETIC
Lambda_var <- BoxCox.lambda(df$Avg_Wv_Dir_N,lower=-5,upper=5)
#Transforma o dado
df$Avg_Wv_Dir_N_t <- (df$Avg_Wv_Dir_N^(Lambda_var) - 1)/Lambda_var   # (var^lambda - 1/lambda)

# Calcular valor de lambda  6 TURBIDEZ
Lambda_var <- BoxCox.lambda(df$Avg_Turb,lower=-5,upper=5)
#Transforma o dado
df$Avg_Turb_t <- (df$Avg_Turb^(Lambda_var) - 1)/Lambda_var   # (var^lambda - 1/lambda)

# Calcular valor de lambda  7 CLOROFILA
Lambda_var <- BoxCox.lambda(df$Avg_Chl,lower=-5,upper=5)
#Transforma o dado
df$Avg_Chl_t <- (df$Avg_Chl^(Lambda_var) - 1)/Lambda_var   # (var^lambda - 1/lambda)

# Calcular valor de lambda  8 OXIGENIO
Lambda_var <- BoxCox.lambda(df$Avg_DO,lower=-5,upper=5)
#Transforma o dado
df$Avg_DO_t <- (df$Avg_DO^(Lambda_var) - 1)/Lambda_var   # (var^lambda - 1/lambda)

# Calcular valor de lambda  9  MÉDIA ALTURA ONDA
Lambda_var <- BoxCox.lambda(df$HAvg,lower=-5,upper=5)
#Transforma o dado
df$HAvg_t <- (df$HAvg^(Lambda_var) - 1)/Lambda_var   # (var^lambda - 1/lambda)

#agora salvando as transformações numa matriz específica 
#Adicionando as colunas de datas
df_t <- df[,1:4]  
#Adicionando as colunas transformadas
df_t <- cbind(df_t, df[, 14:22])

#Vamos verificae a normalidade DE NOVO
#Pegamos apenas as colunas de 5 à 15 que são nossas variáveis ambientais
mvn(data = df_t[,5:13], mvnTest = "royston", univariateTest = "SW", univariatePlot = "histogram")  

#Vamos verificar a diferença antes e depois da transformação com uma variável
#Usando o Q-Q Plot
## Clorofila Média
qqnorm(df$Avg_Chl, pch = 20, frame = FALSE, main = "NÃO NORMAL Clorofila Média (ug/L)")
qqline(df$Avg_Chl, col = "steelblue", lwd = 2)

qqnorm(df$Avg_Chl_t, pch = 20, frame = FALSE, main = "NORMAL Clorofila Média (ug/L)")
qqline(df$Avg_Chl_t, col = "steelblue", lwd = 2)


#Com nossos reultados melhorados vamos partir para os PCAs

####################################  PCA  #####################################

#Instalação Pacotes
#install.packages(c("ggbiplot"))     #ERRO# graficos biplot
#install.packages("scatterplot3d")
#install.packages("ggpmisc")         #plotar informa??es extra no grafico de dispers?o"
#install.packages(c("factoextra"))
#install.packages(c("FactoMineR"))
#install.packages ("GGally")  
#install.packages("devtools")
library(devtools)
library(ggplot2)                
library(ggpubr)                 
library(dplyr) 
library(corrplot)
library(ggbiplot)
library(scatterplot3d)
library(ggpmisc)
library(FactoMineR)
library(factoextra)
library(GGally)

# Matriz  de correlação (varios diagramas de disperssão ao mesmo tempo)
#Verificação da correlação de pearson
matriz_c <- round(cor(df_t[,5:13], method = 'pearson'),2) # OU method = 'spearman'  #round = arredonda o numero de decimais


install.packages("FactoMineR")
install.packages("factoextra")

# Carregue o pacote factoextra
library(factoextra)


library(FactoMineR)

## PCA - Analise de Componentes Principais
#Visualização da relação das variáveis entre sí
#Ordenação de Importancia

#A função PCA neste pacote automatica padroniza os valores de entrada ( = matriz de correlacao)

res.pca <- PCA(df_t[,5:13], graph = FALSE,    # gera os resultados que ser?o usados nos gr?ficos
               scale.unit = TRUE,           #scale.unit: se tRUE os dados s?o padronizados previamente (matriz de correla??o)
               ncp = 6)                     #ncp: numero de dimens?es mantidas no resultado final

res.pca
# Eigenvalues
round(get_eigenvalue(res.pca),3)  # Extract the eigenvalues/variances of principal components

#Screeplot Indix a PORCENTAGEM DE VARIÃNCIA DE CADA VARIÁVEL
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 100)) # ylim ? amplitude do eixo y (percentual de explica??o)

# Carga das vari?veis (loadings)

VAR <-get_pca_var(res.pca)   
VAR[["cor"]]
VAR[["cos2"]]
VAR[["contrib"]]              # expressa em porcentagem


# COMO AS VARIÁVEIS SE RELACIONAM ENTRE SÍ

#plot das vari?veis x PCs  
fviz_pca_var(res.pca, col.var = "contrib", axes = c(1, 2),    # cos2 ? uma medida da qualidade de representa??o das vari?veis no espa?o 2d (plot),       
             gradient.cols = c("gray", "blue", "red"), 
             repel = TRUE ) +  # Avoid text overlapping
  xlim(-1, 1) + ylim (-1, 1)

#Visualizar  tamb?m PC1xPC3 (axes = c(1,3))


#################################  AGRUPAMENTO  ################################

install.packages("cluster")   
install.packages("NbClust")
install.packages("gplots")   # heatmap (cluster biplot)
install.packages("ggdendro")  #dendrogramas
install.packages(ggbiplot)
library(dplyr)
library(ggplot2)
library(ggpubr) 
library(ggbiplot)
library(GGally)         # scatter plot
library(MVN)            # testes de normalidade
library(FactoMineR)
library(factoextra)     # analise de agrupamento (funcoes fviz...)
library(vegan)          # matrizes de associacao
library(cluster)        # analise de agrupamento 
library(ggdendro)       # dendrograma
library(gplots)         # heatmap (cluster biplot)
library(NbClust)        # definicao do numero de clusters


##ANALISE DE AGRUPAMENTO POR "KMEANS" (origem arbitraria)

#Retira as datas e deixa apenas os numéricos
 dz <-df_t[,2:13]


# Definir numero otimo de grupos

# Silhouette method
 ####Indica um possível bom fator que é 4 clusters das estações do ano
 
fviz_nbclust(dz, kmeans, method = "silhouette" )+         # method = "gap_stat"
  labs(subtitle = "Silhouette method")

# Gap statistic
fviz_nbclust(dz, kmeans, method = "gap_stat")+         # method = "gap_stat"
  labs(subtitle = "Gap statistic")

# Elbow
fviz_nbclust(dz, kmeans, method = "wss")+         # wss = elbow method (within sum of square)
  labs(subtitle = "Elbow method")

