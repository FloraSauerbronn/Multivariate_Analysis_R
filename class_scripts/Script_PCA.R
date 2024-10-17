## PPG Oceanografia UFSC 
# Analise Multivariada de Dados Oceanograficos - 2023
# Profa. Carla Bonetti

#PRATICA: PCA


install.packages(c("ggbiplot"))     #ERRO# graficos biplot
install.packages("scatterplot3d")
install.packages("ggpmisc")         #plotar informações extra no grafico de dispersão"
install.packages(c("factoextra"))
install.packages(c("FactoMineR"))
install.packages ("GGally")     # função ggpair (matriz grafica de correlação)

library(devtools)
library(ggplot2)                
library(ggpubr)                 
library(dplyr) 
library(MVN)         #testes de normalidade
library(corrplot)
library(ggbiplot)
library(scatterplot3d)
library(ggpmisc)
library(FactoMineR)
library(factoextra)
library(GGally)


setwd("C:/Users/dell/OneDrive/Desktop/OCN_PPG 2023_Analise Multivariada/Praticas")

df <-read.csv2("Tab_PCA_MDS.csv")

summary(df)

mvn(data = df[,2:7], univariateTest = "SW", mvnTest = "royston", univariatePlot = "scatter")

##CORRELAÇÃO MULTIPLA

# Matriz  de correlação (varios diagramas de dispersão ao mesmo tempo)
matriz_c <- round(cor(df[,2:7], method = 'pearson'),2) # OU method = 'spearman'  #round = arredonda o numero de decimais

matriz_c 

# matriz grafica com dispersao, histograma e valores de r e significancia p
corrplot(cor(df[,2:7,], method ='pearson'), type = "lower",
         method = c("number"), tl.col = "black", number.digits = 1)

#OU
library(PerformanceAnalytics)

chart.Correlation(df[,2:7], histogram=TRUE, pch="+")   #method = c("pearson", "kendall", "spearman") - pearson é o default


## PCA - Analise de Componentes Principais

# a função PCA neste pacote automatica padroniza os valores de entrada ( = matriz de correlacao)

res.pca <- PCA(df[, 2:7], graph = FALSE,    # gera os resultados que serão usados nos gráficos
               scale.unit = TRUE,           #scale.unit: se tRUE os dados são padronizados previamente (matriz de correlação)
               ncp = 6)                     #ncp: numero de dimensões mantidas no resultado final

res.pca
# Eigenvalues
round(get_eigenvalue(res.pca),3)  # Extract the eigenvalues/variances of principal components

#Screeplot
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 100)) # ylim é amplitude do eixo y (percentual de explicação)

# Carga das variáveis (loadings)

VAR <-get_pca_var(res.pca)   
VAR[["cor"]]
VAR[["cos2"]]
VAR[["contrib"]]              # expressa em porcentagem

#plot das variáveis x PCs  
fviz_pca_var(res.pca, col.var = "contrib", axes = c(1, 2),    # cos2 é uma medida da qualidade de representação das variáveis no espaço 2d (plot),       
             gradient.cols = c("gray", "blue", "red"), 
             repel = TRUE ) +  # Avoid text overlapping
              xlim(-1, 1) + ylim (-1, 1)

#Visualizar  também PC1xPC3 (axes = c(1,3))


# Coordenadas dos casos (scores)
CASO <- get_pca_ind(res.pca)    
CASO[["coord"]]

#plot dos casos x PCs 
fviz_pca_ind(res.pca, axes = c(1, 2), col.ind = "contrib", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),   
             pointshape = 20, pointsize = 2, labelsize = 4,
             repel = FALSE,
             legend.title = "contrib")  +
             xlim(-4, 4) + ylim (-4, 4)


#identificar casos por gradiente ou por grupos
fviz_pca_ind(res.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = df[,1], # coluna de referencia para classes ou gradiente de cores 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), #ou igual ao numero de grupos (classes)
             addEllipses = FALSE, # Concentration ellipses
             legend.title = "ID")+
             xlim(-4, 4) + ylim (-4, 4)


# BIPLOT - plotar var e casos por PCs

fviz_pca_biplot(res.pca, col.var = "contrib", 
                pointshape = 20, pointsize = 1, labelsize = 4,
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                repel = TRUE ) +
                xlim(-4, 4) + ylim (-4, 4)
#OU
fviz_pca_biplot(res.pca, col.ind = "contrib", 
                pointshape = 20, pointsize = 2,
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                repel = TRUE ) +
                xlim(-4, 4) + ylim (-4, 4)



# salvar todo o ambiente de trabalho
save.image() #ou
save.image(file = "Pratica_PCA.RData")

#Carregando arquivos salvos
load("teste.RData")

#REFERENCIAS
# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials
# http://www.sthda.com/english/wiki/correlation-test-between-two-variables-in-r
# http://www.sthda.com/english/articles/40-regression-analysis/165-linear-regression-essentials-in-r
# http://www.sthda.com/english/articles/39-regression-model-diagnostics/161-linear-regression-assumptions-and-diagnostics-in-r-essentials
# https://rstudio-pubs-static.s3.amazonaws.com/251265_9f65a6ca663f4eff814ae8b1f5a04ffe.html
# http://www.sthda.com/english/wiki/fviz-pca-quick-principal-component-analysis-data-visualization-r-software-and-data-mining