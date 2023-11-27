## PPG Oceanografia UFSC 
# Analise Multivariada de Dados Oceanograficos - 2023
# Profa. Carla Bonetti

#PRATICA: Analise de Agrupamento

install.packages("cluster")   
install.packages("NbClust")
install.packages("gplots")   # heatmap (cluster biplot)
install.packages("ggdendro")  #dendrogramas

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
                        # https://www.rdocumentation.org/packages/NbClust/versions/3.0.1/topics/NbClust



setwd("C:/Users/dell/OneDrive/Desktop/OCN_PPG 2023_Analise Multivariada/Praticas")

df <-read.csv2("Tab_Cluster.csv")

summary(df)

# Testar normalidade
mvn(data = df[,2:7], univariateTest = "SW", mvnTest = "royston", univariatePlot = "histogram")

# qdo necessario padronizar dados 

df.z <- scale(df[, 2:7])    # padronizar tudo


##COEFIENCTES DE ASSOCIACAO

#Matriz de Distancia [including "euclidean", "manhattan", "maximum", "canberra", "binary", "minkowski"]
matriz.dist <- get_dist(df.z, method = "euclidian")   #esta função aceita apenas dados "numeric" 

#Matriz de Dissimilaridade de Bray Curtis (ou  "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "morisita", "mahalanobis"...)
matriz.bc_d <- vegdist(df[, 2:7], method = "bray", binary=FALSE ) 

# Matriz de Correlacao ["pearson", "spearman" or "kendall"]
#distancias baseadas em correlacao  sao definidas pela subtracao 1 - coef. de correlacao 
matriz.cor <- get_dist(df[, 2:7], method = "pearson")    #method = "spearman"


##ANALISE DE AGRUPAMENTO HIERARQUICA

# Estrategias de agrupamento
dendr_euc <- hclust(matriz.dist, method = "ward.D2")   # "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA)
#ou 
dendr_bc <- hclust(matriz.bc_d, method = "average")


# Dendrograma Modo Q

par(mfrow = c(1, 1))

fviz_dend(dendr_euc, cex = 0.5, horiz = FALSE, 
          k = 3, k_colors = "black",                                # k = numero de grupos a serem formados
          rect = TRUE, rect_border = "jco", rect_fill = TRUE,             # rect_border = c("#2E9FDF", "#c066c0", "#E7B800", "#FC4E07")
          main = "",
          xlab = "", ylab = "Distance Eucl (Ward)", sub = "")

#Coeficiente cofenetico
c.cof <- cophenetic(dendr_euc)   # cophenetic(resultado da funcao hclust)
cor(matriz.dist , c.cof)            #cor(matriz de associacao, c.cof) 


# Dendrograma Modo R

#transpor df (MODO R)   
dfz.rmode <- t(df.z)
head(dfz.rmode)

matriz.rmode <- get_dist(dfz.rmode, method = "euclidian")  #criar a matriz

dend_r <- hclust(matriz.rmode, method = "ward.D2")   # "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA)

fviz_dend(dend_r, cex = 0.5, horiz = FALSE, 
          main = "",
          xlab = "", ylab = "Spearman (UPGMA)", sub = "")

# OU matriz de correlacao

matriz.rmode <- get_dist(dfz.rmode, method = "spearman")  #criar a matriz

dend_r <- hclust(matriz.rmode, method = "average")   # "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA)

fviz_dend(dend_r, cex = 0.5, horiz = FALSE, 
          main = "",
          xlab = "", ylab = "Spearman (UPGMA)", sub = "")


### Representação BIPLOT (heatmap)

library("gplots")

heatmap.2(df.z, scale = "none", col = bluered(100), 
          trace = "none", density.info = "none")


##ANALISE DE AGRUPAMENTO POR "KMEANS" (origem arbitraria)

# Definir numero otimo de grupos

# Silhouette method
fviz_nbclust(df.z, kmeans, method = "silhouette" )+         # method = "gap_stat"
  labs(subtitle = "Silhouette method")

# Gap statistic
fviz_nbclust(df.z, kmeans, method = "gap_stat")+         # method = "gap_stat"
  labs(subtitle = "Gap statistic")

# Elbow
fviz_nbclust(df.z, kmeans, method = "wss")+         # wss = elbow method (within sum of square)
  labs(subtitle = "Elbow method")


#OU função NbClust para calcular varios métodos simultâneamente

#data: matrix
#diss: dissimilarity matrix to be used. By default, diss=NULL, but if it is replaced by a dissimilarity matrix, distance should be "NULL"
#distance: the distance measure to be used to compute the dissimilarity matrix. Possible values include "euclidean", "manhattan" or "NULL".
#min.nc, max.nc: minimal and maximal number of clusters, respectively
#method: The cluster analysis method to be used including "ward.D", "ward.D2", "single", "complete", "average", "kmeans" and more.
#To compute NbClust() for kmeans, use method = "kmeans".
#To compute NbClust() for hierarchical clustering, method should be one of c("ward.D", "ward.D2", "single", "complete", "average").

NbClust(df.z, diss = NULL, distance = "euclidean",
        min.nc = 2, max.nc = 6, method = "kmeans")

# Agrupamento

# a partir de dados padronizados
gr.km <- kmeans(df.z, 3, nstart = 25)   # with k = 3; nstar = numero de centroides testados inicialmente, recomendado >25

gr.km$cluster

#Visualizar plot (grupos sobre PCA)
fviz_cluster(gr.km, df.z,
             palette = c("steelblue", "darkgray", "darkgreen", "orange", "tomato", "sienna4" ),  
             geom = c("point", "text"),pointsize = 1.5, labelsize = 12,
             ellipse.type = "convex", ellipse.alpha = 0.2,        # se elipse.type "confidence" definir ellipse.level = 0.95 , 
             main = "",
             ggtheme = theme_bw())



# salvar todo o ambiente de trabalho
save.image() #ou
save.image(file = "AnalAgrupamento.RData")


### EXTRA
# Dendrograma no ggplot

dend <- as.dendrogram(dendr_euc)
dend_gr <- dendro_data(dend, type = "rectangle")

ggplot(dend_gr$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = dend_gr$labels, aes(x, y, label = label),
            hjust = 1, angle = 0, size = 3) +
  ylim(-2, 20) + 
  xlab("ESTAÇÕES") +
  ylab("Dist Euclidiana (Ward)") +
  theme_dendro()



#REFERENCIAS
#http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning
#https://www.datanovia.com/en/blog/cluster-analysis-in-r-simplified-and-enhanced/
#https://www.datanovia.com/en/blog/clustering-example-4-steps-you-should-know/
#http://www.sthda.com/sthda/ebooks/clustering_english_edition1_preview.pdf
#https://www.datanovia.com/en/lessons/determining-the-optimal-number-of-clusters-3-must-know-methods/
#https://stat.ethz.ch/R-manual/R-devel/library/stats/html/cophenetic.html
#https://www.datanovia.com/en/lessons/heatmap-in-r-static-and-interactive-visualization/
#https://biocorecrg.github.io/CRG_RIntroduction/heatmap-2-function-from-gplots-package.html
# https://rpubs.com/cleviab/877722