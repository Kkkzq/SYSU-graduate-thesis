###########   生信作业   #########

######  获取数据，提取表达基因  #######
.libPaths("C:\\Anaconda3\\Lib\\site-packages\\rpy")
setwd("E:\\R-project\\homework-bioinfo")
library(GEOquery)

#######  下载现成的表达矩阵
options( 'download.file.method.GEOquery' = 'libcurl' )
downGSE <- function(studyID = "GSE1009", destdir = ".") {
  eSet <- getGEO(studyID, destdir = destdir, getGPL = F)
  exprSet = exprs(eSet[[1]])
  pdata = pData(eSet[[1]])
  write.csv(exprSet, paste0(studyID, "_exprSet.csv"))
  write.csv(pdata, paste0(studyID, "_metadata.csv"))
  return(eSet)
}

a <- downGSE(studyID = "GSE1643")  ###肺癌
b <- downGSE(studyID = "GSE3280")  ###白血病
c <- downGSE(studyID = "GSE4587")  ###黑色素瘤

########  进行聚类  #########
lung_data <- read.csv('GSE1643_exprSet.csv')
rownames(lung_data) <- lung_data[,1]
lung <- lung_data[,-1]

white <- read.csv("GSE3280_exprSet.csv")
rownames(white) <- white[,1]
white <- white[,-1]

black <- read.csv("GSE4587_exprSet.csv")
rownames(black) <- black[,1]
black <- black[,-1]

data <- cbind(black,lung,white)
data0 <- apply(t(data),1,scale)
rownames(data0) <- rownames(data)



library(cluster)
library(fpc)

fit <- kmeans(data0,centers = 3,iter.max = 100, nstart = 25)

fit_cluster <- fit$cluster
clusplot(data0,fit_cluster,shade = T)


fit_pam <- pamk(data2,krange = 2:6,critout = T)







BiocManager :: install("impute")
install.packages("WGCNA")
install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival")) 
source("http://bioconductor.org/biocLite.R") 
biocLite(c("GO.db", "preprocessCore", "impute")) 
library(WGCNA)



samples=read.csv('GSE1643_metadata.csv')
dim(samples)

expro=read.csv('GSE1643_exprSet.csv')
row.names(expro) <- expro[,1]
expro <- expro[,-1]
dim(expro)

m.vars=apply(expro,1,var)
expro.upper=expro[which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[4]),]

datExpr=as.data.frame(t(expro.upper));
gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK
sampleTree = hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers"
     , sub="", xlab="")











#########   画图练习  ###########
library(igraph)
set.seed(1)
# 先构建一个有向网络临近矩阵

# create a 5x5 adjacency matrix
adj <- matrix(sample(c(0,1), 25, replace = T), nrow = 5)

# set diagonal to zero (no self-loops)
diag(adj) <- 0

# plotting with igraph
g <- graph.adjacency(adj, mode = "directed")
plot(g)

#再构建一个无向网络
adj[upper.tri(adj)] <- 0 #就是取矩阵的下半部分
g <- graph.adjacency(adj, mode = "undirected")
plot(g)

#再构建加权网络 weighted network
set.seed(1)
adj <- matrix(rnorm(25, mean = 3.5, sd = 5), nrow = 5)
adj[upper.tri(adj, diag = TRUE)] <- 0

# note that igraph ignores edges with negative weights
g <- graph.adjacency(adj, mode = "undirected", weighted = T)
plot(g, edge.width = E(g)$weight)
