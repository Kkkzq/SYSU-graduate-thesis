BiocManager::install("GEOquery")
library(GEOquery)
library(Biobase)
library(dplyr)
library(limma)
setwd("E:/coding exercise/")
browseVignettes("GEOquery")##获取帮助

options( 'download.file.method.GEOquery' = 'libcurl' ) 
geneset <- getGEO('GSE112004',destdir = ".",
               AnnotGPL = F,
               getGPL = F)
save(geneset,file = 'GSE112004.gset.Rdata')

#更换数据下载方法
geneset <- getGEO("GSE112004", GSEMatrix =TRUE, AnnotGPL=TRUE )
show(geneset)

##输入表达矩阵
# 使用GEOquery
expr_Set <- exprs(geneset[[1]])
# 基于matrix
expr.df <- read.table(file = "GSE112004_series_matrix.txt", header =TRUE,
                      comment.char = "!", row.names=1)

##输入分组对象，实验设计矩阵
pData <- pData(geneset[[1]])
#构建实验设计的数据框
sample <- pData$geo_accession
treat_time <- rep(c("day6","66h","6h","day2","day4","42h","0h","114h","18h","day8"),each=384)
treat_type <- rep(c("OSKM","C/EBPα","OSKM","C/EBPα","OSKM"), c(384,768,768,1536,384))
design_df <- data.frame(sample, treat_time, treat_type)

TS <- paste(design_df$treat_time, design_df$treat_type, sep=".")
TS
TS <- factor(TS, levels = unique(TS))
design <- model.matrix(~0+TS)
fit <- lmFit(expr.df, design)



