source("http://bioconductor.org/biocLite.R")#载入安装工具bioconductor
BiobiocLite("DESeq2")
library(DESeq2)

setwd("E:/coding exercise")

dds1 <- DESeqDataSetFromMatrix(countData = raw_count1, colData = coldata1, design= ~ condition)
#dds1 <- DESeq(dds1)
#过滤低质量的低count数据,这里把参数定为所有样本基因表达量之和小于700(一半(另外加了1344)
dds_filter1 <- dds1[rowSums(counts(dds1))>2044, ]

library("BiocParallel")
register(MulticoreParam(4))
register(SnowParam(4))

dds1$condition <- relevel(dds1$condition, ref = "control")
dds_out1 <- DESeq(dds_filter1)# 这步大概要用3h
res1 <- results(dds_out1)
summary(res1)
table(res1$padj<0.05)
res_deseq1 <- res1[order(res1$padj),]

diff_gene_deseq1 <-subset(res_deseq1, padj < 0.05 & abs(log2FoldChange) > 1)
#或 diff_gene_deseq2 <-subset(res_deseq1,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
dim(diff_gene_deseq1)
head(diff_gene_deseq1)
write.csv(diff_gene_deseq1,file= "diffgene_all.csv")






