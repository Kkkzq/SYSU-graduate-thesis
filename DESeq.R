#import data
setwd("E:/coding exercise/other")
mycounts<-read.csv("countdata.csv")

head(mycounts)
rownames(mycounts)<-mycounts[,1]
#把带gene的列删除
mycounts<-mycounts[,-1]
head(mycounts)

condition <- factor(c(rep("control",2),rep("treat",2)))
coldata <- data.frame(row.names=colnames(mycounts), condition)
colData <- read.csv("coldata.csv")
rownames(coldata)<-colData[,1]


dds <- DESeqDataSetFromMatrix(countData = mycounts, colData = coldata, design= ~ condition)
dds <- DESeq(dds)
# 查看一下dds的内容
dds

res = results(dds, contrast=c("condition", "control", "treat"))
res = res[order(res$pvalue),]
head(res)#总体结果查看

summary(res)

diff_gene_deseq2 <-subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
#或
#> diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
dim(diff_gene_deseq2)
head(diff_gene_deseq2)
###write.csv(diff_gene_deseq2,file= "DEG_treat_vs_control.csv")###提取差异表达基因列表


