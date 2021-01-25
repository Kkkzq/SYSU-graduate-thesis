#这一步是把差异表达分析的结果可视化

# MA图，也叫 mean-difference plot或者Bland-Altman plot，用来估计模型中系数的分布。
# M表示log fold change，衡量基因表达量变化，上调还是下调。
# A表示每个基因的count的均值
library(ggplot2)
library(BiocGenerics)
plotMA(res1,ylim=c(-4,4))
topGene <- rownames(res1)[which.min(res1$padj)]
with(res1[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})

#lfcShrink 收缩log2 fold change(好像耗时相当久，和差异表达分析耗时差不多)
resLFC <- lfcShrink(dds_out1,coef = 2,res1=res1)
plotMA(resLFC, ylim=c(-5,5))
topGene <- rownames(resLFC)[which.min(res1$padj)]
with(resLFC[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})
idx <- identify(res1$baseMean, res1$log2FoldChange)

# gene 聚类热图
library(genefilter)
library(pheatmap)

rld <- rlogTransformation(dds_out1,blind = F)
write.csv(assay(rld),file="mm.DESeq2.pseudo.counts-1.csv")

topVarGene <- head(order(rowVars(assay(rld)),decreasing = TRUE),20)
mat  <- assay(rld)[ topVarGene, ]
### mat <- mat - rowMeans(mat) 减去一个平均值，让数值更加集中。第二个图
anno <- as.data.frame(colData(rld)[,c("condition","sizeFactor")])
pheatmap(mat, annotation_col = anno)

# 火山图(不打算再使用了)
### 火山图
data <-read.table(file="volcano.txt",header = TRUE, row.names =1,sep = "\t")
volcano <-ggplot(data = volcano_data,aes(x=log2FoldChange,y= -1*log10(padj)))+geom_point(aes(color=significant))+scale_color_manual(values = c("red","grey","blue")) + labs(title="Volcano_Plot",x=expression((log[2](FC)), y=expression(-log[10](padj)) ))+geom_hline(yintercept=1.3,linetype=4)+geom_vline(xintercept=c(-1,1),linetype=4)
volcano

