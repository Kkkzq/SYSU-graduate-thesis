#GO/KEGG分析及GSEA
# Bioconductor的包，安装都是一个套路，source一下，bioLite一下，就差不多了。
source("https://bioconductor.org/biocLite.R")
BiocManager::install("clusterProfiler")
library(clusterProfiler)
BiocManager::install("DOSE")
require(DOSE)
library(DO.db)
library(GO.db)
library(topGO)
library(GSEABase)
#安装构建自己的基因组注释数据
# 我们是小鼠数据，所以直接安装载入就可以了，当然人类的也是一样。
# 小鼠的注释数据
BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)

# GO（Gene Ontology）分析
# 看一下数据库的ID类型
keytypes(org.Mm.eg.db)
# Jimmy推荐的是使用select()函数进行ID的转换
# keys是原始的ID，columns是转换之后的ID，keytype是要指定的原始ID类型
row.names(diffgene) <- diffgene[,1]
gene <- row.names(diffgene)

# 直接从这里开始
tansid <- select(org.Mm.eg.db,keys = gene,columns = c("GENENAME","SYMBOL","ENTREZID"),keytype = "ENSEMBL")
write.csv(tansid,file= "transgene_all.csv")
## ENSEMBL                                                GENENAME SYMBOL ENTREZID
## 1 ENSMUSG00000003309 adaptor protein complex AP-1, mu 2 subunit  Ap1m2    11768
## 2 ENSMUSG00000046323    developmental pluripotency-associated 3  Dppa3    73708
# 此外还有bitr()函数可以转换ID，得到的结果都是一样的
# anyid <- bitr(gene,fromType = "ENSEMBL",toType = c("GENENAME","SYMBOL","ENTREZID"),OrgDb = org.Mm.eg.db)

# enrichGO()函数进行GO分析及画图
# 进行go分析
ego.BP <- enrichGO(gene = tansid[,1],OrgDb = org.Mm.eg.db,keyType = "ENSEMBL",ont = "BP") #时间长
ego.MF <- enrichGO(gene = tansid[,1],OrgDb = org.Mm.eg.db,keyType = "ENSEMBL",ont = "MF")
# ont:主要的分为三种，三个层面来阐述基因功能，生物学过程（BP），细胞组分（CC），分子功能（MF）
# 气泡图和bar图，记得改id
# 在画图前去GOKEGG可视化脚本里运行限制输出路径的代码
dotplot(ego.BP,showCategory=20)
barplot(ego.BP,showCategory=20)
dotplot(ego.MF,showCategory=20)
barplot(ego.MF,showCategory=20)
# 网络图
enrichMap(ego)
# GO图需要安装额外的包
BiocManager::install("topGO")
BiocManager::install("Rgraphviz")
require(Rgraphviz)
plotGOgraph(ego) #直接画会很丑


# KEGG（pathway）分析
# 转换ID适合KEGG
x=bitr(tansid[,1],fromType = "ENSEMBL",toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
# 获取keggID
kegg <- x[,2]
# KEGG分析，在KEGG官网中，物种都有对应的缩写，小鼠mmu，其他的缩写看官网：http://www.genome.jp/kegg/catalog/org_list.html
ekk <- enrichKEGG(kegg, keyType = "kegg",organism = "mmu", pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05)
head(summary(ekk))
# 气泡图
dotplot(ekk, showCategory=20)
# 将GO/KEGG结果转换成CSV格式输出
write.csv(as.data.frame(ekk),"KEGG-enrich-all.csv",row.names =F)
write.csv(as.data.frame(ego.BP),"GO.all.BP.csv",row.names =F)
write.csv(as.data.frame(ego.MF),"GO.all.MF.csv",row.names =F)
write.csv(as.data.frame(gsemf),"GSEA-enrich-BP.csv",row.names =F)

# GSEA分析
# 基因集富集分析 (Gene Set Enrichment Analysis, GSEA) 的基本思想是使用预定义的基因集（通常来自功能注释或先前实验的结果)
# 将基因按照在两类样本中的差异表达程度排序，然后检验预先设定的基因集合是否在这个排序表的顶端或者底端富集。
# 基因集合富集分析检测基因集合而不是单个基因的表达变化
# Gene Set Enrichment Analysis（GSEA）
# 获取按照log2FC大小来排序的基因列表
genelist <- diff_gene_deseq1$log2FoldChange
names(genelist) <- rownames(diffgene)
genelist <- sort(genelist, decreasing = TRUE)
# GSEA分析（具体参数参考：https://mp.weixin.qq.com/s/p-n5jq5Rx2TqDBStS2nzoQ）
gsemf <- gseGO(genelist,
               OrgDb = org.Mm.eg.db,
               keyType = "ENSEMBL",
               ont="BP"
)# 会很久
# 查看大致信息
head(gsemf)
# 画出GSEA图
gseaplot(gsemf, geneSetID="GO:0003676")


names(genelist) <- kegg
genelist <- sort(genelist, decreasing = TRUE)
kk <- gseKEGG(geneList = genelist,
              organism = 'mmu',
              nPerm = 1000,
              minGSSize = 10,
              maxGSSize = 500,
              pvalueCutoff = 1,
              verbose = FALSE)
go_result_kk <- as.data.frame(kk)
write.csv(as.data.frame(go_result_kk),"GSEA_result_kk.csv",row.names =F)
gseaplot(kk, geneSetID = "mmu00830")

