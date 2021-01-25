#####需要的包
#加载细胞轨迹分析需要的包
library("Seurat")
library("dplyr")
install.packages("ggplot2")
library("ggplot2")
library("monocle")



######细胞提取
#细胞提取和合并
#加载数据
load("/data1/genomes/Right2.RData")
#细胞提取
#dtm
DTM.new <- RenameIdents(DTM.new, `Villous cytotrophoblast` = 'VCT.DTM', 
                    `Extravillous trophoblast` = "EVT.DTM", 
                    `Syncytiotrophoblast` = "SCT.DTM") 
VCT.DTM <- DTM.new[,names(DTM.new@active.ident[DTM.new@active.ident == 'VCT.DTM'])]
EVT.DTM <- DTM.new[,names(DTM.new@active.ident[DTM.new@active.ident == 'EVT.DTM'])]
SCT.DTM <- DTM.new[,names(DTM.new@active.ident[DTM.new@active.ident == 'SCT.DTM'])]
#dtp
DTP.new <- RenameIdents(DTP.new, `Villous cytotrophoblast` = 'VCT.DTP', 
                    `Syncytiotrophoblast` = "SCT.DTP") 
VCT.DTP <- DTP.new[,names(DTP.new@active.ident[DTP.new@active.ident == 'VCT.DTP'])]
SCT.DTP <- DTP.new[,names(DTP.new@active.ident[DTP.new@active.ident == 'SCT.DTP'])]
#xtm
XTM.new <- RenameIdents(XTM.new, `Villous cytotrophoblast` = 'VCT.XTM', 
                    `Extravillous trophoblast` = "EVT.XTM", 
                    `Syncytiotrophoblast` = "SCT.XTM") 
VCT.XTM <- XTM.new[,names(XTM.new@active.ident[XTM.new@active.ident == 'VCT.XTM'])]
EVT.XTM <- XTM.new[,names(XTM.new@active.ident[XTM.new@active.ident == 'EVT.XTM'])]
SCT.XTM <- XTM.new[,names(XTM.new@active.ident[XTM.new@active.ident == 'SCT.XTM'])]
#xtp
XTP.new <- RenameIdents(XTP.new, `Villous cytotrophoblast` = 'VCT.XTP', 
                    `Extravillous trophoblast` = "EVT.XTP", 
                    `Syncytiotrophoblast` = "SCT.XTP") 
VCT.XTP <- XTP.new[,names(XTP.new@active.ident[XTP.new@active.ident == 'VCT.XTP'])]
EVT.XTP <- XTP.new[,names(XTP.new@active.ident[XTP.new@active.ident == 'EVT.XTP'])]
SCT.XTP <- XTP.new[,names(XTP.new@active.ident[XTP.new@active.ident == 'SCT.XTP'])]
#合并
#别忘记加双引号！！！
#整体的
ZYC <- merge(EVT.DTM, y = c(EVT.XTM, EVT.XTP, SCT.DTM, SCT.DTP, SCT.XTM, SCT.XTP,
                            VCT.DTM, VCT.DTP, VCT.XTM, VCT.XTP), merge.data = T, 
             add.cell.ids = c("EVT.DTM", "EVT.XTM", "EVT.XTP", "SCT.DTM", "SCT.DTP", "SCT.XTM", "SCT.XTP",
                              "VCT.DTM", "VCT.DTP", "VCT.XTM", "VCT.XTP"), project = "monocle")
#小胎的
XZYC <- merge(EVT.XTM, y = c(EVT.XTP, SCT.XTM, SCT.XTP, VCT.XTM, VCT.XTP), 
              merge.data = T, project = "monocle",
           add.cell.ids = c("EVT.XTM", "EVT.XTP", "SCT.XTM", "SCT.XTP", "VCT.XTM", "VCT.XTP"))
#大胎的
DZYC <- merge(EVT.DTM, y = c(SCT.DTM, SCT.DTP, VCT.DTM, VCT.DTP), 
              merge.data = T, project = "monocle",
              add.cell.ids = c("EVT.DTM", "SCT.DTM", "SCT.DTP","VCT.DTM", "VCT.DTP"))
#保存数据，这样之后就不用再做前面的提取细胞+合并了
saveRDS(ZYC, "/home/tangzeh/lst/monocle/seurat_data/ZYC.rds")
saveRDS(XZYC, "/home/tangzeh/lst/monocle/seurat_data/XZYC.rds")
saveRDS(DZYC, "/home/tangzeh/lst/monocle/seurat_data/DZYC.rds")


#数据加载和cds构建

######数据加载和cds构建
#加载数据，每次做选其中一个，后面都是一样的
ZYC <- readRDS("/home/tangzeh/lst/monocle/seurat_data/ZYC.rds")
ZYC <- readRDS("/home/tangzeh/lst/monocle/seurat_data/XZYC.rds")
ZYC <- readRDS("/home/tangzeh/lst/monocle/seurat_data/DZYC.rds")

#建立celldataset
#1.表达矩阵
matrix <- as.matrix(ZYC@assays[["RNA"]]@counts)
#2.pd-表型，目前来看应该是细胞ming
cell <- as.data.frame(colnames(matrix))
row.names(cell) <- colnames(matrix)
colnames(cell) <- "cell"
#3.基因名
gene <- as.data.frame(rownames(matrix))
colnames(gene) <- "gene_short_name"
row.names(gene) <- row.names(matrix)
#转换类
pd <- new("AnnotatedDataFrame", data = cell)
fd <- new("AnnotatedDataFrame", data = gene)
#创建celldataset(cds)
#####或许可以尝试lowerdetectionlimit为1
ZYC.MON <- newCellDataSet(as(matrix, "sparseMatrix"),
                          phenoData = pd,
                          featureData = fd,
                          lowerDetectionLimit = 0.5,
                          expressionFamily = negbinomial.size())

cds <- ZYC.MON


#流程分析

#####流程分析
#可能需要好好过滤一下，以减少细胞量
cds <- detectGenes(cds, min_expr = 0.1)
#设置一个基因表达量的过滤阈值
#结果会在cds@featureData@data中新增一列num_cells_expressed
#记录这个基因在多少细胞中有表达
print(head(cds@featureData@data))
#利用上面的结果过滤基因（表达该基因的细胞数目少）
expressed_genes <- row.names(subset(cds@featureData@data,
                                    num_cells_expressed >= 10))
length(expressed_genes)
cds <- cds[expressed_genes,]

#进行聚类
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
disp_table <- dispersionTable(cds) # 挑有差异的
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1) # 挑表达量不太低的
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)  # 准备聚类基因名单
plot_ordering_genes(cds) 
# 图中黑色的点就是被标记出来一会要进行聚类的基因

#选主成分
plot_pc_variance_explained(cds, return_all = F) # norm_method='log'
#根据上图选主成分
#这次的ZYC选9
#DZYC选9，XZYC选
cds <- reduceDimension(cds, max_components = 2, num_dim = 9,
                       reduction_method = 'tSNE', verbose = T)
#聚类
cds <- clusterCells(cds, num_clusters = 4) 
#这里的类好像会比num_clusters设置的值少1
plot_cell_clusters(cds, 1, 2, color = "Cluster")
diff_test_res <- differentialGeneTest(cds,
                                      fullModelFormulaStr = "~Cluster")

#注意修改保存路径
write.csv(diff_test_res, "/home/tangzeh/lst/monocle/DZYC_diff_test_res.csv")
#得到差异基因
sig_genes <- subset(diff_test_res, qval < 0.01)
write.csv(sig_genes, "/home/tangzeh/lst/monocle/DZYC_sig_genes.csv")


#把差异基因保留下来

#可以查看差异基因
head(sig_genes[,c("gene_short_name", "pval", "qval")] )
#做差异基因的图
#这只是个示范，没有真的选。如果要看的话，可以选几个明显的差异基因
#要把基因名转换成character
cg <- as.character(head(sig_genes$gene_short_name))
plot_genes_jitter(cds[cg,], 
                  grouping = "Cluster", ncol= 2,
                  color_by = "Cluster")

#推断发育轨迹
#选择合适基因
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds)
#降维
#默认使用DDRTree的方法 
cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')
#细胞排序
cds <- orderCells(cds)
plot_cell_trajectory(cds, color_by = "Cluster")
plot_genes_in_pseudotime(cds[cg,],
                         color_by = "Cluster")

#还有几个其它可视化函数
plot_cell_trajectory(cds, color_by = "State")
plot_cell_trajectory(cds, color_by = "Pseudotime")
plot_cell_trajectory(cds, color_by = "State") +
  facet_wrap(~State, nrow = 1)

