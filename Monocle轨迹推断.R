#安装monocle3
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install()

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor'))
install.packages("devtools")
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')

library(monocle3)
library(ggplot2)
library(dplyr)

######################################  读入表达矩阵数据
data932S <- read.csv("E:/coding exercise/data932S.csv")
data933S <- read.csv("E:/coding exercise/data933S.csv")
data672T <- read.csv("E:/coding exercise/data672T.csv")
data673T <- read.csv("E:/coding exercise/data673T.csv")
data936S <- read.csv("E:/coding exercise/data936S.csv")
data937S <- read.csv("E:/coding exercise/data937S.csv")
data678T <- read.csv("E:/coding exercise/data678T.csv")
data679T <- read.csv("E:/coding exercise/data679T.csv")
data032U <- read.csv("E:/coding exercise/data032U.csv")
data033U <- read.csv("E:/coding exercise/data033U.csv")
data934S <- read.csv("E:/coding exercise/data934S.csv")
data935S <- read.csv("E:/coding exercise/data935S.csv")
data674T <- read.csv("E:/coding exercise/data674T.csv")
data675T <- read.csv("E:/coding exercise/data675T.csv")
data676T <- read.csv("E:/coding exercise/data676T.csv")
data677T <- read.csv("E:/coding exercise/data677T.csv")
data030U <- read.csv("E:/coding exercise/data030U.csv")
data031U <- read.csv("E:/coding exercise/data031U.csv")
data938S <- read.csv("E:/coding exercise/data938S.csv")
data939S <- read.csv("E:/coding exercise/data939S.csv")

data0h <- full_join(data932S,data933S,by = "gene_id")
data6h <- full_join(data672T,data673T,by = "gene_id")
data18h <- full_join(data936S,data937S,by = "gene_id")
matrix1 <- full_join(full_join(data0h,data6h,by = "gene_id"),data18h,by = "gene_id")
data42h <- full_join(data678T,data679T,by = "gene_id")
data66h <- full_join(data032U,data033U,by = "gene_id")
data114h <- full_join(data934S,data935S,by = "gene_id")
matrix2 <- full_join(full_join(data42h,data66h,by = "gene_id"),data114h,by = "gene_id")
dataD2 <- full_join(data674T,data675T, by = "gene_id")
dataD4 <- full_join(data676T,data677T, by = "gene_id")
dataD6 <- full_join(data030U,data031U, by = "gene_id")
dataD8 <- full_join(data938S,data939S, by = "gene_id")
matrix3 <- full_join(full_join(full_join(dataD2,dataD4,by = "gene_id"),dataD6,by = "gene_id"),dataD8,by = "gene_id")

expression_matrix <- full_join(full_join(matrix1,matrix2,by = "gene_id"),matrix3,by = "gene_id")
# 把前面的geneid取整
ENSEMBL <- gsub("\\.\\d*", "", expression_matrix$gene_id)
# 将ENSEMBL重新添加到raw_count_filt1矩阵
row.names(expression_matrix) <- ENSEMBL
head(expression_matrix)
expression_matrix<-expression_matrix[,-1]
##替换NA为0
expression_matrix[is.na(expression_matrix)] <- 0
# 导出
write.csv(expression_matrix,file = "expression_matrix.csv")

####################################### 建立gene metadata数据框
library(clusterProfiler)
require(DOSE)
library(DO.db)
library(GO.db)
library(topGO)
library(GSEABase)
library(org.Mm.eg.db)

gene <- row.names(expression_matrix)
write.csv(gene,file = "geneid.csv")
# tansid <- select(org.Mm.eg.db,keys = gene,columns = c("SYMBOL"),keytype = "ENSEMBL")
# 选SYMBOL那一栏,但是这样转换完id对应不上，多出来了
# 考虑用网站转换，手动建数据框，在网页上转换完添加进数据框，记得命名行名
geneid<- read.csv("geneid.csv")
gene_short_name <- geneid[,2]
gene_metadata <- data.frame(row.names=geneid[,1], gene_short_name)

################################  建立cell_metadata数据框
# 例子
cell_metadata <- readRDS(url("http://staff.washington.edu/hpliner/data/cao_l2_colData.rds"))
# 决定自己建立一个包含time、treatment、induction和fate的数据框

time <- factor(c(rep("0h",384),rep("6h",384),rep("18h",384),rep("42h",384),
                      rep("66h",384),rep("114h",384),rep("D2",384),rep("D4",384),
                      rep("D6",384),rep("D8",384)))
treatment <- factor(c(rep("n/a",1152),rep("transdifferentiation",1152),rep("reprogramming",1536)))
induction <- factor(c(rep("c/EBPa",2304),rep("OSKM",1536)))
fate <- factor(c(rep("n/a",1152),rep("iMac",1152),rep("iPSCs",1536)))
cell_metadata <- data.frame(row.names=colnames(expression_matrix),time,treatment,induction,fate)


### 建立cds矩阵
### expression_matrix <- as.matrix(expression_matrix)
expression_matrix <- as(expression_matrix, "sparseMatrix")
# cds <- new_cell_data_set(expression_matrix,
#                         cell_metadata = cell_metadata,
#                         gene_metadata = gene_metadata)
# 稀疏矩阵
library(Matrix)
cds <- new_cell_data_set(as(expression_matrix, "sparseMatrix"),
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_metadata)
cds

######################## 细胞聚类与分类
###### 预处理数据
cds <- preprocess_cds(cds, num_dim = 100)
# 检查一下是否使用了足够多的PC来捕获数据集中所有细胞中基因表达的大部分变异
plot_pc_variance_explained(cds) 

###
library(Cairo)
png_path = "./monocle.png"
CairoPNG(png_path, width = 10, height = 7, units='in',dpi=600)
### 这两步是自动输出图像的

### 非线性降维
cds <- reduce_dimension(cds, preprocess_method = "PCA",reduction_method = c("UMAP"))
plot_cells(cds)
### 线性降维
cds <- reduce_dimension(cds, reduction_method="tSNE")
# 然后对tSNE结果可视化
plot_cells(cds, reduction_method="tSNE", color_cells_by="time")
plot_cells(cds, reduction_method="tSNE",color_cells_by="time", label_cell_groups=FALSE)

#### 实际上，您可以在同一cds对象上使用UMAP和t-SNE-一个对象不会覆盖另一个对象的结果。
#### 但是您必须在下游函数（如）中指定所需的对象plot_cells。

### 将细胞分组
### 将cell分组到cluster是识别数据中表示的cell类型的重要步骤。
cds = cluster_cells(cds,resolution=1e-5,reduction_method="tSNE")
cds@clusters$tSNE$clusters # 看一下有多少类-18类；UMAP是24类
# 检查有多少partition，partition是将所述细胞分化成更大，更良好分离的组
cds@clusters$tSNE$partitions # 结果只有1个，那就不用了，以后用cluster

plot_cells(cds,reduction_method="tSNE",label_cell_groups=FALSE, color_cells_by="cluster")
plot_cells(cds,label_cell_groups=FALSE, color_cells_by="cluster")
# P.S.color_cells_by的参数可以去掉

# 查找每个簇表达的标记基因，记得标注用tsne降维
marker_test_res <- top_markers(cds, reduction_method="tSNE", group_cells_by="cluster", 
                               reference_cells=1000, cores=8)
marker_test_res[1:4,1:4]
top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(3, pseudo_R2)
top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))
# 可视化
plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    reduction_method="tSNE",
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=3)

# 对cell进行注释
# 先将partitions的分组由因子型转为字符型
colData(cds)$assigned_cell_type = as.character(cds@clusters$tSNE$clusters)
# 再对字符型重新定义
colData(cds)$assigned_cell_type = dplyr::recode(colData(cds)$assigned_cell_type, reduction_method="tSNE",
                                                "1"="pre-B cell fate",
                                                "2"="pre-B cell fate",
                                                "3"="iMac fate",
                                                "4"="iMac fate",
                                                "5"="iPSc fate",
                                                "6"="Unknown cell type",
                                                "7"="pre-B cell fate",
                                                "8"="pre-B cell fate",
                                                "9"="iMac fate",
                                                "10"="iPSc fate",
                                                "11"="iMac fate",
                                                "12"="Unknown cell type",
                                                "13"="iPSc fate",
                                                "14"="Failed QC",
                                                "15"="iPSc fate",
                                                "16"="iMac fate",
                                                "17"="iPSc fate",
                                                "18"="iMac fate")
plot_cells(cds, group_cells_by="cluster",reduction_method="tSNE", color_cells_by="assigned_cell_type")

# 基于Garnett的自动化注释，但是对于我的分离不太完善的数据集而言似乎并不好用
# step-1：首先根据上面得到的细胞类型assigned_cell_type，找top_marker
assigned_type_marker_test_res = top_markers(cds,
                                            reduction_method="tSNE",
                                            group_cells_by="cluster",
                                            reference_cells=1000,
                                            cores=8)
# step-2：过滤（阈值自定义）
garnett_markers = assigned_type_marker_test_res %>%
  filter(marker_test_q_value < 0.01 & specificity >= 0.5) %>%
  group_by(cell_group) %>%
  top_n(5, marker_score)
# step-3：去重复
garnett_markers = garnett_markers %>% group_by(gene_short_name) %>%
  filter(n() == 1)
# step-4：生成marker文件
generate_garnett_marker_file(garnett_markers, file="./marker_file.txt")


################################## 轨迹推断（使用UMAP降维）
# 加载数据集。覆盖原来的
cds <- new_cell_data_set(as(expression_matrix, "sparseMatrix"),
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_metadata)
cds <- preprocess_cds(cds, num_dim = 50)
cds <- align_cds(cds, alignment_group = "batch", residual_model_formula_str = "~ bg.300.loading + bg.400.loading + bg.500.1.loading + bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading")

cds <- reduce_dimension(cds)
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "time", label_cell_groups=FALSE)

# 可以看某些基因的表达模式
ciliated_genes <- c("Hsp90aa1",
                    "Ccl9",
                    "Ccl6",
                    "Rac2")
# Rac2、F13a1、Sparc、C3、Itgb2、Hspb1、Eif2s2、Mylpf、Vim、Tyrobp
ciliated_genes <- c("F13a1",
                    "Sparc",
                    "C3",
                    "Itgb2")
ciliated_genes <- c("Hspb1",
                    "Eif2s2",
                    "Mylpf",
                    "Vim",
                    "Tyrobp")

plot_cells(cds,
           genes=ciliated_genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

# 细胞cluster
cds <- cluster_cells(cds)
cds@clusters$UMAP$clusters
cds@clusters$UMAP$partitions
plot_cells(cds, color_cells_by = "cluster")

cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "fate",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_cell_groups=FALSE,
           label_branch_points=FALSE)

# 在伪时间内对单元进行排序
plot_cells(cds,
           color_cells_by = "time",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)
# 指定根节点(手动)
cds <- order_cells(cds)

# 写函数找根节点
# 官方给出了一个函数，这里定义了一个time_bin，选择了最早的时间点区间。
get_earliest_principal_node <- function(cds, time_bin = "0h"){
  # 首先找到出现在最早时间区间的细胞ID
  cell_ids <- which(colData(cds)[, "time"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds = order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

# 通过伪时间对其进行着色可以显示cell的排序方式
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

#################################### 对得到的分支进行分析
# 测试轨迹上相似位置的细胞是否具有相关表达
ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))

plot_cells(cds, genes=c("F13a1",
                        "Sparc",
                        "C3",
                        "Itgb2"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)

# 可以将轨迹可变基因收集到模块中
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=c(0,10^seq(-6,-1)))
# 绘制了每组细胞类型内的模块总得分
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=colData(cds)$cell.type)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat,
                   scale="column", clustering_method="ward.D2")

# Monocle提供了另一种绘图功能，有时可以沿一条路径更清晰地观察基因的动力学。
iPSCs_genes <- c("Mylpf",
               "Sparc",
               "Hspb1",
               "Eif2s2",
               "Hsp90aa1")
iPSCs_lineage_cds <- cds[rowData(cds)$gene_short_name %in% iPSCs_genes,
                       colData(cds)$fate %in% c("n/a","iPSCs")]
plot_genes_in_pseudotime(iPSCs_lineage_cds,
                         color_cells_by="time",
                         min_expr=0.5)

iMac_genes1 <- c("Ccl6",
                "Ccl9",
                "Rac2")
iMac_lineage1_cds <- cds[rowData(cds)$gene_short_name %in% iMac_genes1,
                         colData(cds)$fate %in% c("n/a","iMac")]
plot_genes_in_pseudotime(iMac_lineage1_cds,
                         color_cells_by="time",
                         min_expr=0.5)

iMac_genes2 <- c("F13a1",
                 "C3",
                 "Itgb2",
                 "Vim",
                 "Tyrobp")
iMac_lineage2_cds <- cds[rowData(cds)$gene_short_name %in% iMac_genes2,
                         colData(cds)$fate %in% c("n/a","iMac")]
plot_genes_in_pseudotime(iMac_lineage2_cds,
                         color_cells_by="time",
                         min_expr=0.5)














