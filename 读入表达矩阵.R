
######2020.2.3实验代码
data932S <- read.table("E:/coding exercise/GSE112004_counts.932S.tsv")
data933S <- read.table("E:/coding exercise/GSE112004_counts.933S.tsv")
data672T <- read.table("E:/coding exercise/GSE112004_counts.672T.tsv")
data673T <- read.table("E:/coding exercise/GSE112004_counts.673T.tsv")
data936S <- read.table("E:/coding exercise/GSE112004_counts.936S.tsv")
data937S <- read.table("E:/coding exercise/GSE112004_counts.937S.tsv")
data678T <- read.table("E:/coding exercise/GSE112004_counts.678T.tsv")
data679T <- read.table("E:/coding exercise/GSE112004_counts.679T.tsv")
data032U <- read.table("E:/coding exercise/GSE112004_counts.032U.tsv")
data033U <- read.table("E:/coding exercise/GSE112004_counts.033U.tsv")
data934S <- read.table("E:/coding exercise/GSE112004_counts.934S.tsv")
data935S <- read.table("E:/coding exercise/GSE112004_counts.935S.tsv")
data674T <- read.table("E:/coding exercise/GSE112004_counts.674T.tsv")
data675T <- read.table("E:/coding exercise/GSE112004_counts.675T.tsv")
data676T <- read.table("E:/coding exercise/GSE112004_counts.676T.tsv")
data677T <- read.table("E:/coding exercise/GSE112004_counts.677T.tsv")
data030U <- read.table("E:/coding exercise/GSE112004_counts.030U.tsv")
data031U <- read.table("E:/coding exercise/GSE112004_counts.031U.tsv")
data938S <- read.table("E:/coding exercise/GSE112004_counts.938S.tsv")
data939S <- read.table("E:/coding exercise/GSE112004_counts.939S.tsv")

write.csv(data938S,file = "data938S.csv")
write.csv(data939S,file = "data939S.csv")
# 照着这样转换数据格式,导出文件
# 之后手动在第一列加gene_id

#构建表达矩阵的另一个方法
raw_count1 <- merge(merge(merge(data672T,data673T,by = "gene_id"),data932S,by = "gene_id"),data933S,by = "gene_id")


#colData是一个dataframe，第一列是样品名称，第二列是样品的处理情况，可以自己建立
condition1.1 <- rep("rep",times=384) 
condition1.2 <- rep("condition",times=384)
condition1 <- c(condition1.1,condition1.2)#样品的处理情况

count_data1 <- raw_count1[,2:769]
row.names(count_data1) <- raw_count1[,1]
col_data1 <- data.frame(row.names = colnames(count_data1),condition1)#得到的矩阵看上去没什么毛病但是没办法成dds矩阵



###2020.2.4代码
data938S <- read.csv("E:/coding exercise/data938S.csv")
data939S <- read.csv("E:/coding exercise/data939S.csv")
data030U <- read.csv("E:/coding exercise/data030U.csv")
data031U <- read.csv("E:/coding exercise/data031U.csv")
raw_count1 <- merge(merge(merge(data030U,data031U,by = "gene_id"),data938S,by = "gene_id"),data939S,by = "gene_id")

#注意geneid要取整，不能带小数点，不然比对不到（这步我还没研究，等等再说）
# 第一步将匹配到的.以及后面的数字连续匹配并替换为空，并赋值给ENSEMBL
ENSEMBL <- gsub("\\.\\d*", "", raw_count1$gene_id)
# 将ENSEMBL重新添加到raw_count_filt1矩阵
row.names(raw_count1) <- ENSEMBL
head(raw_count1)

raw_count1<-raw_count1[,-1]
#然后需要给数据集的每项加1，在后续的DESeq分析中如果某个基因的表达量含有0就会被筛掉，避免全部被筛所以+1
#之后需要筛掉原本表达量全部为零的基因，也就是加完后表达量总和是是768的基因。
plus <- matrix(1,nrow=18150,ncol=768)
raw_count1 <- raw_count1+plus


#建立coldata矩阵，存储信息
condition <- factor(c(rep("control",384),rep("treat",384)))
cell_name1 <- factor(colnames(raw_count1))
coldata1 <- data.frame(row.names=colnames(raw_count1), condition)

#建立dds矩阵(去DESeq2差异表达分析脚本)

#####2020.2.21
data678T <- read.csv("E:/coding exercise/data678T.csv")
data679T <- read.csv("E:/coding exercise/data679T.csv")
data032U <- read.csv("E:/coding exercise/data032U.csv")
data033U <- read.csv("E:/coding exercise/data033U.csv")
data934S <- read.csv("E:/coding exercise/data934S.csv")
data935S <- read.csv("E:/coding exercise/data935S.csv")

data42h <- full_join(data678T,data679T,by = "gene_id")
data66h <- full_join(data032U,data033U,by = "gene_id")
data114h <- full_join(data934S,data935S,by = "gene_id")
dataMA <- full_join(full_join(data032U,data678T,by = "gene_id"),data934S, by = "gene_id")

data674T <- read.csv("E:/coding exercise/data674T.csv")
data675T <- read.csv("E:/coding exercise/data675T.csv")
data676T <- read.csv("E:/coding exercise/data676T.csv")
data677T <- read.csv("E:/coding exercise/data677T.csv")
data030U <- read.csv("E:/coding exercise/data030U.csv")
data031U <- read.csv("E:/coding exercise/data031U.csv")
data938S <- read.csv("E:/coding exercise/data938S.csv")
data939S <- read.csv("E:/coding exercise/data939S.csv")

dataD2 <- full_join(data674T,data675T, by = "gene_id")
dataD4 <- full_join(data676T,data677T, by = "gene_id")
dataD6 <- full_join(data030U,data031U, by = "gene_id")
dataD8 <- full_join(data938S,data939S, by = "gene_id")
dataIPSC <- full_join(full_join(full_join(data674T,data676T,by = "gene_id"),data030U, by = "gene_id"),data938S, by = "gene_id")

dataMA <- read.csv("E:/coding exercise/dataMA.csv")
dataIPSC <- read.csv("E:/coding exercise/dataIPSC.csv")
raw_count1 <- merge(dataMA,dataIPSC,by = "gene_id")
# 第一步将匹配到的.以及后面的数字连续匹配并替换为空，并赋值给ENSEMBL
ENSEMBL <- gsub("\\.\\d*", "", raw_count1$gene_id)
# 将ENSEMBL重新添加到raw_count_filt1矩阵
row.names(raw_count1) <- ENSEMBL
head(raw_count1)

raw_count1<-raw_count1[,-1]
#然后需要给数据集的每项加1，在后续的DESeq分析中如果某个基因的表达量含有0就会被筛掉，避免全部被筛所以+1
#之后需要筛掉原本表达量全部为零的基因，也就是加完后表达量总和是是768的基因。
plus <- matrix(1,nrow=21780,ncol=1344)
raw_count1 <- raw_count1+plus


#建立coldata矩阵，存储信息
condition <- factor(c(rep("control",576),rep("treat",768)))
cell_name1 <- factor(colnames(raw_count1))
coldata1 <- data.frame(row.names=colnames(raw_count1), condition)

#建立dds矩阵(去DESeq2差异表达分析脚本)
##替换NA为1
raw_count1[is.na(raw_count1)] <- 1


###2020.3.7 进行共同高表达基因筛选
library(dplyr)

data678T <- read.csv("E:/coding exercise/data678T.csv")
data679T <- read.csv("E:/coding exercise/data679T.csv")
data032U <- read.csv("E:/coding exercise/data032U.csv")
data033U <- read.csv("E:/coding exercise/data033U.csv")
data934S <- read.csv("E:/coding exercise/data934S.csv")
data935S <- read.csv("E:/coding exercise/data935S.csv")

data42h <- full_join(data678T,data679T,by = "gene_id")
data66h <- full_join(data032U,data033U,by = "gene_id")
data114h <- full_join(data934S,data935S,by = "gene_id")
dataMA <- full_join(full_join(data42h,data66h,by = "gene_id"),data114h, by = "gene_id")

data674T <- read.csv("E:/coding exercise/data674T.csv")
data675T <- read.csv("E:/coding exercise/data675T.csv")
data676T <- read.csv("E:/coding exercise/data676T.csv")
data677T <- read.csv("E:/coding exercise/data677T.csv")
data030U <- read.csv("E:/coding exercise/data030U.csv")
data031U <- read.csv("E:/coding exercise/data031U.csv")
data938S <- read.csv("E:/coding exercise/data938S.csv")
data939S <- read.csv("E:/coding exercise/data939S.csv")

dataD2 <- full_join(data674T,data675T, by = "gene_id")
dataD4 <- full_join(data676T,data677T, by = "gene_id")
dataD6 <- full_join(data030U,data031U, by = "gene_id")
dataD8 <- full_join(data938S,data939S, by = "gene_id")
dataIPSC <- full_join(full_join(full_join(dataD2,dataD4,by = "gene_id"),dataD6, by = "gene_id"),dataD8, by = "gene_id")

#取并集
raw_count1 <- inner_join(dataMA,dataIPSC,by = "gene_id") #join data, retain only rows in both sets.
# 把前面的geneid取整
ENSEMBL <- gsub("\\.\\d*", "", raw_count1$gene_id)
# 将ENSEMBL重新添加到raw_count_filt1矩阵
row.names(raw_count1) <- ENSEMBL
head(raw_count1)
raw_count1<-raw_count1[,-1]
##替换NA为0
raw_count1[is.na(raw_count1)] <- 0

raw_count1 <- as.data.frame(raw_count1)
raw_filter1 <- raw_count1[rowSums(raw_count1)>21504, ]

final <- cbind(raw_filter1, sum = rowSums(raw_filter1))
geneid <- row.names(final)
final.1 <- data.frame(geneid, final[,2689])
tansid <- select(org.Mm.eg.db,keys = geneid, columns = c("ENTREZID"),keytype = "ENSEMBL")

write.csv(tansid, file = "high-gene-transid.csv")
write.csv(final.1, file = "high-gene.csv")
