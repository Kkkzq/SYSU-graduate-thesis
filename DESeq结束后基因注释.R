#这一步是进行基因注释之前生成的csv文件需要在第一列添加gene_id列名，然后读进来
diffgene <- read.csv("E:/coding exercise/result/diffgene_all.csv")
#用bioMart对差异表达基因进行注释，可以生成基因列表

#注意geneid要取整，不能带小数点，不然比对不到（这步我还没研究，等等再说）
# 第一步将匹配到的.以及后面的数字连续匹配并替换为空，并赋值给ENSEMBL
ENSEMBL <- gsub("\\.\\d*", "", diffgene$gene_id) 
# 将ENSEMBL重新添加到raw_count_filt1矩阵
row.names(diffgene) <- ENSEMBL
head(diffgene)

#在线工具也可以解决
#这一步是通过R包解决
BiocManager::install("biomaRt")
library("biomaRt")
library("curl")
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
my_ensembl_gene_id<-diffgene[,1]
#listAttributes(mart)
mms_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"),
                    filters = 'ensembl_gene_id', values = my_ensembl_gene_id, mart = mart)
write.csv(mms_symbols,file= "descripgene_all.csv")


