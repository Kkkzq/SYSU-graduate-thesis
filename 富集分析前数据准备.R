library(dplyr)
# 富集分析前读入差异基因，把它们合在一起
# 我的GSE112004分成三个阶段。0-18h，18-114h,18-d8
diffgene1 <- read.csv("E:/coding exercise/result/diffgene_0h6h.csv")
diffgene2 <- read.csv("E:/coding exercise/result/diffgene_6h18h.csv")

genelist1 <- diffgene1[,1]
genelist2 <- diffgene2[,1]

Gene <- union(genelist1,genelist2)#取一下基因id的并集
# 这一步做完后去做tansid转换一下id，这是准备工作，之后再接着做GO分析

diffgene1 <- read.csv("E:/coding exercise/result/diffgene_18h42h.csv")
diffgene2 <- read.csv("E:/coding exercise/result/diffgene_42h66h.csv")
diffgene3 <- read.csv("E:/coding exercise/result/diffgene_66h114h.csv")

genelist1 <- diffgene1[,1]
genelist2 <- diffgene2[,1]
genelist3 <- diffgene3[,1]
Gene <- union(union(genelist1,genelist2),genelist3)

diffgene1 <- read.csv("E:/coding exercise/result/diffgene_18hd2.csv")
diffgene2 <- read.csv("E:/coding exercise/result/diffgene_d2d4.csv")
diffgene3 <- read.csv("E:/coding exercise/result/diffgene_d4d6.csv")
diffgene4 <- read.csv("E:/coding exercise/result/diffgene_d6d8.csv")
genelist1 <- diffgene1[,1]
genelist2 <- diffgene2[,1]
genelist3 <- diffgene3[,1]
genelist4 <- diffgene4[,1]
Gene <- union(union(union(genelist1,genelist2),genelist3),genelist4)






