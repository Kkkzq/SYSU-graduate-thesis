library(ggplot2)
library(Cairo)

setwd("E:/coding exercise")
GO_1_BP <- read.table("E:/coding exercise/GO.1.BP.tsv", header = T, sep = "\t")
###
png_path = "./GO_1_BP.png"
CairoPNG(png_path, width = 12, height = 7, units='in',dpi=600)
###这两步是自动输出图像的，还挺好看，应该下面的ggplot用不上了

ggplot(data=GO_1_BP)+
  geom_bar(aes(x=reorder(Term,Count),y=Count, fill=-log10(PValue)), stat='identity') +
  coord_flip() +
  scale_fill_gradient(expression(-log["10"](P.value)),low="red", high = "blue") +
  xlab("") +
  ylab("Gene count") +
  scale_y_continuous(expand=c(0, 0))+
  theme(
    axis.text.x=element_text(color="black",size=rel(1.5)),
    axis.text.y=element_text(color="black", size=rel(1.6)),
    axis.title.x = element_text(color="black", size=rel(1.6)),
    legend.text=element_text(color="black",size=rel(1.0)),
    legend.title = element_text(color="black",size=rel(1.1))
    # legend.position=c(0,1),legend.justification=c(-1,0)
    # legend.position="top",
  )



