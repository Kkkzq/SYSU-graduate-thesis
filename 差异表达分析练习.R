##################尝试示例##########
library(GEOquery)
gset <- getGEO("GSE13535", GSEMatrix =TRUE, AnnotGPL=TRUE )
show(gset)

exprSet <- exprs(gset[[1]])

pData <- pData(gset[[1]])

sample <- pData$geo_accession
treat_time <- rep(c("2h","18h"),each=11)
treat_type <- rep(rep(c("vehicle_control","PE1.3_embolized","PE2.0_embolized"), c(3,4,4)),
                  times=2)
design_df <- data.frame(sample, treat_time, treat_type)

TS <- paste(design_df$treat_time, design_df$treat_type, sep=".")
TS
TS <- factor(TS, levels = unique(TS))
design <- model.matrix(~0+TS)
fit <- lmFit(exprSet, design)

cont.matrix <- makeContrasts(
  vs1  = TS18h.vehicle_control-TS2h.vehicle_control, # 对照组在前后的差异表达基因
  vs2  = TS18h.PE2.0_embolized-TS2h.PE2.0_embolized, # PE2.0处理前后的差异基因
  vs3  = TS18h.PE1.3_embolized-TS2h.PE1.3_embolized, # PE1.3在处理前后差异基因
  # 处理18小时候，PE2.0相对于对照变化的基因再与PE1.3与对照的差异比较
  diff = (TS18h.PE2.0_embolized-TS18h.vehicle_control)-(TS18h.PE1.3_embolized-TS18h.vehicle_control),
  levels = design
)

fit2 <- contrasts.fit(fit, cont.matrix)
results <- decideTests(fit2)

vennDiagram(results)
###################################