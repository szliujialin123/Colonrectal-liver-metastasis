##################Batch effect detection and data integration###################
# library(dplyr)
# library(Seurat)
# library(harmony)
# 
# rm(list = ls())
# 
# setwd(dir = "D:/Projects/Bioinformatics/ZColorectal_liver_metastasis/GSE178318/Analyses/")
# load(file = "D:/Projects/Bioinformatics/ZColorectal_liver_metastasis/GSE178318/Preprocessing/Preprocessed_scRNA.Rdata")
# 
# All.data <- Reduce(merge, CRC.list)
# 
# All.data <- All.data %>% NormalizeData() %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
#   ScaleData() %>% RunPCA() %>% RunUMAP(reduction = "pca", dims = 1:50)
# 
# png(filename = "Seurat Before batch correction.png", width = 2500, height = 2000, res = 300)
# DimPlot(All.data, reduction = "umap", group.by = "orig.ident")
# dev.off()
# 
# All.data <- RunHarmony(All.data, group.by.vars = "orig.ident", plot_convergence = T)
# 
# All.data <- RunUMAP(All.data, reduction = "harmony", dims = 1:50)
# 
# png(filename = "Seurat after batch correction.png", width = 2500, height = 2000, res = 300)
# DimPlot(All.data, reduction = "umap", group.by = "orig.ident")
# dev.off()

##############################Cluster annotation################################
setwd(dir = "D:/Projects/Bioinformatics/ZColorectal_liver_metastasis/GSE178318/Analyses/")

# png(filename = "ElbowPlot.png", width = 2500, height = 2000, res = 500)
# ElbowPlot(All.data, ndims = 50)
# dev.off()
# 
# All.data <- FindNeighbors(All.data, reduction = "pca", dims = 1:30)
# All.data <- FindClusters(All.data)
# All.data <- RunUMAP(All.data, reduction = "harmony", dims = 1:30)
# 
# All.data$Group <- ""
# All.data$Group[grep(All.data$orig.ident, pattern = "CRC")] <- "CRC"
# All.data$Group[grep(All.data$orig.ident, pattern = "LM")] <- "LM"
# 
# save(All.data, file = "scRNA_batch_corrected.Rdata")

rm(list = ls())

library(dplyr)
library(Seurat)
library(ggplot2)
library(aggregateBioVar)

load(file =  "D:/Projects/Bioinformatics/ZColorectal_liver_metastasis/GSE178318/Analyses/scRNA_batch_corrected.Rdata")

png(filename = "Cluster umap.png", width = 2000, height = 2000, res = 300)
DimPlot(All.data, reduction = "umap", label = T) + NoLegend()
dev.off()

png(filename = "FeaturePlotNK.png", width = 1000, height = 1000, res = 300)
FeaturePlot(All.data, features = c("KLRF1"))
dev.off()

png(filename = "FeaturePlotMAST.png", width = 1000, height = 1000, res = 300)
FeaturePlot(All.data, features = c("TPSAB1"))
dev.off()

png(filename = "FeaturePlotMASTfeature.png", width = 2000, height = 1000, res = 300)
FeaturePlot(All.data[,All.data$seurat_clusters == "20"], features = c("TPT1"), split.by = "Group")
dev.off()

Reanno.cluster.ids <- levels(All.data@active.ident)

Reanno.cluster.ids[c(0, 1, 2, 3, 4, 5, 7, 12, 15, 16)+1] <- "T cells"
Reanno.cluster.ids[c(9)+1] <- "NK cells"
Reanno.cluster.ids[c(4, 6, 14, 21, 24)+1] <- "Myeloid cells"
Reanno.cluster.ids[c(20)+1] <- "Mast cells"
Reanno.cluster.ids[c(8)+1] <- "B cells"
Reanno.cluster.ids[c(21)+1] <- "pDC"
Reanno.cluster.ids[c(19, 22)+1] <- "CAF"
Reanno.cluster.ids[c(23)+1] <- "Endothelial cells"
Reanno.cluster.ids[c(10, 13 ,17, 18)+1] <- "Cancer cells" #EPCAM
Reanno.cluster.ids[c(11)+1] <- "Plasma cells" #JCHAIN

names(Reanno.cluster.ids) <- levels(All.data@active.ident)
All.data <- RenameIdents(All.data, Reanno.cluster.ids)

png(filename = "Main cell dim.png", width = 1500, height = 1500, res = 300)
DimPlot(All.data, label = T) + theme(legend.position = "none")
dev.off()

svg(filename = "Main cell markers.svg", width = 8, height = 6)
DotPlot(All.data, features = c("CD3E", "KLRF1", "LYZ", "CD79A", "TPSAB1", "JCHAIN", "LILRA4", "COL1A2", "PECAM1", "EPCAM")) + RotatedAxis() +
  scale_x_discrete("")+scale_y_discrete("")
dev.off()

pseudobulk <- countsBySubject(scExp = as.SingleCellExperiment(All.data), subjectVar = "orig.ident")
pseudobulk <- apply(as.data.frame(pseudobulk), 2, FUN = function(x) {
  return(log1p(x/sum(x)*1000000))
})

AVECRC <- as.data.frame(apply(pseudobulk[,grep(colnames(pseudobulk), pattern = "CRC")], 1, mean))
AVELM <- as.data.frame(apply(pseudobulk[,grep(colnames(pseudobulk), pattern = "LM")], 1, mean))

# Example of ~
# wilcox.test(Ozone ~ Month, data = airquality,
#             subset = Month %in% c(5, 8), exact = F)
# iris_long <- iris %>%
#   dplyr::mutate(id = dplyr::row_number(x = Species)) %>%
#   tidyr::gather(
#     data = .,
#     key = "condition",
#     value = "value",
#     Sepal.Length:Petal.Width,
#     convert = TRUE,
#     factor_key = TRUE
#   ) %>%
#   tidyr::separate(
#     col = "condition",
#     into = c("part", "measure"),
#     sep = "\\.",
#     convert = TRUE
#   ) %>%
#   tibble::as_data_frame(x = .)
# 
# 
# iris_long$value
# table(iris_long$part)
# 
# groupedstats::grouped_wilcox(
#   data = pseudobulk,
#   dep.vars = Gene, # dependent variable
#   indep.vars = Value, # independent variable
#   grouping.vars = c(Species, measure), # for each Species and for each measurement
#   paired = TRUE # paired Wilcoxon signed rank test with continuity correction
# )

pseudobulktest <- reshape2::melt(pseudobulk)
colnames(pseudobulktest) <- c("Gene", "Group", "Value")
pseudobulktest$Group <- substring(pseudobulktest$Group, first = 7, last = 9)

WilcoxTEST <- groupedstats::grouped_wilcox(
  data = as_tibble(pseudobulktest),
  dep.vars = Value,
  indep.vars = Group,
  grouping.vars = Gene,
  paired = T, correct = T
)

WilcoxTEST <- cbind((AVELM - AVECRC), WilcoxTEST)
 
# 
# Wilcox <- function(x,y){
#   statistic <- c()
#   pvalue <- c()
#   for (i in 1:nrow(x)){
#     X <- as.numeric(x[i,])
#     Y <- as.numeric(y[i,])
#     wilcox1 <- wilcox.test(X, Y, alternative = "two.sided", exact = F, paired = T)
#     statistic <- rbind(statistic, wilcox1$statistic)
#     pvalue <- rbind(pvalue, wilcox1$p.value)
#   }
#   tmp <- cbind(statistic, pvalue)
#   rownames(tmp) <- rownames(x)
#   colnames(tmp) <- c("statistic", "pvalue")
#   tmp <- as.data.frame(tmp)
#   return(tmp)
# }
# WilcoxTEST1 <- cbind((AVECRC - AVELM), Wilcox(pseudobulk[,grep(colnames(pseudobulk), pattern = "CRC")],
#                       pseudobulk[,grep(colnames(pseudobulk), pattern = "LM")]))
# i = length(WilcoxTEST$pvalue):1  # The reverse rank order
# o <- order(WilcoxTEST$pvalue, decreasing = TRUE)
# ro <- order(o)
# WilcoxTEST$p.adjusted <- pmin(1, cummin(length(WilcoxTEST$pvalue)/i * WilcoxTEST$pvalue[o]))[ro]

colnames(WilcoxTEST)[1] <- "FoldChange"

BulkDiff <- WilcoxTEST[abs(WilcoxTEST$FoldChange) > 0.2 & WilcoxTEST$p.value < 0.05,]

DEG <- read.csv(file = "D:/Projects/Bioinformatics/ZColorectal_liver_metastasis/Results/DEG/MvsP.csv", row.names = 1)

Significantpurple <- read.table(file = "D:/Projects/Bioinformatics/ZColorectal_liver_metastasis/Results/WGCNA/Significantpurple.txt")$V1

ConsensusUP <- intersect(rownames(DEG)[DEG$logFC > 0], rownames(BulkDiff)[BulkDiff$FoldChange > 0])
ConsensusDown <- intersect(rownames(DEG)[DEG$logFC < 0], rownames(BulkDiff)[BulkDiff$FoldChange < 0])
  
ImmuneSignature <- intersect(c(ConsensusUP, ConsensusDown), Significantpurple)

display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, disable.logging = TRUE, ...)
  grid.draw(venn_object)
}

png(filename = "Up-regulated.png", width = 1000, height = 1000, res = 300)
display_venn(x= list('GEO Up' = rownames(DEG)[DEG$logFC > 0], 'scRNA Up' = rownames(BulkDiff)[BulkDiff$FoldChange > 0]), 
             col ="transparent", fill = c("pink","orange"), alpha = 0.5, 
             cex = 1,fontfamily = "serif", fontface = "bold", cat.cex = 0.9, cat.pos = 0, cat.dist = 0.03, cat.fontfamily = "serif", 
             rotation.degree = 360)
dev.off()

png(filename = "Down-regulated.png", width = 1000, height = 1000, res = 300)
display_venn(x= list('GEO Down' = rownames(DEG)[DEG$logFC < 0], 'scRNA Down' = rownames(BulkDiff)[BulkDiff$FoldChange < 0]), 
             col ="transparent", fill = c("turquoise","yellowgreen"), alpha = 0.5, 
             cex = 1,fontfamily = "serif", fontface = "bold", cat.cex = 0.9, cat.pos = 0, cat.dist = 0.03, cat.fontfamily = "serif", 
             rotation.degree = 360)
dev.off()

png(filename = "ImmuneSignatures.png", width = 1000, height = 1000, res = 300)
display_venn(x= list('CommonUp' = ConsensusUP, 
                     'CommonDown' = ConsensusDown, 
                     'PurpleModule' = Significantpurple),
             col ="transparent", fill = c("#FFC373", "#9DDF8D", "purple"), alpha = 0.5, 
             cex = 1,fontfamily = "serif", fontface = "bold", cat.cex = 1, cat.pos = 0, cat.dist = 0.1, cat.fontfamily = "serif", 
             rotation.degree = 360)
dev.off()

Trait <- read.csv(file = "D:/Projects/Bioinformatics/ZColorectal_liver_metastasis/Results/Preprocessing/ColonliverAlltrait.csv", stringsAsFactors = F)
DataExpr <- read.csv(file = "D:/Projects/Bioinformatics/ZColorectal_liver_metastasis/Results/Preprocessing/ColonliverNormalized.csv", row.names = 1, check.names = F)

setwd(dir = "D:/Projects/Bioinformatics/ZColorectal_liver_metastasis/Results/LASSO/")

library(glmnet)
library(ggplot2)

x <- t(DataExpr[rownames(DataExpr) %in% ImmuneSignature,])
y <- Trait$Status[match(Trait$Accession, colnames(DataExpr))]
y[y == "M"] <- 1
y[y == "P"] <- 0
y <- as.numeric(y)

fit <- glmnet(x, y, family = 'binomial')

set.seed(1)

png(filename = "Coefficient.png", width = 1500, height = 1500, res = 300)
plot(fit, lwd = 2, xvar = 'lambda')
dev.off()

png(filename = "CV.png", width = 1500, height = 1500, res = 300)
cv <- cv.glmnet(x, y, family = 'binomial')
plot(cv)
dev.off()

coefficients <- as.matrix(coef(cv, s = cv$lambda.min))
coeffi <- as.data.frame(sort(coefficients[coefficients[, 1]!=0, ][-1], decreasing = TRUE))
colnames(coeffi) <-'weight'

coef <- rownames(coeffi)

print(coef)

png(filename = " coeffi_weight.png", width = 1500, height = 1500, res = 300)
weight <- coeffi$weight
ggplot(coeffi, aes(coef, weight)) + coord_flip() + 
  ggtitle('Coefficient Weights') + 
  geom_bar(aes(fill = factor((weight < 0) + 1)), stat="identity") + 
  theme(plot.title = element_text(hjust = 0.4, face='bold', size = 20, family = 'serif'), 
        axis.text.x = element_text(hjust = 0.4, face = 'bold', size = 20, colour = 'black', family = 'serif'), 
        axis.text.y = element_text(hjust = 0.5, face = 'bold', size = 16, colour = 'black', family = 'serif'), 
        axis.title.x = element_text(face = 'bold', size = 20, family = 'serif'), 
        axis.title.y = element_text(face = 'bold', size=20, family = 'serif'), 
        legend.position = 'none')
dev.off()

train <- DataExpr[rownames(coeffi), ]
colnames(train) <- paste(Trait$Status, colnames(train), sep = "")
write.csv(train, "D:/Projects/Bioinformatics/ZColorectal_liver_metastasis/Results/MachineLearning/train.csv")

test <- pseudobulk[rownames(coeffi),]
write.csv(test, "D:/Projects/Bioinformatics/ZColorectal_liver_metastasis/Results/MachineLearning/test.csv")

################################################################################

setwd(dir = "D:/Projects/Bioinformatics/ZColorectal_liver_metastasis/Results/scRNADownStream/")

library(ComplexHeatmap)
library(circlize)

svg(filename = "HeatMap.svg", width = 7, height = 5)
scaled <- t(apply(test, 1, scale))
colnames(scaled) <- substring(colnames(test),7,9)
Status <- colnames(scaled)
set.seed(3)

draw(Heatmap(scaled, cluster_rows = FALSE, row_dend_reorder = TRUE,
        row_names_gp = gpar(fontsize = 10), show_column_names = FALSE,
        column_title = NULL, top_annotation = 
          HeatmapAnnotation(Malignancy = Status, show_annotation_name = F,  col = 
                              list(Malignancy = c("CRC" = "#FFCC66", "LM" = "#CC00CC")), 
                            annotation_legend_param = list(Malignancy = list(at = c("CRC", "LM")), nrow = 1)),
        heatmap_legend_param = list(
          title = "PseudoBulk.Scaled.Exp", at = c(-2, 2),direction = "horizontal",
          labels = c("-2", "2"), legend_width = unit(3, "cm")), 
        column_km = 2),  heatmap_legend_side = "bottom", annotation_legend_side = "top")
dev.off()

################################################################################

for (gene in rownames(coeffi)) {
  png(filename = paste(gene, "split.png"), width = 3000, height = 1500, res = 300)
  print(FeaturePlot(All.data, features = gene, split.by = "Group", label = T))
  dev.off()
}

################################ReSetIdentity###################################

for (i in levels(All.data@active.ident)) {
  for (Malignancy in c("LM", "CRC")) {
    cells.use <- WhichCells(All.data[,All.data$Group == Malignancy], idents = i)
    All.data <- SetIdent(All.data, cells = cells.use, value = paste(Malignancy, i))
  }
}

################################################################################

DiffMastAll <- FindMarkers(All.data, ident.1 = "LM Mast cells", ident.2 = "CRC Mast cells", logfc.threshold = 0)
DiffNKAll <- FindMarkers(All.data, ident.1 = "LM NK cells", ident.2 = "CRC NK cells", logfc.threshold = 0)

# Genome-wide gene expression profiling of human mast cells stimulated by IgE or Fc¦ÅRI-aggregation reveals a complex network of genes involved in inflammatory responses
DiffMast <- DiffMastAll[rownames(DiffMastAll) %in% c("CCL4", "CXCL1", "CCL7",
                                   "CXCL3", "CCL5", "CCRL2",
                                   "IL8", "BL34", "LIF", "BL34",
                                   "LIF", "CD69", "LAT", "CD83",
                                   "ADORA2A", "CRIP1", "HLA-DQB1",
                                   "HKE2", "TRAF1", "NKG7", "PTX3",
                                   "FCGR2B", "TNFAIP6", "PTGER2", "CLECSF5",
                                   "EBI2", "TREM1", "BIRC3", "FCAR", "TLR2",
                                   "HIVEP1", "FOSB", "EGR3", "PHLDA2", "SERPINB2",
                                   "FGFR1", "EGR2", "INSIG1", "PDGFA", "IER3",
                                   "PDGFB", "BTG2", "TIEG", "CNK", "PBEF1",
                                   "TNFAIP3", "TNFAIP8", "MLP",
                                   "CRABP2", "IRF2", "RASAL1",
                                   "FLRT2", "SMARCD3", "ARHE",
                                   "KAL1", "FLNB", "CD151", "ARF6",
                                   "ALCAM", "MAFF", "MSC",
                                   "NFKB1", "BCL6", "NFATC1", "ATF3",
                                   "NFKBIA", "NFKBIE", "MYC", "EMP1",
                                   "ELL2", "TOP2A", "GTF2H2", "DUSP2",
                                   "THBD", "NR4A2", "GEM", "OLR1",
                                   "GCH1", "SPHK1", "NDUFA7",
                                   "HBEGF", "JAG1", "LDLR",
                                   "MADH7", "MALT1", "SPRY2",
                                   "DUSP6", "MAP2K3", "CREM",
                                   "DUSP1", "MAP3K14", "FUT4",
                                   "JUN", "PSCDBP", "FYN",
                                   "PGGT1B", "VRK2", "TTK", "ENC1",
                                   "SLC16A6", "HIST2H2AA", "CYP3A4",
                                   "HIST1H1C", "HEC", "STK17A",
                                   "PELI1", "KCNAB1", "B4GALT4"                                  
) & DiffMastAll$p_val < 0.1,]

#DUSP1 PMID: 23190643

DiffNK <- DiffNKAll[rownames(DiffNKAll) %in% c("GZMA", "GZMB", "GZMK","GZMH", "IFNG", "KLRK1", "PRF1",
                               "CD69", "TNFRSF18", "KLRF1", "LAMP1", "NKG2D", "NCR1", "NKG2") 
       & DiffNKAll$p_val < 0.1,]

svg(filename = "GroupSplit 8 key genes and NK mast.svg", width = 12, height = 8)
DotPlot(All.data, features = unique(c(rownames(coeffi), 
                                      rownames(DiffNK), 
                                      rownames(DiffMast)))) + 
  RotatedAxis() +
  scale_x_discrete("")+scale_y_discrete("")
dev.off()

LMpseudobulk <- pseudobulk[,grep(colnames(pseudobulk), pattern = "LM")]

res.pca <- FactoMineR::PCA(t(LMpseudobulk[rownames(coeffi),]), graph = FALSE)

library(ggpubr)

df <- rbind(res.pca[["ind"]][["coord"]][,1],
            LMpseudobulk[rownames(DiffNK),],
            LMpseudobulk[rownames(DiffMast),]
      )

rownames(df)[1] <- "EigenGene"

Correlation <- cor(t(df))[1,]

for (gene in names(Correlation)[-1]) {
  df <- t(rbind(res.pca[["ind"]][["coord"]][,1], LMpseudobulk[gene,]))
  colnames(df) <- c("Eigengene", gene)
  df <- data.frame(df)

  r <- cor.test(df[,1], df[,2])
  
  if (r[["p.value"]] < 0.1) {
    svg(filename = paste("Correlation/", gene, ".svg", sep = ""), width = 4, height = 4)
    print(ggscatter(df, x = "Eigengene", y = gene,
                    add = "reg.line",  # Add regressin line
                    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                    conf.int = TRUE # Add confidence interval
    ) + stat_cor(method = "pearson"))
    
    dev.off()
  }
}

DiffNKAll[rownames(DiffNKAll) %in%rownames(coefficients),]
DiffMastAll[rownames(DiffMastAll) %in%rownames(coefficients),]


write.csv(DiffNKAll, file = "DiffNKAll.csv")
write.csv(DiffMastAll, file = "DiffMastAll.csv")
