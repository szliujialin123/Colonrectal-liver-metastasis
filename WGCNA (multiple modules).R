rm(list = ls())
library(readxl)
library(WGCNA)
library(flashClust)
library(reshape2)
library(stringr)
enableWGCNAThreads()

# setwd(dir = "D:/Projects/Bioinformatics/ZColorectal_liver_metastasis/Results/Preprocessing/")
# DataExpr <- read.csv(file = "ColonliverNormalized.csv", row.names = 1, check.names = F)
# Trait <- read.csv(file = "ColonliverAlltrait.csv", stringsAsFactors = F)
# 
# Trait <- Trait[Trait$Status == "M",]
# DataExpr <- DataExpr[,match(Trait$Accession, colnames(DataExpr))]
# identical(colnames(DataExpr), Trait$Accession)
# 
# setwd(dir = "D:/Projects/Bioinformatics/ZColorectal_liver_metastasis/Results/WGCNA/")
# 
# # ################### The following is an alternative step###################
# # #Screen for top 75% genes with Median
# # #absolute deviation > 0.01
# # m.mad <- apply(DataExpr, 1, mad)
# # DataExprVar <- DataExpr[which(m.mad > max(quantile(m.mad, probs = seq(0, 1, 0.25))[2],0.01)),]
# # ###########################################################################
# # #If the previous step was executed
# # #DataExprVar, instead of DataExpr, should be fed
# # #into the following function
# 
# DataExpr <- as.data.frame(t(DataExpr))
# 
# dim(DataExpr)
# head(DataExpr)[,1:8]
# 
# sampleTree <- hclust(dist(DataExpr), method = "average")
# par(cex = 0.6)
# par(mar = c(0,4,2,0))
# plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, labels = F)
# cutHeight <- 250
# # Plot a line to show the cut
# abline(h = cutHeight, col = "red");
# # Determine cluster under the line
# clust = cutreeStatic(sampleTree, cutHeight = cutHeight, minSize = 10)
# table(clust)
# # clust 1 contains the samples we want to keep.
# keepSamples = (clust == 1)
# 
# ################################################################################
# NumericTrait <- Trait[,c(1,3)]
# 
# NumericTrait$outlier <- 1
# NumericTrait$outlier[keepSamples] <- 0
# 
# Row <- Trait$Accession
# Col <- colnames(NumericTrait)
# 
# NumericTrait <- as.data.frame(NumericTrait[,3])
# 
# rownames(NumericTrait) <- Row
# colnames(NumericTrait) <- Col[3]
# 
# traitColors <- as.data.frame(numbers2colors(NumericTrait, signed = FALSE))
# 
# rownames(traitColors) <- Row
# colnames(traitColors) <- Col[3]
# 
# svg(file = paste("Detection of outliers",".svg", sep = ""), height = 4, width = 6)
# plotDendroAndColors(sampleTree, traitColors, dendroLabels = F, abHeight = cutHeight,
#                     main = "Sample dendrogram and trait heatmap")
# dev.off()
# 
# DataExpr <- DataExpr[NumericTrait[,"outlier"] == 0,]
# 
# #Pick a threshold to build sale-free-network
# networktype <- "signed"
# powers <- c(c(1:10), seq(from = 12, to = 30, by = 2))
# sft <- pickSoftThreshold(DataExpr, powerVector = powers,
#                          networkType = networktype, verbose = 5)
# 
# svg(file = paste("Scale independence",".svg", sep = ""), height = 10, width = 5)
# cex1 = 0.9
# plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      xlab = "Soft Threshold (power)",
#      ylab = "Scale Free Topology Model Fit, signed R^2",type = "n",
#      main = paste("Scale independence"))
# text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      labels = powers,cex = cex1,col = "red")
# abline(h = 0.85, col = "red")
# dev.off()
# 
# svg(file = paste("Mean connectivity",".svg", sep = ""), height = 10, width = 5)
# plot(sft$fitIndices[,1], sft$fitIndices[,5],
#      xlab = "Soft Threshold (power)",ylab = "Mean Connectivity", type = "n",
#      main = paste("Mean connectivity"))
# text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers,
#      cex = cex1, col = "red")
# dev.off()
# 
# power <- sft$powerEstimate
# 
# ###############################Dynamic Cut######################################
# ############################Calculated in Linux#################################
# adjacency <- adjacency(DataExpr, power = power, type = networktype)
# TOM <- TOMsimilarity(adjacency, TOMType = networktype)
# save(TOM, file = "TOM.Rdata")
# dissTOM <- 1-TOM
# 
# hierTOM <- hclust(as.dist(dissTOM), method = "average")
# 
# sizeGrWindow(12,9)
# plot(hierTOM, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
# 
# dynamicMods <- cutreeDynamic(dendro = hierTOM, distM = dissTOM,
#                             deepSplit = 2, pamRespectsDendro = FALSE,
#                             minClusterSize = 30)
# table(dynamicMods)
# 
# dynamicColors <- labels2colors(dynamicMods)
# table(dynamicColors)
# # Plot the dendrogram and colors underneath
# svg(file = paste("Mean connectivity Cluster Dendrogram and Color (Dynamic)",".svg", sep = ""), height = 4, width = 6)
# 
# plotDendroAndColors(hierTOM, dynamicColors, "Dynamic Tree Cut",
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05,
#                     main = "Gene dendrogram and module colors")
# 
# dev.off()
# 
# # Calculate eigengenes
# MEList <- moduleEigengenes(DataExpr, colors = dynamicColors)
# MEs <- MEList$eigengenes
# # Calculate dissimilarity of module eigengenes
# MEDiss <- 1-cor(MEs)
# # Cluster module eigengenes
# METree <- hclust(as.dist(MEDiss), method = "average")
# 
# #Calculate MEs with merged color labels
# MEList <- moduleEigengenes(DataExpr, colors = dynamicColors)
# MEs <- MEList$eigengenes
# 
# GenenamesColor <- colnames(DataExpr)
# names(GenenamesColor) <- dynamicColors
# 
# save(GenenamesColor,file = "dynamicColors.Rdata")
# save(dynamicColors, DataExpr, NumericTrait, hierTOM, MEs, file = "WGCNA_ResultsDynamic.Rdata")

################################################################################
rm(list = ls())
library(WGCNA)

CIBERSORT <- read.csv(file = "D:/Projects/Bioinformatics/ZColorectal_liver_metastasis/Results/CIBERSORT/CIBERSORT Colon Liver.csv", check.names = F)
CIBERSORT <- CIBERSORT[,which(colnames(CIBERSORT) == "Input Sample"):which(colnames(CIBERSORT) == "Neutrophils")]

setwd(dir = "D:/Projects/Bioinformatics/ZColorectal_liver_metastasis/Results/WGCNA/")

load(file = "D:/Projects/Bioinformatics/ZColorectal_liver_metastasis/Results/WGCNA/dynamicColors.Rdata" )
load(file = "D:/Projects/Bioinformatics/ZColorectal_liver_metastasis/Results/WGCNA/WGCNA_ResultsDynamic.Rdata" )

# Calculate the dissimilarity between merged module eigengene
MEDiss = 1-cor(MEs)
# Cluster merged module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")

svg(file = paste("Merged Eigengene adjacency heatmap",".svg", sep = ""), height = 4, width = 8)
plotEigengeneNetworks(MEs, "Merged eigengene adjacency heatmap",
                      marHeatmap = c(3,4,2,2), plotDendrograms = F, xLabelsAngle = 90)
dev.off()

#The larger MEDissThres lead to less modules (more similar modules will be merged)
MEDissThres <- 0.3

# Plot the result
png(file = paste("Merge similar modules",".png"), height = 1200, width = 2000, res = 300)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
# Plot the cut line into the dendrogram
abline(h = MEDissThres, col = "red")
dev.off()

# Call an automatic merging function
merge <- mergeCloseModules(DataExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors <- merge$colors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

png(file = paste("Merged color and dendrogram",".png"), height = 1200, width = 2000, res = 300)
plotDendroAndColors(hierTOM, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#Calculate MEs with merged color labels
MEList <- moduleEigengenes(DataExpr, colors = mergedColors)
MEs <- MEList$eigengenes

CIBERSORT <- CIBERSORT[match(rownames(MEs), CIBERSORT$`Input Sample`),]
identical(CIBERSORT$`Input Sample`, rownames(MEs))
rownames(CIBERSORT) <- CIBERSORT$`Input Sample`
CIBERSORT <- CIBERSORT[,-1]

moduleTraitCor <- cor(MEs, CIBERSORT, use = "p")

signifTraitCor <- t(signif(moduleTraitCor, 2))
signifTraitCor[signifTraitCor >= 0.3] <- "\U2731"

moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(MEs))
signifmoduleTraitPvalue <- t(signif(moduleTraitPvalue, 1))
textMatrix <- paste(signifTraitCor, "\n(", signifmoduleTraitPvalue, ")", sep = "")

#heavy asterisk, as well as other unicode symbols, are not supported by default svg() function
#therefore, svglite should be used 

svglite::svglite(filename = paste("Module-trait relationships",".svg", sep = ""), height = 12, width = 10)
par(mfrow = c(1,1), mai = c(2,2.5,1,1))
labeledHeatmap(Matrix = t(moduleTraitCor),
               xLabels = names(MEs),
               yLabels = colnames(CIBERSORT),
               xSymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.65,
               zlim = c(-1,1),
               colors.lab.x = "black",
               colors.lab.y = "black",
               xLabelsAngle = 45,
               main = paste("Module-trait relationships"))
par(mfrow = c(1,1), mai = c(1,1,1,1))
dev.off()

#Calculate correlations between genes and modules (module membership, MM)
nSamples <- nrow(DataExpr)
modNames <- substring(names(MEs), 3)
geneModuleMembership <- as.data.frame(cor(DataExpr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "")

#Calculate correlations between genes and trait (gene significance, GS)
#And associate MM with GS

ScatterplotMMGS <- function(status, module, Genesignificance = NULL, Modulemembership = NULL){
        Phenotype <- as.data.frame(CIBERSORT[,which(colnames(CIBERSORT) == status)])
        names(Phenotype) <- status
        geneTraitSignificance <- as.data.frame(cor(DataExpr, Phenotype, use = "p"))
        GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

        names(geneTraitSignificance) <- paste("GS.", names(Phenotype), sep = "")
        names(GSPvalue) <- paste("p.GS.", names(Phenotype), sep = "")

        column <- match(module, modNames)
        
        # moduleGenes <- mergedColors == module
        moduleGenes <- dynamicColors == module

        par(mfrow = c(1,1))
        
        svg(file = paste("MM GS scatter plot of", status, "and", module, "module",".svg", sep = ""), height = 7, width = 7)
        verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                           abs(geneTraitSignificance[moduleGenes, 1]),
                           xlab = paste("Module Membership in", module, "module"),
                           ylab = paste("Gene significance for", status),
                           main = paste("Module membership vs. gene significance\n"),
                           cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
        if (!Genesignificance == 0){
                abline(h = Genesignificance, col = module) 
        }
        if (!Modulemembership == 0){
        abline(v = Modulemembership, col = module)
        }
        dev.off()
        
        ModuleHighMM <- intersect(rownames(geneModuleMembership[abs(geneModuleMembership[,grep(colnames(geneModuleMembership), pattern = module)]) > Modulemembership,]),
                                  rownames(geneModuleMembership)[dynamicColors == module])
        ModuleHighGS <- intersect(rownames(geneTraitSignificance)[which((abs(geneTraitSignificance[,1]) > Genesignificance) == TRUE)],
                                  rownames(geneTraitSignificance)[dynamicColors == module])
        
        
        geneModuleMembership[moduleGenes, column]
        par(mfrow = c(1,1), mai = c(1,1,1,1))
        
        return(intersect(ModuleHighMM, ModuleHighGS))
}

Significantpurple <- ScatterplotMMGS(status = "NK cells activated", module = "purple", Genesignificance = 0, Modulemembership = 0)
write.table(Significantpurple,file = "Significantpurple.txt",row.names = F,quote = F,col.names = F)

library(clusterProfiler)
library(org.Hs.eg.db)

EnrichmentAnalysis <- function(x, GO = FALSE, KEGG = FALSE, ont = FALSE, keyType = FALSE, showCategory = 10){
        filename <- deparse(substitute(x))
        if (GO == TRUE & KEGG == FALSE) {
                ego_ALL <- enrichGO(gene = x,
                                    OrgDb = org.Hs.eg.db,
                                    keyType = keyType,
                                    ont = ont,
                                    pAdjustMethod = "BH",
                                    pvalueCutoff = 1,
                                    qvalueCutoff = 1)
                for (i in 1:showCategory){
                        if (nchar(ego_ALL@result$Description[i]) >= 100){
                                CutPoint <- which(charToRaw(ego_ALL@result$Description[i]) == charToRaw(','))[floor(length(which(charToRaw(ego_ALL@result$Description[i]) == charToRaw(',')))/2)+1]
                                ego_ALL@result$Description[i] <- paste(substr(ego_ALL@result$Description[i], 1, CutPoint), "\n", substr(ego_ALL@result$Description[i], CutPoint+1, nchar(ego_ALL@result$Description[i])))
                        }
                }
                return(ego_ALL)
        }
        if (GO == FALSE & KEGG == TRUE & ont == FALSE) {
          EG_IDs <- mget(x, revmap(org.Hs.egSYMBOL), ifnotfound = NA)
          EG_IDs <- EG_IDs[!EG_IDs == "NA"]
          KEGG_input <- as.character(EG_IDs[match(x, table = names(EG_IDs))])
          KEGG_input <- KEGG_input[!KEGG_input == "NULL"]
          ego_ALL <- enrichKEGG(KEGG_input,
                                keyType = keyType,
                                pAdjustMethod = "BH",
                                pvalueCutoff = 1,
                                qvalueCutoff = 1)
          for (i in 1:length(ego_ALL@result$geneID)){
            ego_ALL@result$geneID[i] <- paste(names(EG_IDs)[match(strsplit(ego_ALL@result$geneID[i], split = "/")[[1]], table = EG_IDs)], collapse = "/")
          }
                for (i in 1:showCategory){
                        if (nchar(ego_ALL@result$Description[i]) >= 100){
                                CutPoint <- which(charToRaw(ego_ALL@result$Description[i]) == charToRaw(','))[floor(length(which(charToRaw(ego_ALL@result$Description[i]) == charToRaw(',')))/2)+1]
                                ego_ALL@result$Description[i] <- paste(substr(ego_ALL@result$Description[i], 1, CutPoint), "\n", substr(ego_ALL@result$Description[i], CutPoint+1, nchar(ego_ALL@result$Description[i])))
                        }
                }
                return(ego_ALL)
        }
}

GO <- EnrichmentAnalysis(Significantpurple, GO = T, keyType = "SYMBOL", ont = "all", showCategory = 100)
KEGG <- EnrichmentAnalysis(Significantpurple, KEGG = TRUE, keyType = "kegg", showCategory = 100)

MergeGOKEGG <- function(x, y){
        x <- x@result
        y <- y@result
        y <- cbind(rep("KEGG", nrow(y)), y)
        colnames(y)[1] <- "ONTOLOGY"
        All <- rbind(x, y)
        All <- dplyr::arrange(All, All$qvalue)
        All <- All[c(which(All$ONTOLOGY == "BP")[1:30],
                     which(All$ONTOLOGY == "MF")[1:30],
                     which(All$ONTOLOGY == "CC")[1:30],
                     which(All$ONTOLOGY == "KEGG")[1:30]),]
}

All <- MergeGOKEGG(GO, KEGG)

library(ggplot2)
library(dplyr)

pick <- All[grep(All$Description, pattern = "cancer"),]

pick <- arrange(pick, Count)

pick$Description <- factor(pick$Description, levels = pick$Description)
TotalGene <- as.numeric(strsplit(pick$GeneRatio, split = "/")[[1]][2])
pick$GeneRatio <- pick$Count/TotalGene

p <- ggplot(pick) + geom_point(aes(x = GeneRatio, y = Description, size = Count, color = -log10(qvalue))) + scale_color_gradient(low = "#E3BDFB",high = 'purple')


pr <- p #+ scale_color_gradient(low = "#ffffff", high = "#A020F0")

pr <- pr + labs(color = expression(-log[10](qvalue)),size = "Count",
                x = "GeneRatio", title = paste("Pathway enrichment", sep = " "))

png(file = "cancer_enrichment.png", height = 1200, width = 1400, res = 350)
print(pr + theme(axis.text.x = element_blank(), axis.title.y.left = element_blank()
                 , plot.background = element_rect(fill = "white")))
dev.off()
 

