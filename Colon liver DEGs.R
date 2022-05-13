rm(list = ls())

setwd(dir = "D:/Projects/Bioinformatics/ZColorectal_liver_metastasis/Results/Preprocessing/")

Trait <- read.csv(file = "ColonliverAlltrait.csv", stringsAsFactors = F)
DataExpr <- read.csv(file = "Mergedcleanmatrix_original_scale.csv", row.names = 1, check.names = F)
DataExpr <- log2(DataExpr)

setwd(dir = "D:/Projects/Bioinformatics/ZColorectal_liver_metastasis/Results/DEG/") 

library(limma)

DEGanalysis <- function(x){
  
  x <- DataExpr
  
  status <- Trait$Status[match(colnames(DataExpr), table = Trait$Accession)]
  status <- as.factor(status)
  batch <- as.factor(Trait$Batch)
  
  ###The following step is crucial! The level of status (whether M priors to P)
  #determines how different traits were compared in limma, default level follows 
  #an alphabetic order, therefore, "M" ranks above "P", and representing a 
  #comparison of "M" against "P", the resulting logFC and up-regulated genes
  #were up-regulated in "M" as compared to "P". In other circumstances when
  #the intended comparison was "P" against "M", the levels should be reversed
  #levels(status) <- rev(levels(status))
  ##############################
  design <- model.matrix(~0+status+batch)
  
  png(filename = paste(deparse(substitute(x)), "before Voom.png"), res = 300, width = 3000, height = 3300)
  print(boxplot(x))
  dev.off()

  png(filename = paste(deparse(substitute(x)), "Voom.png"), res = 300, width = 3000, height = 3300)
  print(v <- voom((x+abs(min(x))), design, normalize = 'quantile', plot = TRUE))
  dev.off()

  png(filename = paste(deparse(substitute(x)), "after Voom.png"), res = 300, width = 3000, height = 3300)
  print(boxplot(as.data.frame(v)))
  dev.off()
  
  # fit <- lmFit(v, design)
  fit <- lmFit(x, design)
  status <- as.factor(Trait$Status[match(colnames(x), table = Trait$Accession)])
  cont.matrix <- makeContrasts("statusM-statusP", levels = design)
  fit2 <- contrasts.fit(fit,cont.matrix)
  fit2 <- eBayes(fit2)  ## default no trend !!!
  ##eBayes() with trend = TRUE
  
  tempOutput <- topTable(fit2, coef = 1, n = Inf, adjust = "BH")
  nrDEG <- na.omit(tempOutput)
  return(nrDEG)
}

MvsP <- DEGanalysis(DataExpr)
write.csv(MvsP, file = "MvsP.csv")

library(ggplot2)
library(ggrepel)

Volcanoplot <- function(x, sig = NULL, FC = NULL, ylimit = NULL, xlimit = NULL){
  x$significance[(x$adj.P.Val > sig|x$adj.P.Val=="NA")|(x$logFC < FC)& x$logFC > -FC] <- "no"
  x$significance[x$adj.P.Val <= sig & x$logFC >= FC] <- "up"
  x$significance[x$adj.P.Val <= sig & x$logFC <= -FC] <- "down"
  x$ID <- rownames(x)
  theme_set(theme_bw())
  p <- ggplot(x,aes(logFC, -1*log10(adj.P.Val),
                    color = significance)) + geom_point(alpha = 0.3, size = 1) +
    xlim(-xlimit,xlimit) + ylim(0,ylimit) + labs(x = "LogFoldChange",y = "-Log10(q-value)")
  p <- p + scale_color_manual(values = c("#0072B5","grey","#BC3C28"))+
    geom_hline(yintercept = -log10(sig), linetype = 5)+
    geom_vline(xintercept = c(-FC, FC),linetype = 5)
  p <- p + theme(panel.grid = element_blank())+
    theme(axis.line = element_line(size = 0))
  p <- p + guides(colour = FALSE)
  p <- p + theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20))
  
  print(subset(x, x$adj.P.Val <= sig & abs(x$logFC) >= FC))
  DEG <- subset(x, x$adj.P.Val <= sig & abs(x$logFC) >= FC)
  
  show <- c(rownames(dplyr::arrange((x),desc(logFC)))[1:5], rownames(dplyr::arrange((x),logFC))[1:5])
  p <- p + geom_text_repel(data = subset(x, rownames(x) %in% show & rownames(x) %in% DEG$ID), aes(label = ID))
  print(p)
  return(DEG)
}

svg(filename = "MvsP.svg", width = 5, height = 5)
MvsPVol <- Volcanoplot(MvsP, sig = 0.05, FC = 0.2, ylimit = 350, xlimit = 12)
dev.off()

TestM <- DataExpr[rownames(DataExpr) == "RBM10", Trait$Status == "M"]
TestP <- DataExpr[rownames(DataExpr) == "RBM10", Trait$Status == "P"]

mean(as.numeric(TestM))
mean(as.numeric(TestP))

TestM <- DataExpr[rownames(DataExpr) == "TUBA1C", Trait$Status == "M"]
TestP <- DataExpr[rownames(DataExpr) == "TUBA1C", Trait$Status == "P"]

mean(as.numeric(TestM))
mean(as.numeric(TestP))

write.csv(MvsPVol, file = "MvsPVol.csv")