rm(list = ls())
# library(devtools)
# install_github("jmzeng1314/GEOmirror")
# install_github("jmzeng1314/idmap1")
# 
# library(GEOmirror)
# GSE18549 <- geoChina('GSE18549')
# GSE14095 <- geoChina('GSE14095')
# GSE68468 <- geoChina('GSE68468')
# save(GSE18549, GSE14095, GSE68468, file = "Colon liver datasets.Rdata")

setwd(dir = "D:/Projects/Bioinformatics/ZColorectal_liver_metastasis/Datasets/")

library(idmap1)
load(file = "Colon liver datasets.Rdata")

GSE18549 <- GSE18549$GSE18549_series_matrix.txt.gz
GSE14095 <- GSE14095$GSE14095_series_matrix.txt.gz
GSE68468 <- GSE68468$GSE68468_series_matrix.txt.gz
################################################################################
GSE18549ftr <- GSE18549@phenoData@data
GSE18549ftr <- GSE18549ftr[grep(GSE18549ftr$title, pattern = "Colon_to_Liver"),]
GSE18549ftr$title[grep(GSE18549ftr$source_name_ch1, pattern = "to Liver")] <- "M"
GSE18549ftr <- GSE18549ftr[,c(1,2,8)]
colnames(GSE18549ftr) <- c("Status","Accession","Description")

GSE14095ftr <- readxl::read_xlsx(path = "GSE14095 sample info.xlsx")
GSE14095ftr <- GSE14095ftr[grep(GSE14095ftr$metastasis_status, pattern = "YES"),]
GSE14095ftr$sample_label[grep(GSE14095ftr$metastasis_status, pattern = "YES")] <- "M"
GSE14095ftr <- GSE14095ftr[,c(7,1,6)]
colnames(GSE14095ftr) <- c("Status","Accession","Description")

GSE68468ftr <- GSE68468@phenoData@data
GSE68468ftr <- GSE68468ftr[grep(GSE68468ftr$`histology:ch1`, pattern = "colon carcinoma metastatic to the liver"),]
GSE68468ftr$title[grep(GSE68468ftr$`histology:ch1`, pattern = "colon carcinoma metastatic to the liver")] <- "M"
GSE68468ftr <- GSE68468ftr[,c(1,2,40)]
colnames(GSE68468ftr) <- c("Status","Accession","Description")

GSE18549ftr <- cbind(GSE18549ftr, rep("GSE18549", nrow(GSE18549ftr)))
GSE14095ftr <- cbind(GSE14095ftr, rep("GSE14095", nrow(GSE14095ftr)))
GSE68468ftr <- cbind(GSE68468ftr, rep("GSE68468", nrow(GSE68468ftr)))
colnames(GSE18549ftr)[4] <- "Batch"
colnames(GSE14095ftr)[4] <- "Batch"
colnames(GSE68468ftr)[4] <- "Batch"

Alltrait <- rbind(GSE18549ftr, GSE14095ftr, GSE68468ftr)

GSE18549@experimentData@other$platform_id 
"GPL570"
GSE14095@experimentData@other$platform_id 
"GPL570"
GSE68468@experimentData@other$platform_id 
"GPL96"

GPL570 <- getIDs('GPL570')
GPL96 <- getIDs('GPL96')

GSE18549exp <- as.data.frame(GSE18549@assayData$exprs)
GSE14095exp <- as.data.frame(GSE14095@assayData$exprs)
GSE68468exp <- as.data.frame(GSE68468@assayData$exprs)

GSE14095exp <- log2(GSE14095exp+1)
GSE68468exp <- log2(GSE68468exp+1)

GSE18549exp <- cbind(GPL570$symbol[match(rownames(GSE18549exp), table = GPL570$probe_id)], GSE18549exp)
GSE14095exp <- cbind(GPL570$symbol[match(rownames(GSE14095exp), table = GPL570$probe_id)], GSE14095exp)
GSE68468exp <- cbind(GPL96$symbol[match(rownames(GSE68468exp), table = GPL96$probe_id)], GSE68468exp)

colnames(GSE18549exp)[1] <- "Symbol"  
colnames(GSE14095exp)[1] <- "Symbol"  
colnames(GSE68468exp)[1] <- "Symbol"  

###########################Import data sets#####################################

library(data.table)

setwd(dir = "D:/Projects/Bioinformatics/ZColorectal_liver_metastasis/Datasets/")

GSE131418exp <- fread("GSE131418_MCC_prim_met_GE_probe_level.txt.gz", fill = T, data.table = F, header = T)
GSE131418exp[,c(2:ncol(GSE131418exp))] <- 2^GSE131418exp[,c(2:ncol(GSE131418exp))]
colnames(GSE131418exp)[2:ncol(GSE131418exp)] <- colnames(GSE131418exp)[c(1:ncol(GSE131418exp)-1)]
rownames(GSE131418exp) <- GSE131418exp[,1]
GSE131418exp <- GSE131418exp[,-1]
GSE131418exp <- log2(GSE131418exp)

GSE131418_probe <- fread("GSE131418_cleaned_probes_merck_annotation_Entrez_IDs.csv.gz", fill = T, data.table = F, header = T)

#Some gene symbols have been wrongly converted into date by excel, 
#in GSE131418probe, here we extract probe and gene names to tackle this issue
#1-Dec = DELEC1; 1-Mar = MARCHF1; 1-Sep = SEPTIN1...etc

Excelfault <- read.csv(file = "ExcelfaultProbes.csv", row.names = 1)

tmp <- GSE131418_probe
tmp$V2[which(tmp$V2 %in% Excelfault$Wrong.symbols == TRUE)] <- as.character(Excelfault$Corrected[match(tmp$V2[which(tmp$V2 %in% Excelfault$Wrong.symbols == TRUE)], table = Excelfault$Wrong.symbols)])
GSE131418_probe <- tmp

tmp <- GSE131418_probe$V2[match(rownames(GSE131418exp), table = GSE131418_probe$V1)]
GSE131418exp <- cbind(tmp, GSE131418exp)
colnames(GSE131418exp)[1] <- "Symbol"

Dealwithduplicate <- function(x, Function = NULL){
  for (i in 2:ncol(x)){
    average <- tapply(x[,i], INDEX = x[,"Symbol"], FUN = Function)
    x[,i] <- as.numeric(average[match(x$Symbol, names(average))])
  }
  x
}

GSE131418exp <- Dealwithduplicate(GSE131418exp, Function = max)
GSE14095exp <- Dealwithduplicate(GSE14095exp, Function = max)
GSE18549exp <- Dealwithduplicate(GSE18549exp, Function = max)
GSE68468exp <- Dealwithduplicate(GSE68468exp, Function = max)

GSE131418ftr <- fread(file = "GSE131418_series_matrix.txt.gz", sep = "\t", data.table = F, fill = T)
GSE131418ftr <- GSE131418ftr[-c(1:34),]
rownames(GSE131418ftr) <- c(1:nrow(GSE131418ftr))

GSE131418ftr <- GSE131418ftr[c(1,2,17),]
GSE131418ftr <- as.data.frame(t(GSE131418ftr))

colnames(GSE131418ftr) <- c("Status","Accession","Description")
GSE131418ftr <- GSE131418ftr[-1,]

colnames(GSE131418exp)[as.numeric(na.omit(match(GSE131418ftr$Status, colnames(GSE131418exp))))] <- as.character(GSE131418ftr$Accession[as.numeric(na.omit(match(colnames(GSE131418exp), GSE131418ftr$Status)))])

GSE131418ftr$Status <- as.character(GSE131418ftr$Status)

GSE131418ftr$Status[grep(GSE131418ftr$Description, pattern = "site of metastasis: LIVER")] <- "M"
GSE131418ftr$Status[grep(GSE131418ftr$Description, pattern = "site of metastasis: NA")] <- "P"

GSE131418ftr <- GSE131418ftr[GSE131418ftr$Status %in% c("M","P"),]

GSE131418exp <- cbind(GSE131418exp$Symbol, GSE131418exp[,colnames(GSE131418exp) %in% GSE131418ftr$Accession])
GSE18549exp <- cbind(GSE18549exp$Symbol, GSE18549exp[,colnames(GSE18549exp) %in% Alltrait$Accession])
GSE14095exp <- cbind(GSE14095exp$Symbol, GSE14095exp[,colnames(GSE14095exp) %in% Alltrait$Accession])
GSE68468exp <- cbind(GSE68468exp$Symbol, GSE68468exp[,colnames(GSE68468exp) %in% Alltrait$Accession])

colnames(GSE131418exp)[1] <- "Symbol"
colnames(GSE18549exp)[1] <- "Symbol"
colnames(GSE14095exp)[1] <- "Symbol"
colnames(GSE68468exp)[1] <- "Symbol"

GSE131418ftr <- cbind(GSE131418ftr, rep("GSE131418", nrow(GSE131418ftr)))
colnames(GSE131418ftr)[4] <- "Batch"

GSE131418ftr <- GSE131418ftr[GSE131418ftr$Accession %in% colnames(GSE131418exp),]

Alltrait <- rbind(GSE131418ftr, Alltrait)

tmp1 <- GSE131418exp[!duplicated(GSE131418exp$Symbol),]
tmp2 <- GSE14095exp[!duplicated(GSE14095exp$Symbol),]
tmp3 <- GSE18549exp[!duplicated(GSE18549exp$Symbol),]
tmp4 <- GSE68468exp[!duplicated(GSE68468exp$Symbol),]

intersectiontags <- intersect(intersect(intersect(GSE131418exp$Symbol, GSE14095exp$Symbol), GSE18549exp$Symbol), GSE68468exp$Symbol)
intersectiontags <- intersectiontags[-1]

tmp1 <- tmp1[tmp1$Symbol %in% intersectiontags,]
tmp2 <- tmp2[tmp2$Symbol %in% intersectiontags,]
tmp3 <- tmp3[tmp3$Symbol %in% intersectiontags,]
tmp4 <- tmp4[tmp4$Symbol %in% intersectiontags,]

Mergedcleanmatrix <- Reduce(function(x,y) merge(x = x, y = y, by = "Symbol"), list(tmp1, tmp2, tmp3, tmp4))

rownames(Mergedcleanmatrix) <- Mergedcleanmatrix$Symbol

Mergedcleanmatrix <- Mergedcleanmatrix[,-1]

setwd("D:/Projects/Bioinformatics/ZColorectal_liver_metastasis/Results/Preprocessing/")

Alltrait$Status <- as.character(Alltrait$Status)
Alltrait$Accession <- as.character(Alltrait$Accession)
Alltrait$Description <- as.character(Alltrait$Description)
Alltrait$Batch <- as.character(Alltrait$Batch)

Mergedcleanmatrix <- Mergedcleanmatrix[,colnames(Mergedcleanmatrix) %in% Alltrait$Accession]

library(Rtsne)
library(scales)
#scales for alpha transparency
library(limma)
#remove batch effects

batch <- as.character(Alltrait$Batch[match(colnames(Mergedcleanmatrix), table = Alltrait$Accession)])

dottype <- Alltrait$Status
dottype[which(dottype == "P")] <- 20
dottype[which(dottype == "M")] <- 18
dottype <- as.numeric(dottype)

svg(filename = "Colon liver data before batch correction.svg", height = 5, width = 6)
colors <- rainbow(length(unique(batch)))
names(colors) <- unique(batch)
tsne <- Rtsne(t(Mergedcleanmatrix), dims = 2, perplexity = 50, verbose = TRUE, max_iter = 500)
plot(tsne$Y, main = " ", col = alpha(colors[batch], 0.3), type = 'p', pch = dottype, cex = 1.4, cex.axis = 2, cex.main = 2, font.axis = 2, xlab = '', ylab = '')
legend("topright", title = 'Batch', pch = 16, legend = unique(batch), col = colors, ncol = 3, cex = 0.7, text.font = 1, pt.cex = 1.5)
dev.off()

status <- as.factor(Alltrait$Status[match(colnames(Mergedcleanmatrix), table = Alltrait$Accession)])
design <- model.matrix(~0+status)
Normalized <- as.data.frame(removeBatchEffect(Mergedcleanmatrix, batch = batch, design = design))

svg(filename = "Colon liver data after batch correction.svg", height = 5, width = 6)
tsne<-Rtsne(t(Normalized), dims = 2, perplexity = 50, verbose = TRUE, max_iter = 500)
plot(tsne$Y, main = " ", col = alpha(colors[batch], 0.3), type = 'p', pch = dottype, cex = 1.4, cex.axis = 2, cex.main = 2, font.axis = 2, xlab = '', ylab = '')
legend('topright', title = 'Batch', pch = 16, legend = unique(batch), col = colors, ncol = 3, cex = 0.7, text.font = 1, pt.cex = 1.5)
dev.off()

write.csv(Alltrait, file = "ColonliverAlltrait.csv", row.names = F)
write.csv(Normalized, file = "ColonliverNormalized.csv")
write.csv(2^Mergedcleanmatrix, file = "Mergedcleanmatrix_original_scale.csv")