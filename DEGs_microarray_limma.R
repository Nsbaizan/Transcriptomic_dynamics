# Differential expression analysis with limma
###Libraries---- 
library(Biobase)
library(limma)
library(tidyverse)

library(kableExtra)
library(tinytex)
#----

###Data importation----
Mus_Males <- read.csv2("./Mus_mus_mascles.csv", row.names = 1)
#as matrix
exprsM <- as.matrix(Mus_Males, header=TRUE, sep="\t", row.names=T, as.is=TRUE)
#transform to ExpressionSet
gsetM <- ExpressionSet(assayData=exprsM)
#PhenoData
pDataM <- read.delim("./PDataM.txt", row.names=1, header=TRUE, sep="\t")
metadata <- data.frame(labelDescription= c("Age"))
phenoData <- new("AnnotatedDataFrame", data=pDataM, varMetadata=metadata)
#EnterPhenoData to the ExpressionSet
gsetM <- ExpressionSet(assayData=exprsM, phenoData=phenoData)
#----

###Check if expression values are Log2 transfored, if not, transform them----
fvarLabels(gsetM) <- make.names(fvarLabels(gsetM))
gsms <- "000111"            
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
sel <- which(sml != "X")
sml <- sml[sel]
gsetM <- gsetM[ ,sel]
ex <- exprs(gsetM)

qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
qx

LogC <- (qx[5] > 100) ||
  +   (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  +   (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
LogC

if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gsetM) <- log2(ex) }
#----

##DEGs analysis----
fl <- as.factor(sml)
gsetM$description <- fl
design <- model.matrix(~ description + 0, gsetM)
colnames(design) <- levels(fl)
fit <- lmFit(gsetM, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
TT_Mus_mus_XY <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
TT_Mus_mus_XY <- tibble::rownames_to_column(TT_Mus_mus_XY, var = "ID")
#----

###Results annotation----
Mus_Annot <- read.csv2("./Mus_Anotacio.csv")
Mus_Male_Def_TT <- merge(Mus_Annot, TT_Mus_mus_XY, by = "ID")
#----
