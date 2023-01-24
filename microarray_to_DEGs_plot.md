# Differential expression analysis with limma
## Set up 
```{r}
setwd("~VertebrateDynamics/MAanalysis")
```
## Bioconductor nd CRAN libraries used
```{r}
library(Biobase)
library(GEOquery)
library(limma)
library(openxlsx)
```

## Load the data from GEO platform
Write the Accession number of the target dataset. Then we check taht only one dataset is retrieved when using this code.
```{r}
gset <- getGEO("YourACCESSIONNumber", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("YourPlatformID", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
```
## Data pre-processing
Edit column names to match toptable and group the samples (assign control, 0 or treated, 1). The non chosen samples are indicated with X.
```{r}
fvarLabels(gset) <- make.names(fvarLabels(gset))
# assign control, treatment or non-selected sample code.    
gsms <- "XXXXXX000XXX111XXXXXXXXX"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
# Remove unnecessary samples labelled as X.
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]
```
### Log2 transform
```{r}
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }
```

## DEGs analysis
```{r}
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=50000)

write.xlsx(tT, "TT.xlsx", row.names=F, sep="\t")
```


