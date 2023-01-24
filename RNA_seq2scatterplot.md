## Gene-level differential expression analysis using DESeq2

## Setup
```{r}
setwd("~VertebrateDynamics/DEanalysis")
```
## Bioconductor and CRAN libraries used
```{r}
library(tidyverse)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(DEGreport)
```
## Load in data
The data is stored in the working directory inside a folder named 'data' and the metadata including sample information is found in a separate folder within wd named 'meta'.
```{r}
data1 <- read.table("data/Mov10_full_counts.txt", header=T, row.names=1) 
data <- read.csv("data/Chicken_matrix_FC.csv", header=T,sep = ";",row.names=1)

meta1 <- read.table("meta/Mov10_full_meta.txt", header=T, row.names=1)
meta <- read.csv("meta/Sample_info_Chicken.csv", header=T,sep = ";",row.names=1)
```
Check classes of the data
```{r}
class(meta)
class(data)
```
## RNA seq count distribution 
These images illustrate some common features of RNA-seq count data, including a low number of counts associated with a large proportion of genes, 
and a long right tail due to the lack of any upper limit for expression.
```{r}
ggplot(data) +
  geom_histogram(aes(x = F45S1), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")
```
Check which model type should be used: Poisson or negative binomial distribution.
If the variance across replicates tends to be greater than the mean (red line), do not use Poisson. We use negative bionmial distribution.

```{r}
mean_counts <- apply(data[, 3:5], 1, mean)
variance_counts <- apply(data[, 3:5], 1, var)
df <- data.frame(mean_counts, variance_counts)

ggplot(df) +
  geom_point(aes(x=mean_counts, y=variance_counts)) + 
  geom_line(aes(x=mean_counts, y=mean_counts, color="red")) +
  scale_y_log10() +
  scale_x_log10()
``` 
## Count normalization
Check that sample names match in both files and create the DESeq2 object. Normalize counts.

```{r}
all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))


dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)
View(counts(dds))

dds <- estimateSizeFactors(dds)
sizeFactors(dds)

normalized_counts <- counts(dds, normalized=TRUE)
``` 

## Create DESeq object
```{r}
dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)
## Run analysis
dds <- DESeq(dds)
## Check the size factors
sizeFactors(dds)
## Total number of raw counts per sample
colSums(counts(dds))
## Total number of normalized counts per sample
colSums(counts(dds, normalized=T))
## Plot dispersion estimates
plotDispEsts(dds)
``` 

## Define contrasts, extract results table
```{r}
contrast_F <- c("sampletype", "F_late", "F_early")
contrast_M <- c("sampletype", "M_late", "M_early")

resF <- results(dds, contrast = contrast_F)
summary(resF)

resM <- results(dds, contrast = contrast_M)
summary(resM)
``` 

# Make a basic volcano plot
```{r}
par(mfrow=c(1,1))
with(resF, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", ylim=c(0,90)))
# Add colored points: red if log2FC>1 and padj<0.05)
with(subset(resF, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#For males
par(mfrow=c(1,1))
with(resM, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", ylim=c(0,90)))
with(subset(resM, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))


FDEG<-as.data.frame(resF)
MDEG<-as.data.frame(resM)
DEGs<-cbind(FDEG,MDEG)

write.csv(as.data.frame(DEGs), 
          file="results/ChickenDEGs_2022.csv")
``` 

## Make the scatterplot
```{r}
Chicken <- DEGs[ -c(1,3:5,7,9:11) ]

colnames(Chicken)<- c("logFC.x","adj.P.Val.x","logFC.y", "adj.P.Val.y" )

Chicken$logFC.x<-as.numeric(Chicken$logFC.x)
Chicken$logFC.y<-as.numeric(Chicken$logFC.y)
library(dplyr)

Chicken$Diff <- Chicken$logFC.x - Chicken$logFC.y



Chicken <- mutate(Chicken, pattern =
                  ifelse(Chicken$Diff < -0.585, "M",
                         ifelse(Chicken$Diff > 0.585, "F",
                                ifelse( -0.585 <= Chicken$Diff & Chicken$Diff <= 0.585, "Both", NA))))


Chicken <- mutate(Chicken, ExpressionChange=
                  ifelse(Chicken$pattern == "Both", "Both",
                         ifelse(Chicken$pattern == "M" & abs(Chicken$logFC.y) > abs(Chicken$logFC.x), "Up in Male",
                                ifelse(Chicken$pattern == "M" & abs(Chicken$logFC.x) > abs(Chicken$logFC.y), "Down in Female",
                                       ifelse(Chicken$pattern == "F" & abs(Chicken$logFC.x) > abs(Chicken$logFC.y), "Up in Female",
                                              ifelse(Chicken$pattern == "F" & abs(Chicken$logFC.y) > abs(Chicken$logFC.x), "Down in Male", NA))))))





Chicken1<-filter(Chicken, logFC.x>=0.585|logFC.y>=0.585)
Chicken2<-filter(Chicken, logFC.x<=-0.585|logFC.y<=-0.585)
ChickenFC<-rbind(Chicken1,Chicken2)
Chicken<-filter(ChickenFC, adj.P.Val.x<=0.05|adj.P.Val.y<=0.05)



library(plotly)
x <- list(
  title = "Fold Change in Females", tick0 = 0, dtick = 2)
y <- list(
  title = "Fold Cahnge in Males")


p <- plot_ly(Chicken, x = ~logFC.x, y= ~logFC.y, 
             color = ~pattern, 
             colors = (c("darkgrey", "red", "blue")),
             type = "scatter", mode="markers")%>%
  layout(xaxis = x, yaxis = y)

p

write.csv(as.data.frame(Chicken), 
          file="results/ChickenDEGsS2vsS1_2022.csv")
```

## References
The code used for the analysis of DEGs was adapted from: 
https://github.com/hbctraining/DGE_workshop/blob/master/schedule/1.5-day.md
