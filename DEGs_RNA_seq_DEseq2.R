###Working  directory and libraries
setwd("SetYourPATH/DEanalysis/")

library(DESeq2)
library(dplyr)
###----

###Import the data----
#Count matrix
counts <- read.csv2("./data/Chicken_matrix_FC.csv", row.names = 1)

#Experimental design
experiment_des <- read.csv2("./meta/Sample_info_Chicken.csv", sep = ",", row.names = 1)

#Design formula
des_form <- as.formula(~ sampletype)

#Annotation
annot <- read.csv2("./results/annotations/genename_chicken.csv")
###----


###Check all is OK----
all(colnames(counts) == rownames(experiment_des))
###----

###Apply Deseq2----
#Create the DeseqObject
Chicken <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = experiment_des,
                                  design = des_form)

#Filtering out nonimportant genes----
Chicken <- Chicken[rowSums(counts(Chicken)) > 10,]

#Using the DESEQ function
Chicken <- DESeq(Chicken)
vsd <- vst(Chicken, blind = FALSE)
vsd <- assay(vsd)
vsd <- as.data.frame(vsd)

#Female
female <- results(Chicken, contrast = c("sampletype",'F_late','F_early'))
female_res <- as.data.frame(female)
female_res <- na.exclude(female_res)

#Male
male <- results(Chicken, contrast = c("sampletype", "M_early","M_late"))
male_res <- as.data.frame(male)
male_res <- na.exclude(male_res)

#Merging and shaping
res <- merge(female_res, male_res, by=0)
rownames(res) <- res[,1]
res <- res[,-1]
#----


###Gene expression changes, filtering and annotation----
#Expression Change 
res$Diff <- res$log2FoldChange.x - res$log2FoldChange.y
res <- mutate(res, pattern =
                ifelse(res$Diff < -0.5, "M",
                       ifelse(res$Diff > 0.5, "F",
                              ifelse( -0.5 <= res$Diff & res$Diff <= 0.5, "Both", NA))))

res <- mutate(res, ExpressionChange=
                ifelse(res$pattern == "Both", "Both",
                       ifelse(res$pattern == "M" & abs(res$log2FoldChange.y) > abs(res$log2FoldChange.x), "Up in Male",
                              ifelse(res$pattern == "M" & abs(res$log2FoldChange.x) > abs(res$log2FoldChange.y), "Down in Female",
                                     ifelse(res$pattern == "F" & abs(res$log2FoldChange.x) > abs(res$log2FoldChange.y), "Up in Female",
                                            ifelse(res$pattern == "F" & abs(res$log2FoldChange.y) > abs(res$log2FoldChange.x), "Down in Male", NA))))))
#Gene Filtering
res <- res[(res$padj.x <= 0.05 | res$padj.y <= 0.05),]

#logFC
res <- res[(abs(res$log2FoldChange.x) >= 0.585 | abs(res$log2FoldChange.y) >= 0.585),]

#Gene annotation
res_chick <- merge(x = res, y = annot,
                   by.x = 0, by.y = 1)
###----

## References
The code used for the analysis of DEGs was adapted from: 
https://github.com/hbctraining/DGE_workshop/blob/master/schedule/1.5-day.md
