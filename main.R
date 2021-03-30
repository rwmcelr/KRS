## Library calls/ set script variables -------------------------------------
setwd("/Users/rmcelroy/Desktop/RNAseq/")
library(DESeq2)
library(gplots)
source("/Users/rmcelroy/Desktop/RNAseq/R/getGSdata.R")
source("/Users/rmcelroy/Desktop/RNAseq/R/kshvMA.R")
source("/Users/rmcelroy/Desktop/RNAseq/R/kshvSigGenes.R")
source("/Users/rmcelroy/Desktop/RNAseq/R/kshvVolc.R")
source("/Users/rmcelroy/Desktop/RNAseq/R/SGaP.R")

## Import & pre-process ----------------------------------------------------
# Import data from featureCounts, clean and convert to matrix
countdata <- read.table("counts.txt", header=TRUE, row.names=1)
countdata <- countdata[ ,6:ncol(countdata)]
colnames(countdata) <- gsub("\\.[sb]am$", "", colnames(countdata))
countdata <- as.matrix(countdata)

# Assign condition
condition <- factor(c(rep("mut", 3), rep("vec", 3), rep("wt", 3)))

## Analysis with DESeq2 ----------------------------------------------------
# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
coldata <- data.frame(row.names=colnames(countdata), condition)
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds$condition <- relevel(dds$condition, "vec")

# Filter out genes with hits below 10
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run the DESeq pipeline
dds <- DESeq(dds)

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)

## Generate general plots -------------------------------------------------------
# Principal components analysis (need to use normalized data)
png(filename="PCA.png",width=8,height=8,units="in",res=500)
DESeq2::plotPCA(rld, intgroup="condition")
dev.off()

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
colorGradient <- colorRampPalette(c("black","yellow","white"))(n = 299)
png(filename="sampleDists.png",width=8,height=8,units="in",res=500)
heatmap.2(sampleDists, col = colorGradient, trace = "none", key = FALSE)
dev.off()

## Generate significant genes, heatmap, MA plot, and Volcano plot for specified condition ----------------------------
SGaP("Wt Vs Vector p < 0.05", dds, "pvalue", 0.05, "wt", "vec", GS=TRUE)
SGaP("Wt Vs Vector q < 0.05", dds, "padj", 0.05, "wt", "vec")

SGaP("Wt Vs E14A Mutant p < 0.05", dds, "pvalue", 0.05, "wt", "mut", GS=TRUE)
SGaP("Wt Vs E14A Mutant q < 0.05", dds, "padj", 0.05, "wt", "mut")
 
SGaP("E14A Mutant Vs Vector p < 0.05", dds, "pvalue", 0.05, "mut", "vec", GS=TRUE)
SGaP("E14A Mutant Vs Vector q < 0.05", dds, "padj", 0.05, "mut", "vec")