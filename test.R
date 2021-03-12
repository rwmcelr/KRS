## Library calls/ set script variables -------------------------------------
setwd("/Users/rmcelroy/Desktop/RNAseq/")
library(DESeq2)
library(RColorBrewer)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)
library(pheatmap)
library(ggplot2)
library(scales)
library(viridis)
tf <- read.csv("tf.csv",row.names=1)

## Functions ---------------------------------------------------------------
# Identify significant genes (user selected paramaters) and output MA plot, heatmap, and filtered results
sigGenes <- function(name, plots=T, filBy="padj", filVal=0.05, con1, con2) {
  # plots = T (default), type = normal (default) or shrink, filBy = padj (default) or pvalue, filVal = 0.05 (default)
  # Get differential expression results, both shrunken and not shrunken
  
  if (plots) {
    resMA <- results(dds, contrast=c("condition",con1,con2))
    ens.str <- substr(rownames(resMA), 1, 15)
    resMA$symbol <- mapIds(EnsDb.Hsapiens.v86,
                            keys=ens.str,
                            column="SYMBOL",
                            keytype="GENEID",
                            multiVals="first")
    resMA$entrez <- mapIds(org.Hs.eg.db,
                            keys=ens.str,
                            column="ENTREZID",
                            keytype="ENSEMBL",
                            multiVals="first")
    maPoints <- as.data.frame(resMA)
    maPoints <- merge(x = maPoints, y = tf, by = "symbol", all.x = TRUE, no.dups = TRUE)
    select <- maPoints$padj < 0.2
    maPoints[select,"identity"] <- "Significant"
    write.csv(maPoints, file="maPoints.csv")
    png(file="MA_plot.png",width=768,height=768)
    plot2 <- ggplot(maPoints, aes(baseMean.x, log2FoldChange.x, colour=identity)) + geom_point(size=1) + scale_y_continuous(limits=c(-1, 1), oob=squish) + scale_x_log10() + geom_hline(yintercept = 0.5, colour="blue", size=1, linetype="dotted") + geom_hline(yintercept = -0.5, colour="blue", size=1, linetype="dotted") + labs(x="Mean of Normalized Counts", y="Log2 Fold Change") + scale_colour_manual(name="TF", values=("NM23H2 Transcription Factor"="yellow"), na.value="grey50") + theme_bw()
    print(plot2)
    dev.off()
    }
}

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

## Generate significant genes, MA plot, and sig gene heat map for specified conditions ----------------------------
sigGenes("test", T, "pvalue", 0.05, "wt", "vec")
