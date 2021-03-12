## Library calls/ set script variables -------------------------------------
setwd("/Users/rmcelroy/Desktop/RNAseq/")
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(ggrepel)
library(scales)
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
    maPoints$identity[maPoints$pvalue.x < 0.05 & maPoints$identity == "NM23H2 Transcription Factor"] <- "Significant Transcription Factor"
    maPoints$negLog2pvalue <- log2(maPoints$pvalue.x) * -1
    png(file=paste0(name,"_VolcanoPlot.png"),width=1500,height=1500,res=150)
    plot2 <- ggplot(maPoints, aes(log2FoldChange.x, negLog2pvalue, colour=identity, label=symbol)) + 
      ggtitle(paste0(name,"\nVolcano Plot")) +
      geom_point(size=1) +
      geom_label_repel(aes(label=ifelse(identity=="Significant Transcription Factor",as.character(symbol),'')),
                       box.padding = 0.5,
                       point.padding = 0.5,
                       max.overlaps = Inf,
                       segment.color = 'grey50') +
      scale_y_continuous(limits=c(0, 15), oob=squish) + 
      scale_x_continuous(limits=c(-5, 5), breaks=seq(-5,5,1)) + 
      geom_hline(yintercept = 1, colour="yellow", size=1, linetype="longdash") + 
      geom_text(label="P Value 0.5 Cutoff", x=-4.5, y=1.2, color="black") +
      geom_hline(yintercept = c(4.322, 15.1), colour="skyblue", size=1, linetype="longdash") + 
      geom_text(label="P Value 0.05 Cutoff", x=-4.4, y=4.522, color="black") +
      labs(x="Log2FoldChange", y="-Log2pvalue") + 
      scale_colour_manual(name="Key", values=c("NM23H2 Transcription Factor"="orange",
                                               "Significant Transcription Factor"="darkorchid1"), na.value="grey80") + 
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
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
sigGenes("Wild Type vs Vector", T, "pvalue", 0.05, "wt", "vec")
sigGenes("Wild Type vs E14A Mutant", T, "pvalue", 0.05, "wt", "mut")
sigGenes("E14A Mutant vs Vector", T, "pvalue", 0.05, "mut", "vec")