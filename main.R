## Library calls/ set script variables -------------------------------------
setwd("/Users/rmcelroy/Desktop/RNAseq/")
library(DESeq2)
library(RColorBrewer)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)
library(pheatmap)
library(gplots)

## Functions ---------------------------------------------------------------
# Identify significant genes (user selected paramaters) and output MA plot, heatmap, and filtered results
sigGenes <- function(name, plots=T, filBy="padj", filVal=0.05, con1, con2) {
  # plots = T (default), type = normal (default) or shrink, filBy = padj (default) or pvalue, filVal = 0.05 (default)
  # Get differential expression results, both shrunken and not shrunken
    res <- lfcShrink(dds=dds, contrast=c("condition",con1,con2), type="ashr")
    resNS <- results(dds, contrast=c("condition",con1,con2))
  
  # Extract genes for specified value (padj is default), then order by descending log2foldchange
  if (filBy == "pvalue") {
    resSig <- res[ which(res$pvalue < filVal),]
    resNS <- resNS[which(res$pvalue < filVal),]
  } else {
    resSig <- res[ which(res$padj < filVal),]
    resNS <- resNS[which(res$padj < filVal),]
  }
  resSig <- resSig[order(resSig$log2FoldChange, decreasing=TRUE),]
  resNS <- resNS[order(resNS$log2FoldChange, decreasing=TRUE),]

  # Annotate results gene list with gene symbol and Entrez ID, for clarity
  ens.str <- substr(rownames(resSig), 1, 15)
  resSig$symbol <- mapIds(EnsDb.Hsapiens.v86,
                       keys=ens.str,
                       column="SYMBOL",
                       keytype="GENEID",
                       multiVals="first")
  resSig$entrez <- mapIds(org.Hs.eg.db,
                       keys=ens.str,
                       column="ENTREZID",
                       keytype="ENSEMBL",
                       multiVals="first")
  
  # Create directory for results and set to working directory
  current <- getwd()
  if (!dir.exists(name)) {  dir.create(name)  }
  setwd(paste0(current,"/",name,"/"))
  
  # Write csv file of results (significant gene list)
  write.csv(resSig, file=paste0(name,"_DEGlist.csv"))
  
  # Create gene set with symbols and non shrunken metrics for downstream pathway exploration analysis
  pathways<-resSig[,c("symbol","log2FoldChange","pvalue","padj")]
  pathways$log2FoldChange <- resNS$log2FoldChange
  write.csv(pathways, file="GeneSet.csv", row.names=FALSE)
  
  if (plots) {
    # Create and save MA plot
    df <- res
    png(file=paste0(name,"_MAplot.png"),width=768,height=768)
    plotMA(df, ylim=c(-2,2), colSig = "red", main=name)
    dev.off()
    
    # Heatmap generation (needs work)
    # Generating a matrix to be used heatmap creation, labeling with gene symbol as opposed to ensembl id
    de <- rownames(resSig)
    de_mat <- assay(rld)[de,]
    de_mat <- de_mat - rowMeans(de_mat)
    gns <- select(EnsDb.Hsapiens.v86, row.names(de_mat), "SYMBOL", "GENEID")
    row.names(de_mat)[match(gns[,1], row.names(de_mat))] <- gns[,2]
    # Create and fix annotation list, then generate heatmap
    anno <- as.data.frame(colData(rld)[, c("condition")])
    rownames(anno) <- colnames(de_mat)
    names(anno)[1] <- "Condition"
    tallDim <- dim(resSig)[1]
    png(file=paste0(name,"_HeatMap.png"), width=480, height=tallDim*12+40)
    pheatmap(de_mat, annotation_col = anno, cluster_rows=FALSE, cluster_cols=FALSE)
    dev.off()
  }
  setwd(current)
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

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)

## Generate general plots -------------------------------------------------------
# Principal components analysis (need to use normalized)
png(filename="PCA.png",width=8,height=8,units="in",res=500)
DESeq2::plotPCA(rld, intgroup="condition")
dev.off()

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
png(filename="sampleDists.png",width=8,height=8,units="in",res=500)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()

## Generate significant genes, MA plot, and sig gene heat map for specified conditions ----------------------------
sigGenes("WtVsVecP", T, "pvalue", 0.05, "wt", "vec")
sigGenes("WtVsVecFDR", T, "padj", 0.05, "wt", "vec")

sigGenes("WtVsMutP", T, "pvalue", 0.05, "wt", "mut")
sigGenes("WtVsMutFDR", T, "padj", 0.05, "wt", "mut")

sigGenes("VecVsMutP", T, "pvalue", 0.05, "vec", "mut")
sigGenes("VecVsMutFDR", T, "padj", 0.05, "vec", "mut")