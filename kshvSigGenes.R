## Library calls/ set script variables -------------------------------------
setwd("/Users/rmcelroy/Desktop/RNAseq/")
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)

## Functions ---------------------------------------------------------------
# Identify significant genes (user selected paramaters) and output MA plot, heatmap, and filtered results
sigGenes <- function(name, resNoShrink, resShrunk, filterBy = "padj", filterVal = 0.05) {
    res <- resShrunk
    resNS <- resNoShrink
  
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
  
  # Write csv file of results (significant gene list)
  write.csv(resSig, file=paste0(name,"_DEGlist.csv"))
  
  # Create gene set with symbols and non shrunken metrics for downstream pathway exploration analysis
  pathways<-resSig[,c("symbol","log2FoldChange","pvalue","padj")]
  pathways$log2FoldChange <- resNS$log2FoldChange
  write.csv(pathways, file="GeneSet.csv", row.names=FALSE)
  
  # Create ranked gene list (.rnk) for use with GSEA (ranks are directional log10 p value)
  gs <- resSig[,c("symbol","log2FoldChange","pvalue")]
  gs$fcsign <- sign(gs$log2FoldChange)
  gs$logP=-log10(gs$pvalue)
  gs$metric= gs$logP/gs$fcsign
  write.table(gs[,c("symbol","metric")],file="GeneSet.rnk",quote=F,sep="\t",row.names=F)
  setwd(current)
  
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