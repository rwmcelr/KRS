## Library calls/ set script variables -------------------------------------
setwd("/Users/rmcelroy/Desktop/RNAseq/")
library(DESeq2)
library(pheatmap)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)

## Functions ---------------------------------------------------------------
# Identify significant genes (user selected paramaters) and output MA plot, heatmap, and filtered results
kshvSigGenes <- function(name, resNoShrink, resShrunk, filterBy = "padj", filterVal = 0.05) {
    res <- resShrunk
    resNS <- resNoShrink
  
  # Extract genes for specified value (padj is default), then order by descending log2foldchange
  if (filterBy == "pvalue") {
    resSig <- res[ which(res$pvalue < filterVal),]
    resNS <- resNS[which(res$pvalue < filterVal),]
  } else {
    resSig <- res[ which(res$padj < filterVal),]
    resNS <- resNS[which(res$padj < filterVal),]
  }
  resSig <- resSig[order(resSig$log2FoldChange, decreasing=TRUE),]
  resNS <- resNS[order(resNS$log2FoldChange, decreasing=TRUE),]
  
  # Add Entrez ID for downstream pathway analysis
  ens.str <- substr(rownames(resSig), 1, 15)
  resSig$entrez <- mapIds(org.Hs.eg.db,
                          keys=ens.str,
                          column="ENTREZID",
                          keytype="ENSEMBL",
                          multiVals="first")
  
  # Write csv file of results (significant gene list)
  write.csv(resSig, file=paste0(name,"_DEGlist.csv"))
  
  # Create gene set with symbols and non shrunken metrics for downstream pathway exploration analysis
  pathways<-resNoShrink[,c("symbol","log2FoldChange","pvalue","padj")]
  ens.str2 <- substr(rownames(pathways), 1, 15)
  pathways$entrez <- mapIds(org.Hs.eg.db,
                          keys=ens.str2,
                          column="ENTREZID",
                          keytype="ENSEMBL",
                          multiVals="first")
  write.csv(pathways, file="PathwayData.csv", row.names=FALSE)
  
  # Create ranked gene list (.rnk) for use with GSEA (ranks are directional log10 p value)
  gs <- resSig[,c("symbol","log2FoldChange","pvalue","entrez")]
  gs$fcsign <- sign(gs$log2FoldChange)
  gs$logP=-log10(gs$pvalue)
  gs$metric= gs$logP/gs$fcsign
  write.table(gs[,c("symbol","metric")],file="GeneSetSymbol.rnk",quote=F,sep="\t",row.names=F)
  write.table(gs[,c("symbol","entrez")],file="GeneSetEntrez.rnk",quote=F,sep="\t",row.names=F)
  
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
