# geneSet should be a dataframe containing entrez IDs, log2FoldChange, and pvalues for a set of genes

formatGeneList <- function(geneSet) {
  df <- geneSet
  df <- df[,c("entrez","log2FoldChange")]
  df <- df[complete.cases(df), ]
  df <- df[order(-df$log2FoldChange),]
  geneList <- df$log2FoldChange
  names(geneList) <- df$entrez
  return(geneList)
}

formatGenes <- function(geneSet) {
  df <- geneSet
  df <- df[,c("entrez","log2FoldChange","pvalue")]
  df <- df[complete.cases(df), ]
  df <- df[order(-df$log2FoldChange),]
  dfSig <- subset(df, pvalue < 0.05)
  genes <- dfSig$log2FoldChange
  names(genes) <- dfSig$entrez
  return(names(genes))
}