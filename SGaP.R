library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)

SGaP <- function(name, dds, filterBy="padj", filterVal=0.05, con1, con2, GS = FALSE) {
  # Create and annotate shrunken dds results object
  resShrunk <- lfcShrink(dds=dds, contrast=c("condition",con1,con2), type="ashr")
  ens.strS <- substr(rownames(resShrunk), 1, 15)
  resShrunk$symbol <- mapIds(EnsDb.Hsapiens.v86,keys=ens.strS,column="SYMBOL",keytype="GENEID",multiVals="first")
  # Create and annotate non-shrunken dds results object
  res <- results(dds, contrast=c("condition",con1,con2))
  ens.str <- substr(rownames(res), 1, 15)
  res$symbol <- mapIds(EnsDb.Hsapiens.v86,keys=ens.str,column="SYMBOL",keytype="GENEID",multiVals="first")
  
  # Create gene set object if provided in the original working directory
  if (file.exists("GeneSet.csv")) {gs <- read.csv("GeneSet.csv")}
  
  # Create directory for results and set to working directory
  current <- getwd()
  if (!dir.exists(name)) {  dir.create(name)  }
  setwd(paste0(current,"/",name,"/"))
  
  # Find gene level statistics for a target list of genes in RNA seq data
  if (GS = TRUE) {
    getGSdata(res, gs)
  }
  # Find significant DEGs, and generate: lfc shrunk list(csv), a non shrunken list(csv), heatmap, and a .rnk file
  kshvSigGenes(name, res, resShrunk, filterBy, filterVal)
  # Create a MA plot for the DEGs, with target genes highlighted if present
  kshvMA(name, res)
  # Create a Volcano plot for the DEGs, with target genes highlighted if present
  kshvVolc(name, res)
  
  # Set the working directory back to the original before exiting the function
  setwd(current)
}
  