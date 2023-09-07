## Library calls/ set script variables -------------------------------------
source("/Users/rmcelroy/Desktop/KRS/R/format.R")
source("/Users/rmcelroy/Desktop/KRS/R/PEA_ORA.R")
source("/Users/rmcelroy/Desktop/KRS/R/PEA_GSEA.R")

## Main Body --------------------------------------------------------------
PEA <- function(baseName) {
  current <- getwd()
  setwd(paste0(current,"/",baseName,"/"))
  
  dat <- read.csv("PathwayData.csv")
  geneList <- formatGeneList(dat)
  genes <- formatGenes(dat)
  try(PEA_ORA(baseName, geneList, genes))
  try(PEA_ORA(baseName, geneList, genes))
  setwd(current)
}