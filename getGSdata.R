## Main Function -----------------------------------------------------------
getGSdata <- function(res, GeneSet) {
  temp <- res
  gs <- GeneSet
  resGS <- as.data.frame(temp)
  resGS$important <- row.names(temp)
  gsFinal <- merge(x = resGS, y = gs, by = "symbol")
  row.names(gsFinal) <- gsFinal$important
  gsFinal$important <- gsFinal$symbol
  gsFinal$symbol <- NULL
  colnames(gsFinal)[which(names(gsFinal) == "important")] <- "symbol"
  gsFinal$identity <- "Target Gene"
  write.csv(gsFinal, file="GeneSetData.csv")
}