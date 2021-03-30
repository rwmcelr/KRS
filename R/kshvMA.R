## Library calls/ set script variables -------------------------------------
library(ggplot2)
library(ggrepel)
library(scales)

## Functions ---------------------------------------------------------------
# Identify significant genes (user selected paramaters) and output MA plot, heatmap, and filtered results
kshvMA <- function(name, res) {
  maPoints <- as.data.frame(res)
  maPoints$identity[maPoints$pvalue < 0.05] <- "Significant (p < 0.05)"
  
  # Finds genes which are both significant and part of supplied target gene list, and changes their identity appropriately
  if (file.exists("TargetGeneSetData.csv")) {
    tf <- read.csv("TargetGeneSetData.csv",row.names=1) 
    maPoints <- merge(x = maPoints, y = tf, by = "symbol", all.x = TRUE)
    colnames(maPoints)[which(names(maPoints) == "identity.x")] <- "identity"
    maPoints$identity[maPoints$pvalue.x < 0.05 & maPoints$identity.y == "Target Gene"] <- "Significant Target Gene"
    maPoints$identity[is.na(maPoints$identity) & maPoints$identity.y == "Target Gene"] <- "Target Gene"
    p <- ggplot(maPoints, aes(baseMean.x, log2FoldChange.x, colour=identity, label=symbol))
    key <- c("Target Gene"="orange","Significant (p < 0.05)"="skyblue","Significant Target Gene"="darkorchid1")
  } else { 
    p <- ggplot(maPoints, aes(baseMean, log2FoldChange, colour=identity, label=symbol))
    key <- c("Significant (p < 0.05)"="skyblue")
  }
  
  # Adding elements to the plots
  png(file=paste0(name,"_MAplot.png"),width=1500,height=1500,res=150)
  p <- p + ggtitle(paste0(name,"\nMA Plot")) +
    geom_point(size=1) +
    geom_label_repel(aes(label=ifelse(is.na(identity) == FALSE & identity!="Significant (p < 0.05)",as.character(symbol),'')),
                     box.padding = 0.5,
                     point.padding = 0.5,
                     max.overlaps = Inf,
                     segment.color = 'grey50') +
    scale_y_continuous(limits=c(-3, 3), oob=squish) + 
    scale_x_log10() + 
    labs(x="Mean of Normalized Counts", y="Log2 Fold Change") + 
    scale_colour_manual(name="Key", values=key, na.value="grey80") + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  print(p)
  dev.off()
}