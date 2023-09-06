## Library calls/ set script variables -------------------------------------
library(ggplot2)
library(ggrepel)
library(scales)

## Functions ---------------------------------------------------------------
# Identify significant genes (user selected paramaters) and output MA plot, heatmap, and filtered results
VOLCplot <- function(name, res) {
  volcPoints <- as.data.frame(res)
  volcPoints$identity[volcPoints$pvalue < 0.05] <- "Significant (p < 0.05)"
  
  # Finds genes which are both significant and part of supplied target gene list, and changes their identity appropriately
  if (file.exists("TargetGeneSetData.csv")) {
    tgList <- read.csv("TargetGeneSetData.csv",row.names=1) 
    volcPoints <- merge(x = volcPoints, y = tgList, by = "symbol", all.x = TRUE)
    colnames(volcPoints)[which(names(volcPoints) == "identity.x")] <- "identity"
    volcPoints$identity[volcPoints$pvalue.x < 0.05 & volcPoints$identity.y == "Target Gene"] <- "Significant Target Gene"
    volcPoints$identity[is.na(volcPoints$identity) & volcPoints$identity.y == "Target Gene"] <- "Target Gene"
    volcPoints$negLog2pvalue <- log2(volcPoints$pvalue.x) * -1
    p <- ggplot(volcPoints, aes(log2FoldChange.x, negLog2pvalue, colour=identity, label=symbol))
    key <- c("Target Gene"="orange","Significant (p < 0.05)"="skyblue","Significant Target Gene"="darkorchid1")
  } else { 
    volcPoints$negLog2pvalue <- log2(volcPoints$pvalue) * -1
    p <- ggplot(volcPoints, aes(log2FoldChange, negLog2pvalue, colour=identity, label=symbol))
    key <- c("Significant (p < 0.05)"="skyblue")
  }
  
  # Creates labels for p value cuttof lines, if used
  #lineLabels <- data.frame(x=c(-4.5,-4.4),y=c(1.2,4.522),label=c("P Value 0.5 Cutoff","P Value 0.05 Cutoff"),identity=c("P Val Cutoff","P Val Cutoff"))
  
  # Adding elements to the plots
  png(file=paste0(name,"_VolcanoPlot.png"),width=1500,height=1500,res=150)
  p <- p + ggtitle(paste0(name,"\nVolcano Plot")) +
    geom_point(size=1) +
    geom_label_repel(aes(label=ifelse(is.na(identity) == FALSE & identity!="Significant (p < 0.05)",as.character(symbol),'')),
                     box.padding = 0.5,
                     point.padding = 0.5,
                     max.overlaps = Inf,
                     segment.color = 'grey50') +
    scale_y_continuous(limits=c(0, 15), oob=squish) + 
    scale_x_continuous(limits=c(-5, 5), breaks=seq(-5,5,1)) + 
    #geom_hline(yintercept = 1, colour="yellow", size=1, linetype="longdash") + 
    #geom_hline(yintercept = c(4.322, 15.1), colour="skyblue", size=1, linetype="longdash") + 
    #geom_text(data=lineLabels, aes(x=x, y=y, label=label)) +
    labs(x="Log2 Fold Change", y="-Log2 p-value") + 
    scale_colour_manual(name="Key", values=key, na.value="grey80") + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  print(p)
  dev.off()
}
