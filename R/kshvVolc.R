## Library calls/ set script variables -------------------------------------
library(ggplot2)
library(ggrepel)
library(scales)
tf <- read.csv("tf.csv",row.names=1)

## Functions ---------------------------------------------------------------
# Identify significant genes (user selected paramaters) and output MA plot, heatmap, and filtered results
kshvVolc <- function(name, plots=T, filBy="padj", filVal=0.05, con1, con2) {
    resMA <- results(dds, contrast=c("condition",con1,con2))
    maPoints <- as.data.frame(resMA)
    maPoints <- merge(x = maPoints, y = tf, by = "symbol", all.x = TRUE, no.dups = TRUE)
    maPoints$identity[maPoints$pvalue.x < 0.05 & maPoints$identity == "NM23H2 Transcription Factor"] <- "Significant Transcription Factor"
    maPoints$negLog2pvalue <- log2(maPoints$pvalue.x) * -1
    lineLabels <- data.frame(x=c(-4.5,-4.4),y=c(1.2,4.522),label=c("P Value 0.5 Cutoff","P Value 0.05 Cutoff"),identity=c("P Val Cutoff","P Val Cutoff"))
    
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
      geom_hline(yintercept = c(4.322, 15.1), colour="skyblue", size=1, linetype="longdash") + 
      geom_text(data=lineLabels, aes(x=x, y=y, label=label)) +
      labs(x="Log2FoldChange", y="-Log2pvalue") + 
      scale_colour_manual(name="Key", values=c("NM23H2 Transcription Factor"="orange",
                                               "Significant Transcription Factor"="darkorchid1",
                                               "P Val Cutoff"="black"), na.value="grey80") + 
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
    print(plot2)
    dev.off()
}