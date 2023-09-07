## Library calls/ set script variables -------------------------------------
setwd("/Users/rmcelroy/Desktop/PEA/")
library("clusterProfiler")
library("pathview")
library("ggplot2")
library("enrichplot")
library("ReactomePA")

# ## WikiPathways Analysis --------------------------------------------
# # Gene Set Enrichment Analysis
# wiPath2 <- gseWP(wtMut_geneList, organism = "Homo sapiens")
# write.csv(wiPath2, file="WtVsMut_WikiPathways_GSEA.csv")
# dotplot(wiPath2) + ggtitle("Wt Vs Mut WikiPathways GSEA Dotplot")
# wiPath2 <- pairwise_termsim(wiPath2)
# p4 <- emapplot(wiPath2) + ggtitle("Wt Vs Mut WikiPathways GSEA Map")
# p4
# 
# ## Reactome Pathways Analysis --------------------------------------------
# reactome2 <- gsePathway(wtMut_geneList,
#                         pvalueCutoff = 0.2,
#                         pAdjustMethod = "BH",
#                         verbose = FALSE)
# write.csv(reactome2, file="Wt Vs Mut ReactomePathways GSEA.csv")
# dotplot(reactome2) + ggtitle("Wt Vs Mut Reactome GSEA Dotplot")
# reactome2 <- pairwise_termsim(reactome2)
# p6 <- emapplot(reactome2) + ggtitle("Wt Vs Mut ReactomePathways GSEA Map")
# p6