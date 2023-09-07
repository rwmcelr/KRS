## Library calls/ set script variables -------------------------------------
library("clusterProfiler")
library("ggplot2")
library("enrichplot")
library("ReactomePA")

PEA_ORA <- function(name, geneList, genes) {
  current <- getwd()
  if (!dir.exists(paste0(name,"_ORA_Results"))) {  dir.create(paste0(name,"_ORA_Results"))  }
  setwd(paste0(current,"/",name,"_ORA_Results/"))
  new <- getwd()
  
  ## GO (Gene Ontology) Analysis --------------------------------------------
  # BP = Biological Process, MF = Molecular Function
  if (!dir.exists("GO")) {  dir.create("GO")  }
  setwd(paste0(new,"/GO/"))
  goBP <- enrichGO(gene = genes,
                          universe = names(geneList),
                          OrgDb = "org.Hs.eg.db",
                          keyType = 'ENTREZID',
                          readable = T,
                          ont = "BP",
                          pvalueCutoff = 0.2)
  write.csv(goBP, file=paste0(name,"_GO(BP)_ORA.csv"))
  
  if(dim(goBP) > 1) {
    pdf(paste0(name,"_GO(BP)_Dotplot.pdf"), width = 8, height = 11)
    BPdot <- dotplot(goBP) + ggtitle(paste0(name,"\nGO(BP) ORA Dotplot"))
    print(BPdot)
    dev.off()
    
    goBP <- pairwise_termsim(goBP)
    pdf(paste0(name,"_GO(BP)_MAP.pdf"), width = 8, height = 11)
    BPmap <- emapplot(goBP, cex_label_category = 0.5, cex_category = 0.5, min_edge = 0.75) + ggtitle(paste0(name,"\nGO(BP) ORA Map"))
    print(BPmap)
    dev.off()
    
    pdf(paste0(name,"_GO(BP)_GOplot.pdf"), width = 8, height = 11)
    BPgo <- goplot(goBP)
    print(BPgo)
    dev.off()
  }
  
  goMF <- enrichGO(gene = genes,
                          universe = names(geneList),
                          OrgDb = "org.Hs.eg.db",
                          keyType = 'ENTREZID',
                          readable = T,
                          ont = "MF",
                          pvalueCutoff = 0.2)
  write.csv(goMF, file=paste0(name,"_GO(MF)_ORA.csv"))
  
  if(dim(goMF) > 1) {
    pdf(paste0(name,"_GO(MF)_Dotplot.pdf"), width = 8, height = 11)
    MFdot <- dotplot(goMF) + ggtitle(paste0(name,"\nGO(MF) ORA Dotplot"))
    print(MFdot)
    dev.off()
    
    goMF <- pairwise_termsim(goMF)
    pdf(paste0(name,"_GO(MF)_MAP.pdf"), width = 8, height = 11)
    MFmap <- emapplot(goMF, cex_label_category = 0.5, cex_category = 0.5, min_edge = 0.75) + ggtitle(paste0(name,"\nGO(MF) ORA Map"))
    print(MFmap)
    dev.off()
    
    pdf(paste0(name,"_GO(MF)_GOplot.pdf"), width = 8, height = 11)
    MFgo <- goplot(goMF)
    print(MFgo)
    dev.off()
  }
  setwd(new)
  
  ## KEGG Pathways Analysis --------------------------------------------
  if (!dir.exists("KEGG")) {  dir.create("KEGG")  }
  setwd(paste0(new,"/KEGG/"))
  k <- enrichKEGG(gene = genes,
                   organism = 'hsa',
                   pvalueCutoff = 0.2)
  write.csv(k, file=paste0(name,"_KEGG_ORA.csv"))
  
  if(dim(k) > 1) {
    pdf(paste0(name,"_KEGG_Dotplot.pdf"), width = 8, height = 11)
    KEGGdot <- dotplot(k) + ggtitle(paste0(name,"\nKEGG Dotplot"))
    print(KEGGdot)
    dev.off()
    
    k2 <- pairwise_termsim(k)
    pdf(paste0(name,"_KEGG_Map.pdf"), width = 8, height = 11)
    KEGGmap <- emapplot(k2, cex_label_category = 0.5, cex_category = 0.5, min_edge = 0.5) + ggtitle(paste0(name,"\nKEGG ORA Map"))
    print(KEGGmap)
    dev.off()
  }
  setwd(new)
  
  ## WikiPathways Analysis --------------------------------------------
  if (!dir.exists("WikiPathways")) {  dir.create("WikiPathways")  }
  setwd(paste0(new,"/WikiPathways/"))
  wiPath <- enrichWP(genes, organism = "Homo sapiens")
  write.csv(wiPath, file=paste0(name,"_WikiPathways_ORA.csv"))
  
  if(dim(wiPath) > 1) {
    pdf(paste0(name,"_WikiPathways_Dotplot.pdf"), width = 8, height = 11)
    wiDot <- dotplot(wiPath) + ggtitle(paste0(name,"\nWikiPathways ORA Dotplot"))
    print(wiDot)
    dev.off()
    
    wiPath <- pairwise_termsim(wiPath)
    pdf(paste0(name,"_WikiPathways_Map.pdf"), width = 8, height = 11)
    wiMap <- emapplot(wiPath, cex_label_category = 0.5, cex_category = 0.5, min_edge = 0.5) + ggtitle(paste0(name,"\nWikiPathways ORA Map"))
    print(wiMap)
    dev.off()
  }
  setwd(new)
  
  ## Reactome Pathways Analysis --------------------------------------------
  if (!dir.exists("Reactome")) {  dir.create("Reactome")  }
  setwd(paste0(new,"/Reactome/"))
  reactome <- enrichPathway(gene=genes, pvalueCutoff = 0.05, readable=TRUE)
  write.csv(reactome, file=paste0(name,"_Reactome_ORA.csv"))
  
  if(dim(reactome) > 1) {
    pdf(paste0(name,"_Reactome_Dotplot.pdf"), width = 8, height = 11)
    reactDot <- dotplot(reactome) + ggtitle(paste0(name,"\nReactome ORA Dotplot"))
    print(reactDot)
    dev.off()
    
    reactome <- pairwise_termsim(reactome)
    pdf(paste0(name,"_Reactome_Map.pdf"), width = 8, height = 11)
    reactMap <- emapplot(reactome, cex_label_category = 0.5, cex_category = 0.5, min_edge = 0.75) + ggtitle(paste0(name,"\nReactome ORA Map"))
    print(reactMap)
    dev.off()
  }
  setwd(current)
}