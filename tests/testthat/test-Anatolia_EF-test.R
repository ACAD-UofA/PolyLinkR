test_that("Enrichment test and ploting works", {
  output = polylinkr(obj.info = Anatolia_EF_CLR, 
                     set.info = PolyLinkR_SetInfo,
                     set.obj = PolyLinkR_SetObj, 
                     n.cores = 4, 
                     emp.nruns = 10000, 
                     NN = 1000)
  

  genes = PolyLinkR_SetObj[setID %in% output[setP<0.001, unique(setID)]]
  gene_scores = merge(genes, Anatolia_EF_CLR, by = "objID")
  gene_scores = merge(gene_scores, PolyLinkR_SetInfo, by = "setID")
  gene_scores = merge(gene_scores, output, by = "setName")
  
  library(ggplot2)
  ggplot(gene_scores, aes(x = setName, y = objStat)) +
    geom_violin(trim = FALSE, aes(color = setName)) +
    geom_boxplot(width = 0.15, position = position_dodge(0.9), aes(color = setName)) +
    geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = .5, alpha = .35, aes(color = setName, fill = setName)) +
    scale_x_discrete(breaks=unique(gene_scores[,setName]), 
                     labels=paste0(unique(gene_scores[,setName]), "\n[p-value = ", unique(gene_scores[,setP]), "]")) +
    coord_flip() +
    labs(y="Gene-scores distribution", x = "Gene set/Pathwat", color = "Gene set/Pathway", fill = "Gene set/Pathway") +
    theme_bw()
  
})
