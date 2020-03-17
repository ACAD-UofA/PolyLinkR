# Ploting function to generate violoin plots of object (gene) score 
# distribution of enriched pathways above a user-specified cutoff


#violin plotting function
#' Plots a violin plot showing the gene-score distribution of the enriched gene sets/pathways
#' above a give q-value and/or p-value cutoff(s)
#' 
#'
#' @param polylinkr.out data.table: output of polylinkr function
#' @param set.info data.frame: four required fields; gene IDs, associated scores, chromosome, and start position
#' @param obj.info data.frame: two required fields; genomic regions and associated scores
#' @param set.obj data.frame: two required fields; genomic regions and associated scores
#' @param p.cutoff float: desired p-value enrichement cutoff
#' @param q.cutoff float: desired q-value enrichement cutoff
#'
#' @return Returns a ggplot object.
#'
#' @seealso \code{\link{polylinkr}}, \code{\link{PolyLinkR_SetInfo}}, \code{\link{PolyLinkR_SetObj}}, \code{\link{Anatolia_EF_CLR}} 
#'
#'
#' @examples
#' output = (obj.info = Anatolia_EF_CLR, set.info = PolyLinkR_SetInfo, set.obj = PolyLinkR_SetObj, n.cores = 8, emp.nruns = 1000, NN = 10000)
#' p = plot_polyviolinr(polylinkr.out =  output, obj.info = Anatolia_EF_CLR, set.info = PolyLinkR_SetInfo, set.obj = PolyLinkR_SetObj, p.cutoff = 0.001)
#'
#' @export
plot_polyviolinr <- function(polylinkr.out, set.info, obj.info, set.obj, p.cutoff = 0, q.cutoff = 0){
  
  # Get the Gene-scores in the enriched pathways bellow or equal a given q or p-value cutoff
  genes = PolyLinkR_SetObj[setID %in% output[setP <= p.cutoff & setQ <= q.cutoff, unique(setID)]]
  gene_scores = merge(genes, Anatolia_EF_CLR, by = "objID")
  gene_scores = merge(gene_scores, PolyLinkR_SetInfo, by = "setID")
  gene_scores = merge(gene_scores, output, by = "setName")
  
  library(ggplot2)
  plot <- ggplot(gene_scores, aes(x = setName, y = objStat)) +
    geom_violin(trim = FALSE, aes(color = setName)) +
    geom_boxplot(width = 0.15, position = position_dodge(0.9), aes(color = setName)) +
    geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = .5, alpha = .35, aes(color = setName, fill = setName)) +
    scale_x_discrete(breaks=unique(gene_scores[,setName]), 
                     labels=paste0(unique(gene_scores[,setName]), "\n[p-value = ", unique(gene_scores[,setP]), "]")) +
    coord_flip() +
    labs(y="Gene-scores distribution", x = "Gene set/Pathwat", color = "Gene set/Pathway", fill = "Gene set/Pathway") +
    theme_bw() 
  
  return(plot)
}
