#' SweepFinder CLR scores for Anatolian Early Farmers.
#'
#' A data.table containing the log-transformed gene-length standardised
#' CLR scores assigned to genes generated using SweepFinder2 (Huber et al., 2016)
#' ran on the 1240k SNPs dataset of the Anatolian Early Farmers.
#'
#' @format A data.table with 11023 rows and 8 variables:
#' \describe{
#' \item{objID: gene number in the PolyLinkR pathwat set}
#' \item{objStat: log-transformed gene-length standardised CLR score}
#' \item{objName: HGNC gene name}
#' \item{GeneLength: gene length}
#' \item{chr: chromosome}
#' \item{startpos: start position of the gene}
#' \item{endpos: end position of the gene}
#' \item{strand: strand the gene is encoded on}
#' }
"Anatolia_EF_CLR"
