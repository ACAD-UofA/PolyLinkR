#' Anatolia gene expression example data
#'
#' Example data: Gene expression factors from Anatolia population.
#'
#' @docType data
#' @name anatolia_ef_clr
#' @usage data(anatolia_ef_clr)
#' @format A data frame with gene expression values.
#' @source Internal example dataset for testing.
#' @keywords datasets
NULL

#' PolylinkR example gene set information
#'
#' Example gene set information for testing.
#'
#' @docType data
#' @name polylinkr_set_info
#' @usage data(polylinkr_set_info)
#' @format A data frame with setID, setName, setSource columns.
#' @source Internal example dataset for testing.
#' @keywords datasets
NULL

#' PolylinkR example gene-to-set mapping
#'
#' Example gene-to-set mapping for testing.
#'
#' @docType data
#' @name polylinkr_set_obj
#' @usage data(polylinkr_set_obj)
#' @format A data frame with setID and objID columns.
#' @source Internal example dataset for testing.
#' @keywords datasets
NULL

#' @name rr.dt
#' @title rr.dt
#' @description Recombination rate data for testing.
#' @format A data frame with recombination rates.
#' @source Internal example dataset for testing.
NULL

#' Human recombination rate data from Bherer et al. 2017
#'
#' Genome-wide recombination rate data from Bherer et al. 2017.
#' Sex-averaged recombination rates for GRCh37/hg19 across the human genome.
#'
#' @docType data
#' @name human_sexavg_rr_bherer2017
#' @usage data(human_sexavg_rr_bherer2017)
#' @format A data.table with columns:
#' \describe{
#'   \item{chr}{Chromosome number}
#'   \item{pos}{Base pair position}
#'   \item{rate}{Recombination rate (cM/Mb)}
#'   \item{cM}{Cumulative genetic map position (centiMorgans)}
#' }
#' @source Bherer et al. (2017) - refined maps from pedigree and population data
#' @keywords datasets
NULL