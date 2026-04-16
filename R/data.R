#' @name anatolia_ef_clr
#' @title anatolia_ef_clr
#' @description Example data: Gene expression factors from Anatolia population.
#' @format A data frame with gene expression values.
#' @source Internal example dataset for testing.
NULL

#' @name polylinkr_set_info
#' @title polylinkr_set_info
#' @description Example gene set information for testing.
#' @format A data frame with setID, setName, setSource columns.
#' @source Internal example dataset for testing.
NULL

#' @name polylinkr_set_obj
#' @title polylinkr_set_obj
#' @description Example gene-to-set mapping for testing.
#' @format A data frame with setID and objID columns.
#' @source Internal example dataset for testing.
NULL

#' @name rr.dt
#' @title rr.dt
#' @description Recombination rate data for testing.
#' @format A data frame with recombination rates.
#' @source Internal example dataset for testing.
NULL

#' @name human_sexavg_rr_bherer2017
#' @title human_sexavg_rr_bherer2017
#' @description Genome-wide recombination rate data from Bherer et al. 2017.
#'   Sex-averaged recombination rates for GRCh37/hg19 across the human genome.
#' @format A data.table with columns:
#'   \describe{
#'     \item{chr}{Chromosome number}
#'     \item{pos}{Base pair position}
#'     \item{rate}{Recombination rate (cM/Mb)}
#'     \item{cM}{Cumulative genetic map position (centiMorgans)}
#'   }
#' @source Bherer et al. (2017) - refined maps from pedigree and population data
NULL

# Deprecated data object aliases for backward compatibility
#' @name Anatolia_EF_CLR
#' @title Anatolia_EF_CLR (deprecated)
#' @description Deprecated alias for \code{anatolia_ef_clr}.
#' @format See \code{anatolia_ef_clr}.
#' @keywords internal
NULL

#' @name PolyLinkR_SetInfo
#' @title PolyLinkR_SetInfo (deprecated)
#' @description Deprecated alias for \code{polylinkr_set_info}.
#' @format See \code{polylinkr_set_info}.
#' @keywords internal
NULL

#' @name PolyLinkR_SetObj
#' @title PolyLinkR_SetObj (deprecated)
#' @description Deprecated alias for \code{polylinkr_set_obj}.
#' @format See \code{polylinkr_set_obj}.
#' @keywords internal
NULL

#' @name human_sexavg_RR_Bherer2017
#' @title human_sexavg_RR_Bherer2017 (deprecated)
#' @description Deprecated alias for \code{human_sexavg_rr_bherer2017}.
#' @format See \code{human_sexavg_rr_bherer2017}.
#' @keywords internal
NULL