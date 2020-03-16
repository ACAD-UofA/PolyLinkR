#' PolyLinkR Set Information.
#'
#' A data.table containing PolyLinkR Set Information. This data.table contains information
#' about PolyLink's 1,859 unique pathways created by combining pathways
#' definitions from BIOCYC (Romero et al., 2005), KEGG (Kanehisa and Goto, 2000),
#' PID (Schaefer et al., 2009), and REACTOME (Fabregat et al., 2018).
#'
#' @seealso \code{\link{PolyLinkR_SetObj}}
#' @format A data.table with 1868 rows and 3 variables:
#' \describe{
#'   \item{setID}{pathway ID}
#'   \item{setName}{pathwaty name}
#'   \item{setSource}{source database of the pathwaty}
#' }
"PolyLink_SetInfo"
