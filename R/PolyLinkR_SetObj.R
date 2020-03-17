#' PolyLinkR Set Object.
#'
#' A data.table containing PolyLinkR Set Object. This data.table links the gene number
#' to the corresponding PolyLink 1,868 unique pathways created by combining pathways
#' definitions from BIOCYC (Romero et al., 2005), KEGG (Kanehisa and Goto, 2000),
#' PID (Schaefer et al., 2009), and REACTOME (Fabregat et al., 2018).
#'
#' @seealso \code{\link{polylinkr}},\code{\link{PolyLinkR_SetInfo}}, \code{\link{Anatolia_EF_CLR}}
#' @format A data.table with 132382 rows and 2 variables:
#' \describe{
#'   \item{setID}{pathway ID}
#'   \item{objID}{gene ID}
#' }
"PolyLink_SetObj"
