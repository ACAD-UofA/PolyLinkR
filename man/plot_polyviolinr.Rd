% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_polyviolinr.R
\name{plot_polyviolinr}
\alias{plot_polyviolinr}
\title{Plots a violin plot showing the gene-score distribution of the enriched gene sets/pathways
above a give q-value and/or p-value cutoff(s)}
\usage{
plot_polyviolinr(
  polylinkr.out,
  set.info,
  obj.info,
  set.obj,
  p.cutoff = 0,
  q.cutoff = 0
)
}
\arguments{
\item{polylinkr.out}{data.table: output of polylinkr function}

\item{set.info}{data.frame: four required fields; gene IDs, associated scores, chromosome, and start position}

\item{obj.info}{data.frame: two required fields; genomic regions and associated scores}

\item{set.obj}{data.frame: two required fields; genomic regions and associated scores}

\item{p.cutoff}{float: desired p-value enrichement cutoff}

\item{q.cutoff}{float: desired q-value enrichement cutoff}
}
\value{
Returns a ggplot object.
}
\description{
Plots a violin plot showing the gene-score distribution of the enriched gene sets/pathways
above a give q-value and/or p-value cutoff(s)
}
\examples{
output = (obj.info = Anatolia_EF_CLR, set.info = PolyLinkR_SetInfo, 
          set.obj = PolyLinkR_SetObj, n.cores = 8, emp.nruns = 1000, NN = 10000)
p = plot_polyviolinr(polylinkr.out =  output, obj.info = Anatolia_EF_CLR, 
                     set.info = PolyLinkR_SetInfo, set.obj = PolyLinkR_SetObj, p.cutoff = 0.001)

}
\seealso{
\code{\link{polylinkr}}, \code{\link{PolyLinkR_SetInfo}}, \code{\link{PolyLinkR_SetObj}}, \code{\link{Anatolia_EF_CLR}}
}
