\name{goTools}
\alias{ontoCompare}
\alias{ontoPlot}

\title{Wrapper functions}
\description{
  This functions will allow you to describe and compare sets of oligo ids using
  Gene Ontology database
}
\usage{

ontoCompare(obj,  method=c("TGenes", "TIDS", "none"),
            probeType=c("GONames", "GO", "operon", "hgu133a"), goType="All",
            plot=FALSE,endnode, ...)


ontoPlot(objM, beside=TRUE, las=2,legend.text=TRUE, ...)
}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{obj}{list of valid oligos ids.}
  \item{method}{method used to evaluate the percentage of oligos for
    each end-node. 'TGenes' = number of oligo for a given
    end-node/number total number of oligos (default). This includes
    oligos not found in GO database. 'TIDS' = number of oligos for a
    given end-node/number of oligos which have annotations only. 'none'
    = percentage.}
  \item{probeType}{type of input given to the function.}
  \item{goType}{help sort the data by type. If 'ALL' (default), all oligos
    are taken into account. 'BP' restricts information
    to Biological Process, 'CC' to Cellular Component, and 'MF' to
    Molecular Function.}
  \item{plot}{logical: if 'TRUE', results are output as a graph.}
  \item{endnode}{list of GO ids corresponding to end-nodes of interest.}
  
  \item{objM}{result from ontoCompare.}
  \item{...}{extra layout parameters to be passed to ontoPlot.}
  
  }

  \value{

    Returns the percentage of probes children of nodes contained in
    endnode. If 'plot' = TRUE, results are plotted as a pie chart or a bargraph. 
}

\author{Yee Hwa (Jean) Yang, Agnes Paquet}


\examples{
# Examples use the probeID dataset. For description type ? probeID.
data(probeID)

ontoCompare(affylist, probeType="hgu133a", plot=TRUE)
res <- ontoCompare(operonlist["L1"], probeType="operon", method="TIDS")
ontoPlot(res, cex=0.7)
}

\keyword{file}% at least one, from doc/KEYWORDS
\keyword{hplot}

