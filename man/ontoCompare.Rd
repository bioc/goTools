\name{goTools}
\alias{ontoCompare}
\alias{ontoPlot}

\title{Wrapper functions}
\description{
  This functions will allow you to describe and compare sets of oligo ids using
  Gene Ontology database
}
\usage{

ontoCompare(genelist,probeType=c("GO","hgu133a"),
goType="All", endnode, method=c("TGenes", "TIDS", "none"), plot=FALSE,
...)

ontoPlot(objM, names.arg=NULL,beside=TRUE, las=2,legend.text=TRUE, ...)
}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{genelist}{list of list of valid probe ids.}
  \item{method}{method used to evaluate the percentage of oligos for
    each end-node. 'TGenes' = for each end node, return the number of direct
  children found / total number of probe ids. (default). This includes
    oligos which do not have GO annotations. 'TIDS' = for each end node,
    return the number of direct children found / total number of GO ids
    describing the list. 'none' = for each end node, return the number of direct
    children found.}
  \item{probeType}{type of input given to the function.Valid probe types
    include GO ids and any probes ids for which a BioC annotation
    package providing a mapping to GO is available. ontoCompare is
    expecting valid probe ids.}
  \item{goType}{help sort the data by type. If 'All' (default), all oligos
    are taken into account. 'BP' restricts information
    to Biological Process, 'CC' to Cellular Component, and 'MF' to
    Molecular Function.}
  \item{plot}{logical: if 'TRUE', results are output as a graph.}
  \item{endnode}{list of GO ids corresponding to end-nodes of interest.}
  \item{beside}{Logical. If 'TRUE', the bars of the barplot are
    portrayed as juxtaposed bars. See ?barplot for more details.}
  \item{las}{numeric: if las=2, the axis labels are displayed
    perpendicular to the axis. See ?par for more details.}
  \item{legend.text}{ vector of text used to construct a legend for the
    plot. See ?barplot for more details.}
  \item{objM}{results from ontoCompare.}
  \item{names.arg}{Labels to use in ontoPlot.}
  \item{...}{extra layout parameters to be passed to ontoPlot.}

  }

  \value{

    Returns the percentage of probes children of nodes contained in
    endnode. If 'plot' = TRUE, results are plotted as a pie chart or a bargraph.
}

\author{Yee Hwa (Jean) Yang, Agnes Paquet}


\examples{
# Examples use the probeID dataset. For description type ?probeID.
# Not run

#library(GO.db)
#data(probeID)
#ontoCompare(affylist, probeType="hgu133a", plot=TRUE)
#res <- ontoCompare(operonlist["L1"], probeType="operon", method="TIDS")
#ontoPlot(res, cex=0.7)
}

\keyword{file}% at least one, from doc/KEYWORDS
\keyword{hplot}

