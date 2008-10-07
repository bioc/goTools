\name{goTools}
\alias{EndNodeList}
\alias{CustomEndNodeList}

\title{Reference GO nodes list.}
\description{
  The function \code{EndNodeList} builds the default end node list used in
  \code{\link{ontoCompare}}. \code{CustomEndNodeList} helps you build a list
    of children nodes starting from a GO id of interest. 
}
\usage{
EndNodeList()
CustomEndNodeList(id,rank=1)
}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{id}{Valid GO id: "GO:XXXXXXX".}
  \item{rank}{Number of levels of children of the GO DAG under id you
    want to add to your nodes list.}
}

\value{
  \code{EndNodeList} returns a vector of all GO ids 1 level below MF("GO:0003674"),
  BP("GO:0008150") and CC("GO:0005575"). MF, BP and CC nodes are
  included.
  \code{CustomEndNodeList} returns a vector of all GO ids children of
  \code{id}, \code{rank} levels below it. 
}

\author{Yee Hwa (Jean) Yang, Agnes Paquet}

\seealso{\code{\link{ontoCompare}}}

\examples{
# Examples use the probeID dataset. For description type ? probeID.
library(GO.db)
EndNodeList()
MFendnode <- CustomEndNodeList("GO:0003674", rank=2)

## Example (not run)
## data(probeID)
## ontoCompare(affylist, probeType="hgu133a", endnode=MFendnode, goType="MF")

}

\keyword{file}% at least one, from doc/KEYWORDS
\keyword{manip}

