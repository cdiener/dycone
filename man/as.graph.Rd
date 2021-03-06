% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_utility.R
\name{as.graph}
\alias{as.graph}
\title{Converts a reaction list into a graph object.}
\usage{
as.graph(reacts)
}
\arguments{
\item{reacts}{A reaction list.}

\item{...}{other arguments passed to specific methods}
}
\value{
An igraph object representing with the species being the nodes and
 reactions denoting edges.
}
\description{
Converts a reaction list into a graph object.
}
\examples{
data(eryth)
print(as.graph(eryth))
}
\keyword{graph}
\keyword{reactions,}
