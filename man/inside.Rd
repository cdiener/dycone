% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kcone.R
\name{inside}
\alias{inside}
\title{Evaluates whether a given point lies within a given k-cone. Assumes a
pointed cone (all reactions irreversible).}
\usage{
inside(x, s_matrix, m_terms, tol = sqrt(.Machine$double.eps))
}
\arguments{
\item{x}{A single point (vector) or multiple points (matrix with points as
columns) to be checked.}

\item{s_matrix}{The stochiometric matrix of the k-cone to be used.}

\item{m_terms}{The metabolic terms of the k-cone.}

\item{tol}{The numerical accuracy of the check. Defaults to the square
of the double accuracy.}
}
\value{
A single boolean of vector of booleans indicating whether the
 point(s) lie within the k-cone.
}
\description{
Evaluates whether a given point lies within a given k-cone. Assumes a
pointed cone (all reactions irreversible).
}
\examples{
S <- matrix(c(1,0,0,1,-1,0, 0, -1), nrow=2)
rownames(S) <- c('A', 'B')

# Check whether a random point falls within the flux cone
inside(runif(ncol(S)), S, rep(1,ncol(S))) # probably not true
}
\keyword{basis}
