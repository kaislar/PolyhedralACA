\name{ACA}
\alias{ACA}
\title{
Wrapper function of the polyhedral adapative conjoint
}
\description{
Top level function that runs the polyhedral adaptive conjoint.
}
\usage{
ACA(X, A, pref, upperBound, delta = NULL)
}
\arguments{
  \item{X}{The matrix of question designs asked
}
  \item{A}{A vector of responses/ratings to each of the questions asked
}
  \item{pref}{An ordering of the attribute levels
}
  \item{upperBound}{An arbitrary large number to bound the computation.
}
  \item{delta}{Adjustment factor to make the polyhedron described by \code{X}, \code{a} non-empty. NULL by default
}
}
\value{
A list containing:
\item{nextCard}{A vector pointing along the next question vector}
\item{est}{The estimate of the utilities based on the questions and responses}
\item{delta}{Adjustment factor needed. NULL if \code{X}, \code{a} described a full polyhedron.}
}

\seealso{
\code{\link{infesasibleACA}}, \code{\link{fastpolyACA}}
}