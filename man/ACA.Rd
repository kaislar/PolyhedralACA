\name{ACA}
\alias{ACA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Wrapper function of the polyhedral adapative conjoint
}
\description{
Top level function that runs the polyhedral adaptive conjoint.
}
\usage{
ACA(X, A, pref, upperBound, delta = NULL)
}
%- maybe also 'usage' for other objects documented here.
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
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list containing:
\item{nextCard}{A vector pointing along the next question vector}
\item{est}{The estimate of the utilities based on the questions and responses}
\item{delta}{Adjustment factor needed. NULL if \code{X}, \code{a} described a full polyhedron.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{infesasibleACA}}, \code{\link{polyhedralACA}}
}
\examples{

}
