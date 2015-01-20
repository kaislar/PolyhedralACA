\name{fastpolyACA}
\alias{fastpolyACA}
\title{
Polyhedral Adaptive Conjoint algorithm. 
}
\description{
Based on previous questions (X) asked to 
respondents, the responses (a), and an enumeration of the attributed from
least preferrable to most preferable (pref) plus an arbitrary cutoff (upperBound), the
algorithm gives the next question to ask that reduces the phase space of utilities 
the fastest.
}
\usage{
fastpolyACA(X, a, pref, upperBound)
}
\arguments{
  \item{X}{
The matrix of question designs asked.
}
  \item{a}{
Vector of responses/ratings to each of the questions asked
}
  \item{pref}{
An ordering of the attribute levels
}
  \item{upperBound}{
An arbitrary large number to bound the computation.
}
}
\details{
}
\value{
Returns a list containing:
\item{nextCard}{Vector that points in the direction of the next question}
\item{est}{Estimate of the utilities for each of the attributes given the questions and responses}
\item{delta}{NULL. Adjustment factor needed if polyhedron is empty. NULL in this case}
}
\references{
Toubia, Simester, Hauser, Dahan, \emph{Fast Polyhedral Adaptive Conjoint Estimation}, \bold{Marketing Science} (2003)
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{infeasibleACA}}
}
\examples{
X <- matrix( c(0,0,1,-1,-1,1,1,-1,0,0,0,0,0,0,0,0,0,0,0,0, 
               0,0,0,0,-1,1,0,0,0,0,1,-1,1,-1,0,0,0,0,0,0, 
               1,-1,0,0,0,0,0,0,0,0,0,0,1,-1,0,0,0,0,-1,1, 
               0,0,0,0,0,0,0,0,1,-1,-1,1,0,0,0,0,-1,1,0,0, 
               -1,1,0,0,0,0,0,0,0,0,0,0,1,-1,1,-1,0,0,0,0), nrow=5, ncol=20, byrow = TRUE)

a <- matrix( c(-68.6799, -9.4126, -57.002, 64.9769, -79.9109), nrow=1,ncol=5, byrow=TRUE)

pref <- matrix( c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1), nrow=1, ncol=20, byrow=TRUE)

upperBound <- 100

fastpolyACA(X, a, pref, upperBound)
}       