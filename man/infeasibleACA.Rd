\name{infeasibleACA}
\alias{infeasibleACA}
\title{
Polyhedral Adaptive Conjoint algorithm for infeasible/empty polyhedrons. 
}
\description{
Based on previous questions (X) asked to 
respondents, the responses (a), and an enumeration of the attributed from
least preferrable to most preferable (pref) plus an arbitrary cutoff (upperBound), the
algorithm gives the next question to ask that reduces the phase space of utilities 
the fastest. 
}
\usage{
infeasibleACA(X, a, pref, upperBound)
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
\item{delta}{Adjustment factor needed to make the polyhedron non-empty}
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
\code{\link{fastpolyACA}}
}
\examples{
X <- matrix( c( 0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0 ,1,-1,-1,1,
                0,0,1,-1 ,-1,1,0,0,0,0,0,0,0,0,1,-1,0,0,0,0,
                0,0,0,0,0,0,-1,1,0,0,1,-1,0,0,0,0,-1,1,0,0,
                0,0,0,0,1,-1,0,0,1,-1,0,0,-1,1,0,0,0,0,0,0,
                0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,1,-1,1,-1,
                1,-1,0,0,0,0,0,0,1,-1,-1,1,0,0,0,0,0,0,0,0,
                0,0,-1,1,-1,1,0,0,1,-1,0,0,0,0,0,0,0,0,0,0,
                0,0,-1,1,-1,1,0,0,0,0,0,0,1,-1,0,0,0,0,0,0), nrow = 8, ncol=20, byrow=TRUE)

a <- matrix( c(43.8819, -2.6829, 20.7079, 59.4305, -57.5984, -17.7432, 88.7329, 25.0771), nrow=1, ncol=8,byrow=TRUE)

pref <- matrix( c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1), nrow=1, ncol=20, byrow=TRUE)

upperBound <- 100

infeasibleACA(X, a, pref, upperBound)
}
