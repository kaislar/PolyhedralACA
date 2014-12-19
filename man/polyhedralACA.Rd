\name{polyhedralACA}
\alias{polyhedralACA}
%- Also NEED an '\alias' for EACH other topic documented here.
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
polyhedralACA(X, a, pref, upperBound)
}
%- maybe also 'usage' for other objects documented here.
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
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a list containing:
\item{nextCard}{Vector that points in the direction of the next question}
\item{est}{Estimate of the utilities for each of the attributes given the questions and responses}
\item{delta}{NULL. Adjustment factor needed if polyhedron is empty. NULL in this case}
}
\references{
Toubia, Simester, Hauser, Dahan, \it{Fast Polyhedral Adaptive Conjoint Estimation}, \bf{Marketing Science} (2003)
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

polyhedralACA(X, a, pref, upperBound)

## The function is currently defined as
function (X, a, pref, upperBound) 
{
    q <- dim(X)[1]
    p <- dim(X)[2]
    zeros <- matrix(0, nrow = q, ncol = p)
    bigX <- rbind(cbind(X, zeros), cbind(diag(p), diag(p)))
    bigA <- cbind(a, matrix(upperBound, nrow = 1, ncol = p))
    bigPref <- cbind(pref, matrix(1, nrow = 1, ncol = p))
    reducedBigX <- bigX[, which(bigPref > 0)]
    c0 <- rescale(interiorPoint(reducedBigX, bigA)$par)
    center <- findCenter(reducedBigX, bigA, c0)
    axis <- majorAxis(X[, which(pref > 0)], center, pref)
    analytic_center <- reshape(center, pref)
    nextQuesVector <- reshape(axis, pref)
    return(list(nextCard = Re(nextQuesVector), est = Re(analytic_center), 
        delta = NULL))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
