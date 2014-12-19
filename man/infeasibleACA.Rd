\name{infeasibleACA}
\alias{infeasibleACA}
%- Also NEED an '\alias' for EACH other topic documented here.
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
\item{delta}{Adjustment factor needed to make the polyhedron non-empty}
%% ...
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
\code{\link{polyhedralACA}}
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

## The function is currently defined as
function (X, a, pref, upperBound) 
{
    q <- dim(X)[1]
    p <- dim(X)[2]
    A <- rbind(cbind(X, diag(q), matrix(0, nrow = q, ncol = q)), 
        cbind(X, matrix(0, nrow = q, ncol = q), -diag(q)))
    B <- cbind(diag(p), matrix(0, nrow = p, ncol = 2 * q))
    bigX <- rbind(cbind(A, matrix(0, nrow = 2 * q, ncol = p)), 
        cbind(B, diag(p)))
    bigPref <- cbind(pref, matrix(1, nrow = 1, ncol = dim(bigX)[2] - 
        dim(pref)[2]))
    reducedBigX <- bigX[, which(bigPref > 0)]
    deltaX <- rbind(cbind(X, matrix(1, nrow = q, ncol = 1)), 
        cbind(-X, matrix(1, nrow = q, ncol = 1)), cbind(-diag(p), 
            matrix(0, nrow = p, ncol = 1)))
    aTemp <- cbind(a, -a, matrix(-upperBound, nrow = 1, ncol = p))
    deltaPref <- cbind(pref, matrix(1, nrow = 1, ncol = dim(deltaX)[2] - 
        dim(pref)[2]))
    reducedDeltaX <- deltaX[, which(deltaPref > 0)]
    sol <- resizedPolyhedron(reducedDeltaX, aTemp)
    delta <- tail(sol$par, n = 1)
    bigA <- cbind(a + delta + 0.5, a - delta - 0.5, matrix(upperBound, 
        nrow = 1, ncol = p))
    c0 <- rescale(interiorPoint(reducedBigX, bigA)$par)
    center <- findCenter(reducedBigX, bigA, c0)
    axis <- majorAxis(X[, which(pref > 0)], center, pref)
    analytic_center <- reshape(center, pref)
    nextQuesVector <- reshape(axis, pref)
    return(list(nextCard = Re(nextQuesVector), est = Re(analytic_center), 
        delta = delta))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
