\name{infeasibleACA}
\alias{infeasibleACA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
infeasibleACA(X, a, pref, upperBound)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{a}{
%%     ~~Describe \code{a} here~~
}
  \item{pref}{
%%     ~~Describe \code{pref} here~~
}
  \item{upperBound}{
%%     ~~Describe \code{upperBound} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

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
