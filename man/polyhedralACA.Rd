\name{polyhedralACA}
\alias{polyhedralACA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
polyhedralACA(X, a, pref, upperBound)
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
