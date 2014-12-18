\name{interiorPoint}
\alias{interiorPoint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
interiorPoint(X, a)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{a}{
%%     ~~Describe \code{a} here~~
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
function (X, a) 
{
    qr <- dim(X)[1]
    pr <- dim(X)[2]
    length <- function(x) {
        return(-1 * sum(x[(pr + 1):(2 * pr)]))
    }
    constraints1 <- function(x) {
        h <- numeric(3 * pr + 1)
        h[1:pr] <- 1 - x[(pr + 1):(2 * pr)]
        h[(pr + 1):(2 * pr)] <- x[1:pr] - x[(pr + 1):(2 * pr)]
        h[(2 * pr + 1):(3 * pr)] <- x[(pr + 1):(2 * pr)]
        h[3 * pr + 1] <- x[2 * pr + 1] - 1
        return(h)
    }
    constraints2 <- function(x) {
        h <- X \%*\% x[1:pr] - x[(2 * pr + 1)] * t(a)
        return(h)
    }
    x0 <- runif(2 * pr + 1, 0, 1)
    s <- slsqp(x0, fn = length, hin = constraints1, heq = constraints2, 
        control = list(xtol_rel = 1e-16))
    return(s)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
