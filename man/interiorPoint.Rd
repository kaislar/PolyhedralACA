\name{interiorPoint}
\alias{interiorPoint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find an interior point of a polyhedron.
}
\description{
Finds an interior point in the augmented question matrix X and answer vector a.
}
\usage{
interiorPoint(X, a)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Augmented matrix describing the polyhedron/conjoint cards seend
}
  \item{a}{
User responses to each of the conjoint cards.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A vector describing an interior point inside the system of equations described by X and a.
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
