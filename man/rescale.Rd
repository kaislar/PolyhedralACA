\name{rescale}
\alias{rescale}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
Rescales an approximate solution to the center of a polyhedron
}
\usage{
rescale(vec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vec}{
Vector describing the approximate center of the polyehdron
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Approximation of the center that has been rescaled
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
function (vec) 
{
    p <- (length(vec) - 1)/2
    rescaled <- sapply(seq_along(vec[1:p]), function(i) {
        if (vec[p + 1] > 0) {
            return(vec[i]/tail(vec, n = 1))
        }
        else {
            return(0)
        }
    })
    return(c(rescaled, vec[(p + 1):length(vec)]))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
