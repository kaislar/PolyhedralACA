\name{reshape}
\alias{reshape}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
Reshapes the a vector to match attribute level preferences
}
\usage{
reshape(originalVector, pref)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{originalVector}{
Vector to be reshaped
}
  \item{pref}{
Attribute level preferences. 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
originalVector reshaped to match non-zero levels of pref.
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
function (originalVector, pref) 
{
    vec <- matrix(0, nrow = 1, ncol = length(pref))
    vec[, which(pref > 0)] <- originalVector[1:length(which(pref > 
        0))]
    return(vec)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
