\name{findCenter}
\alias{findCenter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
findCenter(X, a, x0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{a}{
%%     ~~Describe \code{a} here~~
}
  \item{x0}{
%%     ~~Describe \code{x0} here~~
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
function (X, a, x0) 
{
    qr <- dim(X)[1]
    pr <- dim(X)[2]
    u <- x0[1:pr]
    matU <- diag(u)
    matU2 <- matU \%*\% matU
    Uinv <- ginv(matU)
    z <- ginv(X \%*\% (matU2 \%*\% t(X))) \%*\% t(a)
    d <- u - (matU2 \%*\% t(X)) \%*\% z
    if (t(Uinv \%*\% d) \%*\% (Uinv \%*\% d) < 0.25) {
        center <- u + d
    }
    else {
        f <- function(x) return(-sum(log(x)))
        gradf <- function(x) return(-1/x)
        step <- linesearch_ww(u, d, fn = f, gr = gradf)
        center <- step$xalpha
    }
    return(center)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
