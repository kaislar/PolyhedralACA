\name{majorAxis}
\alias{majorAxis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
majorAxis(X, center, pref)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{center}{
%%     ~~Describe \code{center} here~~
}
  \item{pref}{
%%     ~~Describe \code{pref} here~~
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
function (X, center, pref) 
{
    qr <- dim(X)[1]
    pr <- dim(X)[2]
    u <- center[1:pr]
    u_complement <- center[(pr + 1):(3 * pr)][which(pref > 0)]
    U2 <- diag(1/u^2 + 1/u_complement^2)
    z <- ginv(X \%*\% t(X))
    XNormed <- (t(X) \%*\% z) \%*\% X
    M <- U2 - XNormed \%*\% U2
    eigenSys <- eigen(M)
    minIndx <- which.min(eigenSys$values[which(Re(eigenSys$values) > 
        1e-12)])
    return(eigenSys$vectors[, minIndx])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
