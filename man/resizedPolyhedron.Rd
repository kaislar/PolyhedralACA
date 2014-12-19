\name{resizedPolyhedron}
\alias{resizedPolyhedron}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Resize the polyhedron describe by X and a to be non-empty
}
\description{
Create a non-empty polyhedron X and a + delta by finding the smallest delta so that the system is non-empty.
}
\usage{
resizedPolyhedron(X, a)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Matrix of question design vectors describing an empty/infeasible polyhedron
}
  \item{a}{
Vector of responses to each conjoint question
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a minimum adjustment factor delta to adjust the polyhedron X,a so that it is non-empty
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
    f <- function(x) {
        return(abs(tail(x, n = 1)))
    }
    constraints1 <- function(x) {
        h <- numeric(qr + pr)
        h[1:qr] <- X \%*\% x - t(a)
        h[(qr + 1):(pr + qr)] <- x
        return(h)
    }
    x0 <- runif(pr, 0, 1)
    s <- slsqp(x0, fn = f, hin = constraints1, control = list(xtol_rel = 1e-16))
    return(s)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
