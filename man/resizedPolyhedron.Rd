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
}
