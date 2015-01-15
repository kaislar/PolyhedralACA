\name{majorAxis}
\alias{majorAxis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find the major axis of an ellipse.
}
\description{
Find the major axis of the ellipse described by the matrix X
with specified center and set of preferred attribute levels.
}
\usage{
majorAxis(X, center, pref)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Matrix describing polyhedron/paired comparision design. This has been reduced to the preferred levels
}
  \item{center}{
The approximate center of the polyhedron descried by matrix X
}
  \item{pref}{
The preferred attribute levels
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns the vector orthogonal to the longest axis of the polyhedron.
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