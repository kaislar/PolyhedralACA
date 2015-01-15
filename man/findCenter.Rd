\name{findCenter}
\alias{findCenter}
\title{Finds center of a polyhedron.
}
\description{
Finds an improved estimate of the polyhedron described by the matrix X of question designs and
vector of answers a starting from an initial guess x0. 
}
\usage{
findCenter(X, a, x0)
}

\arguments{
  \item{X}{
Augmented matrix describing polyhedron/question designs
}
  \item{a}{
Vector of responses for each question row in X
}
  \item{x0}{
Initial point inside the polyhedron
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns the approximate center of the polyhedron described by X and a.
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
\examples{}
