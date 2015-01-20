\name{resizedPolyhedron}
\alias{resizedPolyhedron}
\title{Resize the polyhedron describe by X and a to be non-empty
}
\description{
Create a non-empty polyhedron X and a + delta by finding the smallest delta so that the system is non-empty.
}
\usage{
resizedPolyhedron(X, a)
}
\arguments{
  \item{X}{
Matrix of question design vectors describing an empty/infeasible polyhedron
}
  \item{a}{
Vector of responses to each conjoint question
}
}
\value{
Returns a minimum adjustment factor delta to adjust the polyhedron X,a so that it is non-empty
}
\examples{
}
