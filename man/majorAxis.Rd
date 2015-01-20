\name{majorAxis}
\alias{majorAxis}
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
\value{
Returns the vector orthogonal to the longest axis of the polyhedron.
}
\examples{
}