
\name{polyhedralACA-package}
\alias{polyhedralACA-package}
\alias{polyhedralACA}
\docType{package}
\title{Implements a polyhdral adapative conjoint estimation
}
\description{
This is an R implementation of the fast polyhedral adaptive conjoint as described by 
Toubia, et. al. in \emph{Fast Polyhedral Adaptive Conjoint Estimation}, \bold{Marketing Science} (2003). 
}
\details{
\tabular{ll}{
Package: \tab polyhedralACA\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-12-18\cr
License: \tab LGPL - 2.1\cr
LazyLoad: \tab yes\cr
}

}
\author{
B. Hsu
}
\references{
Toubia, et. al. \emph{Fast Polyhedral Adaptive Conjoint Estimation}, \bold{Marketing Science} (2003).
}
\keyword{ survey }
\seealso{
}
\examples{
X <- matrix( c(0,0,1,-1,-1,1,1,-1,0,0,0,0,0,0,0,0,0,0,0,0, 
               0,0,0,0,-1,1,0,0,0,0,1,-1,1,-1,0,0,0,0,0,0, 
               1,-1,0,0,0,0,0,0,0,0,0,0,1,-1,0,0,0,0,-1,1, 
               0,0,0,0,0,0,0,0,1,-1,-1,1,0,0,0,0,-1,1,0,0, 
               -1,1,0,0,0,0,0,0,0,0,0,0,1,-1,1,-1,0,0,0,0), nrow=5, ncol=20, byrow = TRUE)

a <- matrix( c(-68.6799, -9.4126, -57.002, 64.9769, -79.9109), nrow=1,ncol=5, byrow=TRUE)

pref <- matrix( c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1), nrow=1, ncol=20, byrow=TRUE)

upperBound <- 100

nextQuestion <- ACA(X, a, pref, upperBound)
}