PolyhedralACA
=============

An R implementation of the fast polyhedral adaptive conjoint analysis. 

This is an R implementation of the fast polyhedral adaptive conjoint as described by Toubia, et. al. in Fast
Polyhedral Adaptive Conjoint Estimation, Marketing Science (2003). 

These codes are free software; you can redistribute then and/or modify them under the terms of the GNU Lesser
General Public License as published by the Free Software Foundation; either version 2.1 of the license, or (at
your option) any later version. 

These codes are distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the
implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

Example
-------

```RDoc
X <- matrix( 
c(0,0,1,-1,-1,1,1,-1,0,0,0,0,0,0,0,0,0,0,0,0, 
0,0,0,0,-1,1,0,0,0,0,1,-1,1,-1,0,0,0,0,0,0, 
1,-1,0,0,0,0,0,0,0,0,0,0,1,-1,0,0,0,0,-1,1, 
0,0,0,0,0,0,0,0,1,-1,-1,1,0,0,0,0,-1,1,0,0, 
-1,1,0,0,0,0,0,0,0,0,0,0,1,-1,1,-1,0,0,0,0), nrow=5, ncol=20, byrow = TRUE)

a <- matrix( c(-68.6799, -9.4126, -57.002, 64.9769, -79.9109), nrow=1,ncol=5, byrow=TRUE)

pref <- matrix( c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1), nrow=1, ncol=20, byrow=TRUE)

upperBound <- 100

nextQuestion <- ACA(X, a, pref, upperBound)

```

Depends 
-------

R (>= 1.8.0), nloptr, MASS, adagio

License 
-------

LGPL - 2.1