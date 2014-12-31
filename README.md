PolyhedralACA
=============

An R implementation of the fast polyhedral adaptive conjoint analysis. 

This is an R implementation of the fast polyhedral adaptive conjoint as described by Toubia, et. al. in Fast
Polyhedral Adaptive Conjoint Estimation, Marketing Science (2003). 

Example
-------

```R
X <- matrix( c(0,0,1,-1,-1,1,1,-1,0,0,0,0,0,0,0,0,0,0,0,0, 
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