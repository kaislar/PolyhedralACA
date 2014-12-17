library(MASS); library(nloptr);
library(adagio)

X_example <- matrix( c(0,0,1,-1,-1,1,1,-1,0,0,0,0,0,0,0,0,0,0,0,0, 
                       0,0,0,0,-1,1,0,0,0,0,1,-1,1,-1,0,0,0,0,0,0, 
                       1,-1,0,0,0,0,0,0,0,0,0,0,1,-1,0,0,0,0,-1,1, 
                       0,0,0,0,0,0,0,0,1,-1,-1,1,0,0,0,0,-1,1,0,0, 
                       -1,1,0,0,0,0,0,0,0,0,0,0,1,-1,1,-1,0,0,0,0), nrow=5, ncol=20, byrow = TRUE)

a_example <- matrix( c(-68.6799, -9.4126, -57.002, 64.9769, -79.9109), nrow=1,ncol=5, byrow=TRUE)

pref_example <- matrix( c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1), nrow=1, ncol=20, byrow=TRUE)

upperBound_example <- 100

polyhedralACA(X_example, a_example, pref_example, upperBound_example)