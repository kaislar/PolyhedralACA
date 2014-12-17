library(MASS); library(nloptr);
library(adagio)


X_empty <- matrix( c( 0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0 ,1,-1,-1,1,
                      0,0,1,-1 ,-1,1,0,0,0,0,0,0,0,0,1,-1,0,0,0,0,
                      0,0,0,0,0,0,-1,1,0,0,1,-1,0,0,0,0,-1,1,0,0,
                      0,0,0,0,1,-1,0,0,1,-1,0,0,-1,1,0,0,0,0,0,0,
                      0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,1,-1,1,-1,
                      1,-1,0,0,0,0,0,0,1,-1,-1,1,0,0,0,0,0,0,0,0,
                      0,0,-1,1,-1,1,0,0,1,-1,0,0,0,0,0,0,0,0,0,0,
                      0,0,-1,1,-1,1,0,0,0,0,0,0,1,-1,0,0,0,0,0,0), nrow = 8, ncol=20, byrow=TRUE)

a_empty <- matrix( c(43.8819, -2.6829, 20.7079, 59.4305, -57.5984, -17.7432, 88.7329, 25.0771), nrow=1, ncol=8,byrow=TRUE)

pref_example <- matrix( c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1), nrow=1, ncol=20, byrow=TRUE)

upperBound_example <- 100

infeasibleACA(X_empty, a_empty, pref_example, upperBound_example)