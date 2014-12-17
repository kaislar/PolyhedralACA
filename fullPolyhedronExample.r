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


q <- dim(X_example)[1]
p <- dim(X_example)[2]

zeros <- matrix(0,nrow=q,ncol=p)
bigX <- rbind( cbind(X_example, zeros), cbind(diag(p),diag(p)))
bigA <- cbind( a_example, matrix(upperBound_example, nrow=1, ncol=p))
bigPref <- cbind(pref_example, matrix(1,nrow=1,ncol=p))

reducedBigX <- bigX[,which(bigPref>0)]
reducedX <- X_example[,which(pref_example>0)]

qr <- dim(reducedBigX)[1]
pr <- dim(reducedBigX)[2]



