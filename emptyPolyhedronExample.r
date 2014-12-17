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


q<- dim(X_empty)[1]
p<- dim(X_empty)[2]

A <- rbind( cbind(X_empty, diag(q),matrix(0,nrow=q,ncol=q)), 
            cbind(X_empty, matrix(0,nrow=q,ncol=q),-diag(q)) )

B <- cbind(diag(p), matrix(0,nrow=p, ncol=2*q))

bigX_empty<- rbind( cbind(A,matrix(0,nrow=2*q,ncol=p)), cbind(B,diag(p)) )
bigPref <- cbind(pref_example, matrix(1,nrow=1, ncol=dim(bigX_empty)[2] - dim(pref_example)[2]))
reducedBigX_empty <- bigX_empty[,which(bigPref>0)]

deltaX <- rbind( cbind(X_empty,matrix(1,nrow=q,ncol=1)), cbind(-X_empty, matrix(1,nrow=q,ncol=1)), cbind(-diag(p), matrix(0,nrow=p, ncol=1)))

aTemp <- cbind(a_empty,-a_empty, matrix(-upperBound_example,nrow=1,ncol=p))
deltaPref <- cbind(pref_example, matrix(1,nrow=1,ncol = dim(deltaX)[2]-dim(pref_example)[2]))
reducedDeltaX <- deltaX[,which(deltaPref>0)]

sol <- resizedPolyhedron(reducedDeltaX,aTemp)

delta <- tail(sol$par,n=1)
bigA <- cbind(a_empty+delta+.5, a_empty-delta-.5, matrix(upperBound_example,nrow=1,ncol=p))
u <- sol$par[1:length(sol$par)-1]

qr <- dim(reducedBigX_empty)[1]
pr <- dim(reducedBigX_empty)[2]

u0 <- rescale( interiorPoint(reducedBigX_empty, bigA)$par )
center <- findCenter(reducedBigX_empty, bigA, u0)
