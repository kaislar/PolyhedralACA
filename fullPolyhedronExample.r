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


interiorPoint <- function(X,a){
  qr <- dim(X)[1]; pr <- dim(X)[2];
  
  length <- function(x) {
    return(-1 * sum(x[(pr+1):(2*pr)]))
  }
  
  constraints1 <- function(x){
    h <- numeric(3*pr+1)
    h[1:pr] <- 1 - x[(pr+1):(2*pr)]
    h[(pr+1):(2*pr)] <- x[1:pr] - x[(pr+1):(2*pr)]
    h[(2*pr+1):(3*pr)] <- x[(pr+1):(2*pr)]
    h[3*pr+1] <- x[2*pr]- 1
    
    return(h)
  }
  
  constraints2 <- function(x){
    h <- X %*% x[1:pr] - x[2*pr]*t(a)
    return(h)
  }
  
  x0 <- runif(2*pr,0,1)
  
  s <- slsqp(x0, fn = length, 
             hin = constraints1,
             heq = constraints2,
             control = list(xtol_rel = 1e-16))
  return(s)
}

findCenter <- function(X,a, x0){
  qr <- dim(X)[1]; pr <- dim(X)[2];
  
  u <- x0[1:pr]
  matU <- diag(u)
  matU2 <- matU %*% matU
  Uinv <- ginv(matU)
  
  z <- ginv(X %*% (matU2 %*% t(X))) %*% t(a)
  d <- u - (matU2 %*% t(X)) %*% z
  
  if( t(Uinv %*% d) %*% (Uinv %*% d) <0.25){
    center <- u + d
  }
  else{
    f <-function(x) return(-sum(log(x)))
    gradf <- function(x) return(-1/x)
    
    step <- linesearch_ww(u, d, fn = f, gr=gradf)
    center <- step$xalpha
  }
  
  return(center)  
}

majorAxis <- function(X, center, pref){
  qr <- dim(X)[1]; pr <- dim(X)[2]
  
  u <- center[1:pr]
  u_complement <- center[(pr+1):(3*pr)][which(pref>0)]
  U2 <- diag( 1/u**2 + 1/u_complement**2)

  z <- ginv( X %*% t(X))
  XNormed <- (t(X) %*% z) %*% X
  M <- U2 - XNormed %*% U2
  
  eigenSys <- eigen(M)
  minIndx <- which.min( eigenSys$values[which(eigenSys$values>1e-12)])
  
  return(eigenSys$vectors[,minIndx])
}


