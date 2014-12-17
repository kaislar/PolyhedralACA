
#' An R implementation of the fast polyhedral adapative conjoint analysis
#' 
#' \code{polyhedralACA} implements the fast polyhedral adaptive conjoint analysis for paired comparisions
#' as described by Toubia et. al. (2003)
#' 
#'

library(MASS); library(nloptr);
library(adagio)

interiorPoint <- function(X,a){
  #' Finds an interior point in the augmented question matrix X and answer vector a.
  #' see Appendix 1 in Toubia et. al. (2003)
  #' 
  #' @param X augmented matrix describing the polyhedron/conjoint cards seend
  #' @param a user responses to each of the conjoint cards.
  #' @return an interior point inside the system of equations described by X and a.
  #' 

  qr <- dim(X)[1]; pr <- dim(X)[2];
  
  length <- function(x) {
    return(-1 * sum(x[(pr+1):(2*pr)]))
  }
  
  constraints1 <- function(x){
    h <- numeric(3*pr+1)
    h[1:pr] <- 1 - x[(pr+1):(2*pr)]
    h[(pr+1):(2*pr)] <- x[1:pr] - x[(pr+1):(2*pr)]
    h[(2*pr+1):(3*pr)] <- x[(pr+1):(2*pr)]
    h[3*pr+1] <- x[2*pr+1]- 1
    
    return(h)
  }
  
  constraints2 <- function(x){
    h <- X %*% x[1:pr] - x[(2*pr+1)]*t(a)
    return(h)
  }
  
  x0 <- runif(2*pr+1,0,1)
  
  s <- slsqp(x0, fn = length, 
             hin = constraints1,
             heq = constraints2,
             control = list(xtol_rel = 1e-16))
  return(s)
}

rescale <- function(vec){
  p <- (length(vec)-1)/2
  
  rescaled <- sapply(seq_along(vec[1:p]), function(i){ if(vec[p+1] > 0){return(vec[i]/tail(vec,n=1))}else{return(0)} })
  return( c(rescaled, vec[(p+1):length(vec)]))
}

findCenter <- function(X,a, x0){
  #' Finds an improved estimate of the polyhedron described by the matrix X of question designs and
  #' vector of answers a starting from an initial guess x0. See Appendix 1 of Toubia et. al. (2003)
  #' 
  #' @param X augmented matrix describing polyhedron/question designs
  #' @param a vector of responses for each question row in X
  #' @param x0 initial point inside the polyhedron
  #' @return center the approximate center of the polyhedron described by X and x.
  #' 

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

reshape  <- function(center, pref){
  
  vec <- matrix(0, nrow=1, ncol=length(pref))
  vec[,which(pref >0)] <- center[1:length(which(pref>0))]
  
  return(vec)  
}

majorAxis <- function(X, center, pref){
  #' Find the major axis of the ellipse described by the matrix X
  #' with specified center and set of preferred attribute levels.
  #' See Appendix 1 of Toubia, et. al. (2003)
  #' 
  #' @param X matrix describing polyhedron/paired comparision design. This has been reduced to the preferred levels
  #' @param center the approximate center of the polyhedron descried by matrix X
  #' @param pref the preferred attribute levels
  #' @returns the vector orthogonal to the longest axis of the polyhedron.
  #' 

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

resizedPolyhedron <- function( X, a){
  qr <- dim(X)[1]
  pr <- dim(X)[2]
  
  f <- function(x){ return( abs(tail(x,n=1))) }
  
  constraints1 <- function(x){
    h <- numeric(qr+pr)
    h[1:qr] <-  X %*% x - t(a)
    h[(qr+1):(pr+qr)] <- x
    return(h)
  }
  
  x0 <- runif(pr,0,1)
  
  s <- slsqp(x0, fn = f, 
             hin = constraints1,
             control = list(xtol_rel = 1e-16))
  return(s)
}

polyhedralACA <- function(X, a, pref, upperBound){
  #' Polyhedral Adaptive Conjoint algorithm. Based on previous questions (X) asked to 
  #' respondents, the responses (a), and an enumeration of the attributed from
  #' least preferrable to most preferable (pref) plus an arbitrary cutoff (upperBound), the
  #' algorithm gives the next question to ask that reduces the phase space of utilities 
  #' the fastest. See Toubia et. al (2003) paper.
  #' 
  #' @param X the matrix of question designs asked.
  #' @param a a vector of responses/ratings to each of the questions asked
  #' @param pref an ordering of the attribute levels
  #' @param uppberBound an arbitrary large number to bound the computation.
  #' @return a vector that points in the direction of the next best question.
  #' 
  
  q <- dim(X)[1]
  p <- dim(X)[2]
  
  zeros <- matrix(0, nrow=q, ncol=p)
  bigX <- rbind( cbind(X, zeros), cbind(diag(p),diag(p)))
  bigA <- cbind( a, matrix(upperBound, nrow=1, ncol=p))
  bigPref <- cbind(pref, matrix(1, nrow=1, ncol=p))
  
  reducedBigX <- bigX[,which(bigPref > 0)]
  
  c0 <- rescale( interiorPoint( reducedBigX, bigA)$par)
  
  center <- findCenter(reducedBigX, bigA, c0)
  axis <- majorAxis( X[,which(pref>0)], center, pref)
  
  analytic_center <- reshape( center, pref)
  nextQuesVector <- reshape( axis, pref)
  
  return(list(nextCard = nextQuesVector, est = analytic_center))
}

infeasibleACA <- function(X, a, pref, upperBound){
  q<- dim(X)[1]; p<- dim(X)[2];
  
  A <- rbind( cbind(X, diag(q), matrix(0, nrow=q, ncol=q)), 
              cbind(X, matrix(0, nrow=q, ncol=q), -diag(q)) )
  
  B <- cbind(diag(p), matrix(0,nrow=p, ncol=2*q))
  
  bigX <- rbind( cbind(A, matrix(0, nrow=2*q, ncol=p)),
                 cbind(B,diag(p)) )
  
  bigPref <- cbind(pref , matrix(1, nrow=1, ncol=dim(bigX)[2] - dim(pref)[2]))
  reducedBigX <- bigX[ , which(bigPref>0)]
  
  deltaX <- rbind( cbind(X, matrix(1, nrow=q, ncol=1)), 
                   cbind(-X, matrix(1, nrow=q, ncol=1)), 
                   cbind(-diag(p), matrix(0, nrow=p, ncol=1))
                )
  
  aTemp <- cbind( a, -a, matrix(-upperBound, nrow=1, ncol=p))
  deltaPref <- cbind(pref, matrix(1, nrow=1, ncol = dim(deltaX)[2]-dim(pref)[2]))
  reducedDeltaX <- deltaX[ , which(deltaPref>0)]
  
  sol <- resizedPolyhedron(reducedDeltaX, aTemp)
  delta <- tail(sol$par,n=1)
  
  bigA <- cbind( a + delta + .5, a -delta -.5, matrix(upperBound, nrow=1, ncol=p))
  
  c0 <- rescale( interiorPoint(reducedBigX, bigA)$par )
  center <- findCenter(reducedBigX, bigA, c0)
  axis <- majorAxis( X[,which(pref>0)], center, pref)

  analytic_center <- reshape( center, pref)
  nextQuesVector <- reshape( axis, pref)
  
  return(list(nextCard = nextQuesVector, est = analytic_center))
}
