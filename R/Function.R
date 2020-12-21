########################################
## Loading needed packages
require(foreach)
require(doParallel)
require(MASS)
require(nnet)
require(glmnet)
require(VGAM)

##############################################
## Compute negative log likelihood value
f <- function(Y, X, C) {
  n = dim(X)[1]
  z <- X%*%C
  z <- exp(z)
  z <- log(apply(z,1,sum)+1)
  loglike <- sum(diag(X%*%C%*%t(Y[,-1])))-sum(z)
  return(-loglike/n)
}


## Compute probability matrix P
exception1 <- function(X, C) {
  z <- X%*%C
  z <- exp(z)
  z1 <- as.matrix(z[,-1])
  d <- apply(z1,1,sum)+1
  p1 <- as.matrix(1/d)
  p2 <- z1/d
  P <- cbind(p1,p2)
  return(P)
}

exception <- function(X, B, V) {
  X <- as.matrix(X)
  B <- as.matrix(B)
  if(dim(X)[2] == 1){
    z <- X%*%t(B)%*%t(V)
  } else {z <- X%*%B%*%t(V)}
  z <- exp(z)
  d <- apply(z,1,sum)+1
  P <- z/d
  return(P)
}

#############################################################
# Algorithm for SRRMLR problem
#############################################################

srrmlr <- function(X, Y, r, s, eps = 1e-3, maxit = 100, inner.maxit = 10, V0 = NULL, B0 = NULL) {

  ## Initialization
  n <- dim(X)[1]
  X <- cbind(rep(1, n), X)
  p <- dim(X)[2] - 1
  q <- dim(Y)[2]
  m <- q - 1
  D <- Y[,-1] - matrix(ncol = m, nrow = n, 1/q)
  if(is.null(V0)) V0 <- svd(D)$v[,1:r]
  if(is.null(B0)) B0 <- matrix(ncol = r, nrow = p + 1, 0)
  C0 <- B0 %*% t(V0)
  K <- diag(t(X)%*%X)
  j <- 0
  while(j < maxit){
    k <- 0
    aset0 <- NULL
    u <- exception(X, B0, V0)
    z <- 2*max(u*(1-u))
    G <- t(X)%*%(Y[,-1]-u)%*%V0
    Gamma <- G/(z*K)

    B <- B0
    while(k < inner.maxit){
      u <- exception(X, B, V0)
      z <- 2*max(u*(1-u))
      deta <- 1/2*z*K*diag((B+Gamma)%*%t(B+Gamma))
      aset <- order(deta[-1], decreasing = TRUE)[1:s]             # active set for row selection in B
      aset <- aset + 1

      ## Update B
      aset1 <- c(1,aset)
      B[aset1,] <- B[aset1,] + solve(t(X[,aset1])%*%X[,aset1])%*%t(X[,aset1])%*%(Y[,-1]-u)%*%V0/z
      B[-aset1,] <- 0

      ## Update Gamma
      u <- exception(X, B, V0)
      z <- 2*max(u*(1-u))
      G <- t(X)%*%(Y[,-1]-u)%*%V0
      Gamma[-aset1,] <- G[-aset1,]/(z*K[-aset1])
      Gamma[aset1,] <- 0

      ## If active set is the same then break
      if(setequal(aset, aset0)) break;
      aset0 <- aset
      k <- k + 1
    }

    ## Fixed B, solve V via SVD decomposion ------
    S <- X %*% B
    u <- exception(X, B, V0)
    z <- 2*max(u*(1-u))
    W <- t(Y[,-1] - u) %*% S + z * V0 %*% t(S) %*% S
    W <- svd(W)
    V <- W$u %*% t(W$v)
    C <- B %*% t(V)

    ## Check whether the objective function has converged
    if(abs(f(Y, X, C)- f(Y, X, C0)) < eps) break;
    j <- j + 1
    C0 <- C
    B0 <- B
    V0 <- V
  }
  out <- list(C = C, B = B, V = V)
  class(out) <- 'srrmlr'
  return(out)
}
