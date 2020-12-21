#' Multinomial logistic regression model with SRRMLR
#'
#' Construct a regression model with responses be multinomial variable,
#' predictors be continuous high-dimensional data. In addition, this function
#' provides results obtained from other 4 baseline methods.
#'
#' @param X matrix
#' @param Y matrix
#' @param train_prop numeric
#'
#' @return list
#' @import foreach
#' @import doParallel
#' @import MASS
#' @import nnet
#' @import glmnet
#' @import VGAM
#' @importFrom stats runif var
#'
#' @export
#'
#' @examples
#' \dontrun{
#' cv_srrmlr(HAM$hamimgs, HAM$hamlab)
#' }

cv_srrmlr <- function(X, Y, train_prop=0.9){

  ## initialization
  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- dim(Y)[2]
  m <- q - 1
  n.train = floor(train_prop*n)    # Take 90% number of observations as training sets

  train <- sample(1:n, n.train)
  X.train <- X[train,]
  Y.train <- Y[train,]
  X.test <- X[-train,]
  Y.test <- Y[-train,]

  # Take 5 fold of training samples
  size <- ceiling(n.train/5)    # happened to be integer

  r.max <- min(p, m, n)    # maximum rank
  s.max <- 40
  C.hat <- array(NA, dim = c(5, p+1, m))
  Pred <- rep(NA, 5)
  r.hat <- rep(NA, 5)
  s.hat <- rep(NA, 5)
  aset <- vector("list", 5)

  #######################################################
  ## 1. SRRMLR
  err.bess <- matrix(NA, r.max, s.max)
  possibleError1 <- tryCatch({
    for(i in 1:r.max){
      for(j in 1:s.max){
        sum <- rep(NA, 5)
        for(k in 1:5){
          set <- ((k-1)*size+1):min((k*size),n.train)
          fit.bess <- srrmlr(X = X.train[-set,], Y = Y.train[-set,], r = i, s = j)
          loglike <- f(Y = Y.train[set,], X = X.train[set,], C = fit.bess$C[-1,])
          sum[k] <- loglike
        }
        err.bess[i, j] <- mean(sum)    # Average of log-likelihood
      }
    }
  }, error=function(e) e)
  if(inherits(possibleError1, "error")) {} else {
    # The optimal tunning parameters
    id.bess <- which(err.bess == min(err.bess), arr.ind = TRUE)
    bess.r <- id.bess[1]
    bess.s <- id.bess[2]
    bess.fit <- srrmlr(X = X.train, Y = Y.train, r = bess.r, s = bess.s)
    bess.C <- bess.fit$C
    bess.aset <- which(apply(bess.C, 1, sum)[-1] != 0)
    bess.Pred <- f(Y = Y.test, X = X.test, C = bess.C[-1,])

    ##################################
    Pred[1] <- bess.Pred
    r.hat[1] <- bess.r
    s.hat[1] <- bess.s
    aset[[1]] <- bess.aset
    C.hat[1,,] <- bess.C
  }


  ##################################################
  ## 2. GLMNET(ungrouped)
  possibleError2 <- tryCatch({
    unglm.fit <- cv.glmnet(X.train, Y.train, family = "multinomial", type.measure = "deviance", nfolds = 5, type.multinomial = "ungrouped")
    c1 <- coef(unglm.fit, s = "lambda.min")
    c2 <- NULL
    for(j in 2:q){
      z <- as.matrix(c1[[j]])
      colnames(z) <- rownames(z) <- NULL
      c2 <- cbind(c2, z)
    }

    unglm.C <- c2
    unglm.aset <- which(apply(unglm.C, 1, sum)[-1] != 0)
    unglm.s <- length(unglm.aset)
    unglm.Pred <- f(Y = Y.test, X = X.test, C = unglm.C[-1,])}, error=function(e) e)
  if(inherits(possibleError2, "error")) {} else {
    Pred[2] <- unglm.Pred
    s.hat[2] <- unglm.s
    aset[[2]] <- unglm.aset
    C.hat[2,,] <- unglm.C
  }


  ##################################################
  ## 3. GLMNET(grouped)
  possibleError3 <- tryCatch({
    glm.fit <- cv.glmnet(X.train, Y.train, family = "multinomial", type.measure = "deviance", nfolds = 5, type.multinomial = "grouped")
    c1 <- coef(glm.fit, s = "lambda.min")
    c2 <- NULL
    for(j in 2:q){
      z <- as.matrix(c1[[j]])
      colnames(z) <- rownames(z) <- NULL
      c2 <- cbind(c2, z)
    }

    glm.C <- c2
    glm.aset <- which(apply(glm.C, 1, sum)[-1] != 0)
    glm.s <- length(glm.aset)
    glm.Pred <- f(Y = Y.test, X = X.test, C = glm.C[-1,])}, error=function(e) e)
  if(inherits(possibleError3, "error")) {} else {
    Pred[3] <- glm.Pred
    s.hat[3] <- glm.s
    aset[[3]] <- glm.aset
    C.hat[3,,] <- glm.C
  }

  #################################################
  ## 4. multinom
  possibleError4 <- tryCatch({mul.fit <- multinom(Y.train ~ X.train, MaxNWts = 5000, trace = FALSE)
  mul.C <- t(coef(mul.fit))
  mul.aset <- which(apply(mul.C, 1, sum)[-1] != 0)
  mul.s <- length(mul.aset)
  mul.Pred <- f(Y = Y.test, X = X.test, C = mul.C[-1,])}, error=function(e) e)
  if(inherits(possibleError4, "error")) {} else {
    if(mul.Pred == Inf){} else {
      Pred[4] <- mul.Pred
      s.hat[4] <- mul.s
      aset[[4]] <- mul.aset
      C.hat[4,,] <- mul.C
    }
  }


  ##################################################
  ## 5. rrvglm
  err.vglm <- rep(NA, r.max)
  for(i in 1:r.max){
    sum <- rep(NA, 5)
    possibleError5 <- tryCatch({
      for(k in 1:5){
        set <- ((k-1)*size+1):min((k*size),n.train)
        fit.vglm <- rrvglm(Y.train[-set,] ~ X.train[-set,], multinomial(refLevel = 1), Rank = i)
        loglike <- f(Y = Y.train[set,], X = X.train[set,], C = coef(fit.vglm, matrix = TRUE)[-1,])
        sum[k] <- loglike}
    }, error=function(e) e)
    if(inherits(possibleError5, "error")) {} else {err.vglm[i] <- mean(sum)}
  }

  # The optimal tunning parameters
  if(all(is.na(err.vglm))){} else {
    vglm.r <- which.min(err.vglm)
    possibleError6 <- tryCatch({
      vglm.fit <- rrvglm(Y.train ~ X.train, multinomial(refLevel = 1), Rank = vglm.r)
    }, error=function(e) e)
    if(inherits(possibleError6, "error")) {} else {
      vglm.C <- coef(vglm.fit, matrix = TRUE)
      vglm.aset <- which(apply(vglm.C, 1, sum)[-1] != 0)
      vglm.s <- length(vglm.aset)
      vglm.Pred <- f(Y = Y.test, X = X.test, C = vglm.C[-1,])
      if(vglm.Pred == Inf){} else {
        Pred[5] <- vglm.Pred
        r.hat[5] <- vglm.r
        s.hat[5] <- vglm.s
        aset[[5]] <- vglm.aset
        C.hat[5,,] <- vglm.C
      }}
  }
  #############################
  ## OUTPUT
  tab <- cbind(Pred, r.hat, s.hat)
  list(C.hat = C.hat, tab = tab, aset = aset)

}

parfun <- function(iter, X, Y, train_prop){
  sim <- cv_srrmlr(X = X, Y = Y, train_prop = train_prop)
  return(sim$tab)
}
