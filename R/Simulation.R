######################################
# Generate simulation data
######################################
gendata <- function(n, p, q = 8, s = 5, r = 5, type, C = NULL) {
  ## Generate the covariance matrix of predictors
  Sigma <- array(NA, dim = c(p, p))
  if(type == "CS"){
    for(i in 1:p){
      for(j in 1:p){
        if(i == j){
          Sigma[i,j] <- 1
        } else {
          Sigma[i,j] <- 0.5
        }
      }
    }
  } else {
    for(i in 1:p){
      for(j in 1:p){
        Sigma[i,j] <- 0.5^(abs(i-j))
      }
    }
  }

  ## Generate the predictors X based on the covariance matrix
  X <- mvrnorm(n, rep(0,p), Sigma)
  X <- cbind(rep(1, n), X)

  ## Generate the coefficient matrix C
  B <- NULL
  V <- NULL
  if(is.null(C)) {
    B1 <- matrix(ncol = r, nrow = s + 1, runif((s + 1) * r, -2/q, 2/q))
    B2 <- matrix(ncol = r, nrow = p - s, 0)
    B <- rbind(B1, B2)
    V <- matrix(ncol = r, nrow = q, runif(q * r, 0, q))
    C <- B %*% t(V)
  }
  Y <- t(apply(exception1(X,C), 1, rmultinom, n = 1, size = 1))
  list(Y = Y, X = X, C = C, B = B, V = V)
}

#############################################
# The main function for simulation studies
#############################################
foo <- function(n, p, type, n.val = 1000, n.test = 1000, q = 8, s = 5, r = 5, s.max = s + 5, r.max = min(s.max, q-1)) {
  ## Generate data
  data <- gendata(n = n, p = p, q = q, s = s, r = r, type = type)
  Y <- data$Y
  X <- data$X
  C <- data$C
  valdata <- gendata(n = n.val, p = p, q = q, s = s, r = r, type = type, C = C)
  testdata <- gendata(n= n.test, p = p, q = q, s = s, r = r, type = type, C = C)
  m <- q - 1
  P.set <- 1:s
  N.set <- setdiff(1:p, P.set)
  C.hat <- array(NA, dim = c(5, p+1, m))
  Pred <- rep(NA, 5)
  Est <- rep(NA, 5)
  r.hat <- rep(NA, 5)
  s.hat <- rep(NA, 5)
  Sen <- rep(NA, 5)
  Spe <- rep(NA, 5)

  #######################################################
  ## 1. SRRMLR
  # Initialization
  C.bess <- array(0, dim = c(r.max*s.max, p+1, m))
  err.bess <- matrix(0, r.max, s.max)

  # Tuning parameter determined by loglikehood on validation dataset
  possibleError1 <- tryCatch({for(i in 1:r.max){
    V0 <- NULL
    B0 <- NULL
    for(j in 1:s.max){
      fit.bess <- srrmlr(X = X[,-1], Y = Y, r = i, s = j, V0 = V0, B0 = B0)
      C.bess[(i-1)*s.max + j,,] <- fit.bess$C
      err.bess[i, j] <- f(Y = valdata$Y, X = valdata$X, C = fit.bess$C)
      V0 <- fit.bess$V
      B0 <- fit.bess$B
    }
  }}, error=function(e) e)
  if(inherits(possibleError1, "error")) {} else {
    id.bess <- which(err.bess == min(err.bess), arr.ind = TRUE)

    # The optimal tunning parameters
    bess.r <- id.bess[1]
    bess.s <- id.bess[2]
    bess.C <- C.bess[(bess.r-1)*s.max + bess.s,,]
    bess.aset <- which(apply(bess.C, 1, sum)[-1] != 0)
    bess.inaset <- setdiff(1:p, bess.aset)
    bess.Est <- sum((bess.C - C[, -1])^2)/((p+1)*m)
    bess.Pred <- f(Y = testdata$Y, X = testdata$X, C = bess.C)
    bess.Sen <- length(intersect(bess.aset, P.set))/length(P.set)
    bess.Spe <- length(intersect(bess.inaset, N.set))/length(N.set)

    ##################################
    Pred[1] <- bess.Pred
    Est[1] <- bess.Est
    r.hat[1] <- bess.r
    s.hat[1] <- bess.s
    Sen[1] <- bess.Sen
    Spe[1] <- bess.Spe
    C.hat[1,,] <- bess.C
  }

  ##################################################
  ## 2. GLMNET(ungrouped)
  fit.unglm <- glmnet(X[, -1], Y, family = "multinomial", type.multinomial = "ungrouped")
  C.unglm <- array(0, dim = c(length(fit.unglm$lambda), p+1, m))
  err.unglm <- rep(0, length(fit.unglm$lambda))

  # Tuning parameter determined by loglikehood on validation dataset
  for(i in 1:length(fit.unglm$lambda)){
    c1 <- coef(fit.unglm, fit.unglm$lambda[i])
    c2 <- NULL
    for(j in 2:q){
      z <- as.matrix(c1[[j]])
      colnames(z) <- rownames(z) <- NULL
      c2 <- cbind(c2, z)
    }
    C.unglm[i,,] <- c2
    err.unglm[i] <- f(Y = valdata$Y, X = valdata$X, C = c2)
  }

  # The optimal tunning parameters
  id.unglm <- which.min(err.unglm)
  unglm.C <- C.unglm[id.unglm,,]
  unglm.aset <- which(apply(unglm.C, 1, sum)[-1] != 0)
  unglm.inaset <- setdiff(1:p, unglm.aset)
  unglm.s <- length(unglm.aset)
  unglm.Est <- sum((unglm.C-C[, -1])^2)/((p+1)*m)
  unglm.Pred <- f(Y = testdata$Y, X = testdata$X, C = unglm.C)
  unglm.Sen <- length(intersect(unglm.aset, P.set))/length(P.set)
  unglm.Spe <- length(intersect(unglm.inaset, N.set))/length(N.set)

  #########################################
  if(unglm.Pred == Inf){} else {
    Pred[2] <- unglm.Pred
  }
  Est[2] <- unglm.Est
  s.hat[2] <- unglm.s
  Sen[2] <- unglm.Sen
  Spe[2] <- unglm.Spe
  C.hat[2,,] <- unglm.C


  ##################################################
  ## 3. GLMNET(grouped)
  fit.glm <- glmnet(X[, -1], Y, family = "multinomial", type.multinomial = "grouped")
  C.glm <- array(0, dim = c(length(fit.glm$lambda), p+1, m))
  err.glm <- rep(0, length(fit.glm$lambda))

  # Tuning parameter determined by loglikehood on validation dataset
  for(i in 1:length(fit.glm$lambda)){
    c1 <- coef(fit.glm, fit.glm$lambda[i])
    c2 <- NULL
    for(j in 2:q){
      z <- as.matrix(c1[[j]])
      colnames(z) <- rownames(z) <- NULL
      c2 <- cbind(c2, z)
    }
    C.glm[i,,] <- c2
    err.glm[i] <- f(Y = valdata$Y, X = valdata$X, C = c2)
  }

  # The optimal tunning parameters
  id.glm <- which.min(err.glm)
  glm.C <- C.glm[id.glm,,]
  glm.aset <- which(apply(glm.C, 1, sum)[-1] != 0)
  glm.inaset <- setdiff(1:p, glm.aset)
  glm.s <- length(glm.aset)
  glm.Est <- sum((glm.C - C[, -1])^2)/((p+1)*m)
  glm.Pred <- f(Y = testdata$Y, X = testdata$X, C = glm.C)
  glm.Sen <- length(intersect(glm.aset, P.set))/length(P.set)
  glm.Spe <- length(intersect(glm.inaset, N.set))/length(N.set)

  #########################################
  Pred[3] <- glm.Pred
  Est[3] <- glm.Est
  s.hat[3] <- glm.s
  Sen[3] <- glm.Sen
  Spe[3] <- glm.Spe
  C.hat[3,,] <- glm.C


  #################################################
  ## 4. NNET
  possibleError2 <- tryCatch({fit.mul <- multinom(Y ~ X[, -1], MaxNWts = 6000, trace = FALSE)
  mul.C <- t(coef(fit.mul))
  mul.aset <- which(apply(mul.C, 1, sum)[-1] != 0)
  mul.inaset <- setdiff(1:p, mul.aset)
  mul.s <- length(mul.aset)
  mul.Est <- sum((mul.C - C[, -1])^2)/((p+1)*m)
  mul.Pred <- f(Y = testdata$Y, X = testdata$X, C = mul.C)
  mul.Sen <- length(intersect(mul.aset, P.set))/length(P.set)
  mul.Spe <- length(intersect(mul.inaset, N.set))/length(N.set)}, error=function(e) e)
  if(inherits(possibleError2, "error")) {} else {
    if(mul.Pred == Inf){} else {
      Pred[4] <- mul.Pred
      Est[4] <- mul.Est
      s.hat[4] <- mul.s
      Sen[4] <- mul.Sen
      Spe[4] <- mul.Spe
      C.hat[4,,] <- mul.C
    }
  }

  ##################################################
  ## 5. RRVGLM
  C.vglm <- array(0, dim = c(r.max, p+1, m))
  err.vglm <- rep(0, r.max)

  # Tuning parameter determined by loglikehood on validation dataset
  possibleError3 <- tryCatch({for(i in 1:r.max){
    fit.vglm <- rrvglm(Y ~ X[, -1], multinomial(refLevel = 1), Rank = i)
    C.vglm[i,,] <- coef(fit.vglm, matrix = TRUE)
    err.vglm[i] <- f(Y = valdata$Y, X = valdata$X, C = C.vglm[i,,])
  }

    # The optimal tunning parameters
    id.vglm <- which.min(err.vglm)
    vglm.C <- C.vglm[id.vglm,,]
    vglm.r <- id.vglm
    vglm.aset <- which(apply(vglm.C, 1, sum)[-1] != 0)
    vglm.inaset <- setdiff(1:p, vglm.aset)
    vglm.s <- length(vglm.aset)
    vglm.Est <- sum((vglm.C - C[, -1])^2)/((p+1)*m)
    vglm.Pred <- f(Y = testdata$Y, X = testdata$X, C = vglm.C)
    vglm.Sen <- length(intersect(vglm.aset, P.set))/length(P.set)
    vglm.Spe <- length(intersect(vglm.inaset, N.set))/length(N.set)}, error=function(e) e)
  if(inherits(possibleError3, "error")) {} else {
    if(vglm.Pred == Inf){} else {
      Pred[5] <- vglm.Pred
      Est[5] <- vglm.Est
      r.hat[5] <- vglm.r
      s.hat[5] <- vglm.s
      Sen[5] <- vglm.Sen
      Spe[5] <- vglm.Spe
      C.hat[5,,] <- vglm.C
    }
  }

  #############################
  ## OUTPUT
  tab <- cbind(Pred, Est, r.hat, s.hat, Sen, Spe)
  if(n <= p){
    list(C.hat = C.hat[-5,,], tab = tab[-5,], X = X, Y = Y, C = C)
  } else {
    list(C.hat = C.hat, tab = tab, X = X, Y = Y, C = C)
  }
}
