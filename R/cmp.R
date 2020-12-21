#' Title
#'
#' @param n.rep numeric
#' @param X matrix
#' @param Y matrix
#' @param train_prop numeric
#'
#' @return list
#' @importFrom stats rmultinom
#' @importFrom parallel makeCluster stopCluster
#' @export
#' @examples
#' \dontrun{
#' cmp(X = HAM$hamimgs, Y = HAM$hamlab)
#' }
cmp <- function(n.rep=100, X, Y, train_prop=0.9){

set.seed(1234)
set.bg <- set.mul <- set.vglm <- NULL
tab.hat <- array(0, dim = c(n.rep, 5, 3))

cl = makeCluster(2)
registerDoParallel(cl)
fit <- foreach(i = 1:n.rep,
               .export=c('parfun'),
               .packages = c("glmnet","nnet","VGAM"))%dopar%{
                 parfun(iter = i, X = X, Y = Y, train_prop=train_prop)
               }
stopImplicitCluster()
stopCluster(cl)

for(i in 1:n.rep){
  tab <- fit[[i]]
  tab[2:4,2] <- 0
  if(any(is.na(tab[1:3,])) == F){set.bg <- c(set.bg, i)}
  if(any(is.na(tab[4,])) == F){set.mul <- c(set.mul, i)}
  if(any(is.na(tab[5,])) == F){set.vglm <- c(set.vglm, i)}
  tab.hat[i,,] <- tab
}
set <- list(set.bg, set.mul, set.vglm)
ave.bg <- apply(tab.hat[set.bg,1:3,], c(2,3), mean)
ave.mul <- apply(tab.hat[set.mul,4,], 2, mean)
ave.vglm <- apply(tab.hat[set.vglm,5,], 2, mean)
est.ave <- rbind(ave.bg, ave.mul, ave.vglm)
est.ave[2:4,2] <- NA
rownames(est.ave) <- NULL

var.bg <- apply(tab.hat[set.bg,1:3,], c(2,3), var)
var.mul <- apply(tab.hat[set.mul,4,], 2, var)
var.vglm <- apply(tab.hat[set.vglm,5,], 2, var)
est.var <- rbind(var.bg, var.mul, var.vglm)
est.var[2:4,3] <- NA
rownames(est.var) <- NULL

colnames(est.ave) <- colnames(est.var) <- c("Pred", "r.hat", "s.hat")
rownames(est.ave) <- rownames(est.var) <- c("SRRMLR", "GLMNET(ungrouped)", "GLMNET(grouped)", "NNET", "RRVGLM")

return(list('est.ave' = est.ave, 'est.var' = est.var))
}
