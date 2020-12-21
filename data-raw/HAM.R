## code to prepare `DATASET` dataset goes here

data.all <- read.csv("/Users/dell1/Desktop/2020_RA/BeSS/Supplementary/hmnist/hmnist_8_8_RGB.csv")

n <- dim(data.all)[1]    # Number of observations
p <- dim(data.all)[2] - 1    # Number of columns of predictors
X <- as.matrix(data.all[,1:p])
names <- colnames(X)
X <- scale(X)    # Centering and Scaling

y_label = as.factor(data.all$label)
q <- length(levels(y_label))    # Number of different classes
m <- q - 1
Y <- array(0, dim = c(n, q))    # Producing a matrix with 0-entries
colnames(Y) <- levels(y_label)
for(lab in levels(y_label)){
  set <- which(y_label==lab)
  idx = as.numeric(lab)+1
  Y[set,idx] <- 1
}    # Real classification results

HAM = list(hamimgs = X, hamlab = Y)

usethis::use_data(HAM, overwrite = TRUE)
