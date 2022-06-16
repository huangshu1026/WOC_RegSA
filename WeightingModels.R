####################################################################################
########## This file includes all previous individual weighting models  ############
########## including CWM, SS, TOPn, RP, OW, CLASSIC LASSO, RIDGE REGRESSION ########
########## and our regularized weighting method with leave-one-out CV ##############
####################################################################################

########## Input: (X, y) ###########
########## Output: weights #########

# loading the package and library
if (!require('glmnet')) install.packages('glmnet'); library('glmnet')
if (!require('mgcv')) install.packages('mgcv'); library('mgcv')
if (!require('Matrix')) install.packages('Matrix'); library('Matrix')
if (!require('MBESS')) install.packages('MBESS'); library('MBESS')
if (!require('matrixcalc')) install.packages('matrixcalc'); library('matrixcalc')
if (!require('quadprog')) install.packages('quadprog'); library('quadprog')
if (!require('readr')) install.packages('readr'); library('readr')
if (!require('tseries')) install.packages('tseries'); library('tseries')
if (!require('PACLasso')) install.packages('PACLasso'); library('PACLasso')

#################################################
########## Crowds Selection Methods #############
#################################################
###### Contribution Weighted Model (CWM) ######
score <- function(X, y){
  outcome <- y
  predict <- rowMeans(X)
  s <- sum((outcome - predict)^2)
  return(s)
}
contribution <- function(X, y){
  n <- nrow(X)
  s <- score(X, y)
  sj <- c()
  for(j in 1:ncol(X)){
    sj[j] <- score(X[ , -j], y)
  }
  c <- -(s - sj)/n
  return(c)
}
CWM <- function(X, y){
  M <- ncol(X)
  cont <- contribution(X, y)
  ## CEWM: equal weighting people with positive contribution
  CEWM.select <- which(cont >= 0)
  CEWM.weight <- rep(0, M)
  CEWM.weight[CEWM.select] <- 1/length(CEWM.select)
  ## CWM: weighted average of people with positive contribution
  CWM.select <- which(cont >= 0)
  CWM.weight <- rep(0, M)
  CWM.weight[CWM.select] <- cont[CWM.select]/sum(cont[CWM.select])
  
  return(list(weight.CEWM = CEWM.weight, weight.CWM = CWM.weight))
  # weight is a M-length vector
}
###### Sequential Searching Model ######
SeqSearch <- function(X, y){
  # set-up
  M <- ncol(X)
  MSE.person <- sapply(1:M, function(i) mean((X[ , i] - y)^2))
  rank.person <- rank(MSE.person) # from lowest mse to the highest mse
  ## increasing way 
  SSin.select <- c(which.min(rank.person))
  mse.train.ssin <- c(mean((X[ , SSin.select] - y)^2))
  for(i in 2:M){
    candidate <- c(1:M)[-SSin.select]
    new.mse <- sapply(1:length(candidate), function(j) mean((rowMeans(X[ , c(SSin.select , candidate[j])]) - y)^2))
    SSin.select[i] <- candidate[which.min(new.mse)]
    mse.train.ssin[i] <- min(new.mse)
  }
  SSin.select <- SSin.select[1:which.min(mse.train.ssin)]
  SSin.weight <- rep(0, M)
  SSin.weight[SSin.select] <- 1/length(SSin.select)

  ## decreasing way 
  SSde.select <- c(1:M)
  mse.train.ssde <- c()
  remove <- c()
  for(i in 1:(M-2)){
    new.mse <- sapply(1:length(SSde.select), function(j) mean((rowMeans(X[ , SSde.select[-which(SSde.select == SSde.select[j])]]) - y)^2))
    remove[i] <- SSde.select[which.min(new.mse)]
    mse.train.ssde[i] <- mean((rowMeans(X[ , -remove]) - y)^2)
    SSde.select <- c(1:M)[-remove]
  }
  new.mse1 <- mean((X[ , SSde.select[1]] - y)^2) 
  new.mse2 <- mean((X[ , SSde.select[2]] - y)^2) 
  remove[M-1] <- SSde.select[3 - which.min(c(new.mse1, new.mse2))]
  mse.train.ssde[M-1] <- min(c(new.mse1, new.mse2))
  remove <- remove[1:which.min(mse.train.ssde)]
  SSde.select <- c(1:M)[-remove]
  SSde.weight <- rep(0, M)
  SSde.weight[SSde.select] <- 1/length(SSde.select)
  
  return(list(weight.SSIN = SSin.weight, weight.SSDE = SSde.weight))
}
###### Top N Models ######
TopNModel <- function(X, y, N){
  # set-up
  M <- ncol(X)
  # selecting the top N model
  MSE.train.AllModel <- sapply(1:M, function(i) mean((X[ , i] - y)^2))
  MSEranking <- rank(MSE.train.AllModel, ties.method = "random")
  choose <- which(MSEranking <= N)
  weight <- rep(0, M)
  weight[choose] <- 1/length(choose)
  return(weight)
} 
###### Rank Performance ######
RankPerformance <- function(X, y){
  # set-up
  M <- ncol(X)
  # rank performance model
  MSE.person <- sapply(1:M, function(i) mean((X[ , i] - y)^2))
  rank.person <- rank(MSE.person, ties.method = "random") # from lowest mse to the highest mse
  MSE.RANK <- sapply(2:M, function(i) mean((rowMeans(X[ , which(rank.person <= i)]) - y)^2))
  select.num <- which.min(MSE.RANK)
  RANK.select <- which(rank.person <= select.num + 1)
  weight <- rep(0, M)
  weight[RANK.select] <- 1/length(RANK.select)
  return(weight)
} # only return the selected model's index 

#################################################
####### Theoretical True Optimal Weights#########
#################################################
###### Optimal weights without constraints on weights ######
optw <- function(x, y){
  n <- nrow(x)
  # estimated bias 
  mu.obs <- colMeans(x - y)
  # estimated covariance matrix 
  sigma.obs <- cov(x)
  
  # this function calculates the observed optimal weights 
  # by using estimated bias and covariance matrix 
  mat <- sigma.obs + mu.obs%*%t(mu.obs)
  sigma.xy <- cov(x, y)
  A <- t(rep(1, ncol(x)))
  QPmodel <- solve.QP(Dmat = nearPD(mat, doSym = TRUE)$mat, dvec = sigma.xy, Amat = t(A), bvec = c(1), meq = 1) 
  w <- QPmodel$solution
  return(w)
}
optw.sum1 <- function(x, y){
  M <- ncol(x)
  mat <- t(x) %*% x
  if(!is.positive.definite(mat)){
    mat <- nearPD(mat)$mat
  } # change the optimal weights
  Rinv <- solve(chol(mat))
  C <- cbind(rep(1, M))
  b <- c(1)
  d <- t(y) %*% x
  weights <- rep(NA, M)
  tryCatch({
    qp.model <- solve.QP(Dmat = Rinv, factorized = TRUE, dvec = d, Amat = C, bvec = b, meq = 1)
    weights <- qp.model$solution
  }, error = function(e){
    cat("No solution for OPTW.SUM1!")
  })
  return(weights)
}
optw.sum1pos <- function(x, y){
  M <- ncol(x)
  mat <- t(x) %*% x
  if(!is.positive.definite(mat)){
    mat <- nearPD(mat)$mat
  } # change the optimal weights
  Rinv <- solve(chol(mat))
  C <- cbind(rep(1, M), diag(M))
  b <- c(1, rep(0, M))
  d <- t(y) %*% x
  weights <- rep(NA, M)
  tryCatch({
    qp.model <- solve.QP(Dmat = Rinv, factorized = TRUE, dvec = d, Amat = C, bvec = b, meq = 1)
    weights <- qp.model$solution
  }, error = function(e){
    cat("No solution for OPTW.SUM1POS!")
  })
  return(weights)
}

####################################################
####### Common Correlation Weighting Model #########
####################################################
## Soule et al., (2020): require debias
Calweight_rho_c <- function(rho_c, sd_vec){
  k <- length(sd_vec)
  cor.mat <- matrix(rep(rho_c, k*k), k, k)
  diag(cor.mat) <- 1
  sigma <- cor2cov(cor.mat, sd_vec)
  sigma <- nearPD(sigma)$mat 
  weight <- colSums(solve(sigma))/sum(solve(sigma))
  return(weight)
}
optw_var_1cor_exo <- function(X, y){
  # obtain unbiased judgments 
  n <- nrow(X)
  mu.obs <- colMeans(X - y)
  debias.mat <- matrix(rep(mu.obs, each = n), ncol = ncol(X), nrow = nrow(X))
  X.new <- X - debias.mat
  # model
  k <- ncol(X.new)
  sd_vec <- apply(X.new, 2, sd)
  if(length(which(sd_vec <= 0)) > 0){
    sd_vec[which(sd_vec <= 0)] <- min(sd_vec[which(sd_vec > 0)])*0.001
  }
  rho_c_seq <- seq(1/(1-k)+0.001, 0.999, length.out = 1000)
  Weights <- c()
  for(i in 1:length(rho_c_seq)){
    tempdata <- Calweight_rho_c(rho_c_seq[i], sd_vec)
    Weights <- rbind(Weights, tempdata)
  }
  error_MSE <- colMeans((X.new %*% t(Weights) - y)^2)
  best_rho_c <- rho_c_seq[which.min(error_MSE)]
  weight <- Calweight_rho_c(best_rho_c, sd_vec)
  return(list(rho_c = best_rho_c, weight = weight, bias = mu.obs))
}
optw_var_1cor_BayesMean <- function(X, y){
  # obtain unbiased judgments 
  n <- nrow(X)
  mu.obs <- colMeans(X - y)
  debias.mat <- matrix(rep(mu.obs, each = n), ncol = ncol(X), nrow = nrow(X))
  X.new <- X - debias.mat
  # model
  k <- ncol(X.new)
  delta0 <- k + 1
  rho_c_0 <- 0
  sd_vec <- apply(X.new, 2, sd)
  if(length(which(sd_vec <= 0)) > 0){
    sd_vec[which(sd_vec <= 0)] <- min(sd_vec[which(sd_vec > 0)])*0.001
  }
  cor_mat <- cor(X.new)
  if(length(which(is.na(cor_mat) == TRUE)) > 0){
    k.na <- length(which(is.na(cor_mat) == TRUE))
    rho_bar <- (sum(na.omit(cor_mat)) - k)/(k*(k-1) - k.na)
  } else {
    rho_bar <- (sum(cor_mat)-k)/(k*(k-1))
  }
  rho_c <- ((delta0-k-1)*rho_c_0 + n*rho_bar) / (delta0+n-k-1)
  weight <- Calweight_rho_c(rho_c, sd_vec)
  return(list(rho_c = rho_c, weight = weight, bias = mu.obs))
}
optw_var_1cor_BayesMAP <- function(X, y){
  # obtain unbiased judgments 
  n <- nrow(X)
  mu.obs <- colMeans(X - y)
  debias.mat <- matrix(rep(mu.obs, each = n), ncol = ncol(X), nrow = nrow(X))
  X.new <- X - debias.mat
  # model 
  k <- ncol(X.new)
  delta0 <- k + 1
  rho_c_0 <- 0
  sd_vec <- apply(X.new, 2, sd)
  if(length(which(sd_vec <= 0)) > 0){
    sd_vec[which(sd_vec <= 0)] <- min(sd_vec[which(sd_vec > 0)])*0.001
  }
  cor_mat <- cor(X.new)
  if(length(which(is.na(cor_mat) == TRUE)) > 0){
    k.na <- length(which(is.na(cor_mat) == TRUE))
    rho_bar <- (sum(na.omit(cor_mat)) - k)/(k*(k-1) - k.na)
  } else {
    rho_bar <- (sum(cor_mat)-k)/(k*(k-1))
  }
  rho_c <- ((delta0-k-1)*rho_c_0 + n*rho_bar) / (delta0+n+k+2)
  weight <- Calweight_rho_c(rho_c, sd_vec)
  return(list(rho_c = rho_c, weight = weight, bias = mu.obs))
}
optw_var_1cor_mean <- function(X, y){
  # obtain unbiased judgments 
  n <- nrow(X)
  mu.obs <- colMeans(X - y)
  debias.mat <- matrix(rep(mu.obs, each = n), ncol = ncol(X), nrow = nrow(X))
  X.new <- X - debias.mat
  # model 
  k <- ncol(X.new)
  sd_vec <- apply(X.new, 2, sd)
  if(length(which(sd_vec <= 0)) > 0){
    sd_vec[which(sd_vec <= 0)] <- min(sd_vec[which(sd_vec > 0)])*0.001
  }
  cor_mat <- cor(X.new)
  if(length(which(is.na(cor_mat) == TRUE)) > 0){
    k.na <- length(which(is.na(cor_mat) == TRUE))
    rho_c <- (sum(na.omit(cor_mat)) - k)/(k*(k-1) - k.na)
  } else {
    rho_c <- (sum(cor_mat)-k)/(k*(k-1))
  }
  weight <- Calweight_rho_c(rho_c, sd_vec)
  return(list(rho_c = rho_c, weight = weight, bias = mu.obs))
}


####################################################
####### Classic LASSO and ridge regression #########
####################################################
##### Lasso regression without constraint#####
lasso.uncon <- function(X, y, int, lambda.seq){
  model.cv <- cv.glmnet(X, y, alpha = 1, lambda = lambda.seq, nfolds = 10, intercept = int)
  best.lambda <- model.cv$lambda.min
  model <- glmnet(X, y, alpha = 1, lambda = best.lambda, intercept = int)
  weights <- as.vector(coef(model))
  return(weights)
}

##### Ridge regression without constraint #####
ridge.uncon <- function(X, y, int, lambda.seq){
  model.cv <- cv.glmnet(X, y, alpha = 0, lambda = lambda.seq, nfolds = 10, intercept = int)
  best.lambda <- model.cv$lambda.min
  model <- glmnet(X, y, alpha = 0, lambda = best.lambda, intercept = int)
  weights <- as.vector(coef(model))
  return(weights)
}

####################################################
####### Our regularized weighting method  ##########
####################################################
RegW_GivenOneLambda_lasso <- function(X, y, lambda, wprior){
  M <- ncol(X)
  n <- nrow(X)
  XTX <- t(X) %*% X
  XTY <- t(X) %*% y
  D <- rbind(cbind(XTX, -XTX), cbind(-XTX, XTX))
  if(!is.positive.definite(D)){D <- nearPD(D)$mat}
  sc <- norm(D, "2")
  A <- cbind(c(rep(1, M), rep(-1, M)), diag(2*M))
  b0 <- c(0, rep(0, 2*M)) # different part 
  d <- c(XTY - XTX%*%wprior, -XTY + XTX%*%wprior) - lambda*rep(1, 2*M) # different part 
  result <- solve.QP(D/sc, d/sc, A, b0, meq = 1)$solution
  weights <- result[1:M] - result[(M+1):(2*M)] + wprior # different part
  return(weights)
}
RegW_EstLambda_h1o_lasso <- function(X.train, Y.train, lambda.vec, wprior){
  M <- ncol(X.train)
  rmse <- rsdse <- c()
  for(i in 1:length(lambda.vec)){
    weights <- matrix(unlist(sapply(1:nrow(X.train), function(x) RegW_GivenOneLambda_lasso(X.train[-x, ], Y.train[-x], lambda.vec[i], wprior))), nrow = M, ncol = nrow(X.train))
    h1o_pred <- diag(X.train %*% weights) 
    rmse[i] <- sqrt(mean((h1o_pred - Y.train)^2))
    rsdse[i] <- sqrt(sd((h1o_pred - Y.train)^2))
  }
  best.lambda <- lambda.vec[which(rmse == min(rmse))]
  best.weights <- RegW_GivenOneLambda_lasso(X.train, Y.train, best.lambda, wprior)
  return(list(rmse = rmse, rsdse = rsdse, est_lambda = best.lambda, reg_weights = best.weights))
}
RegW_GivenOneLambda_ridge <- function(X, y, lambda, wprior){
  M <- ncol(X)
  n <- nrow(X)
  XTX <- t(X) %*% X
  XTY <- t(X) %*% y
  mat <- cbind(XTX + lambda*diag(M))
  sc <- norm(mat, "2")
  A <- cbind(rep(1, M), diag(M))
  b0 <- c(0, rep(0, M)) # different part 
  #result <- solve.QP(mat, (XTY-XTX%*%wprior), A, b0, meq = 0)$solution
  result <- solve.QP(mat/sc, (XTY-XTX%*%wprior)/sc, A, b0, meq = 0)$solution
  #mat <- cbind(XTX + lambda*diag(M), rep(1, M))
  #mat <- rbind(mat, t(c(rep(1, M), 0)))
  #result <- solve(mat) %*% c(XTY-XTX%*%wprior, 0) # different part 
  weights <- result[1:M] + wprior # different part 
  return(weights)
}
RegW_EstLambda_h1o_ridge <- function(X.train, Y.train, lambda.vec, wprior){
  M <- ncol(X.train)
  rmse <- rsdse <- c()
  for(i in 1:length(lambda.vec)){
    weights <- matrix(unlist(sapply(1:nrow(X.train), function(x) RegW_GivenOneLambda_ridge(X.train[-x, ], Y.train[-x], lambda.vec[i], wprior))), nrow = M, ncol = nrow(X.train))
    h1o_pred <- diag(X.train %*% weights) 
    rmse[i] <- sqrt(mean((h1o_pred - Y.train)^2))
    rsdse[i] <- sqrt(sd((h1o_pred - Y.train)^2))
  }
  best.lambda <- lambda.vec[which(rmse == min(rmse))]
  best.weights <- RegW_GivenOneLambda_ridge(X.train, Y.train, best.lambda, wprior)
  return(list(rmse = rmse, rsdse = rsdse, est_lambda = best.lambda, reg_weights = best.weights))
}

