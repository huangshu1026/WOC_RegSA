#### code for the simulation section in the regularized weighting method paper ####
source("WeightingModels.R")

#########################################################
######### Bootstrap based on Inflation Rate #############
#########################################################
## Prepare the data 
infdata <- read.csv("inf_filtered.csv", header = T, sep = ",")
infdata <- infdata[ , -1]
M <- 9
X <- as.matrix(infdata[ , 4:12]) # this is for inflation rate data
Y <- infdata[ , 3]
## estimate the bias and covariance matrix
muy.emp <- mean(Y)
sdy.emp <- sd(Y)
mu.emp <- colMeans(X - Y)
sigma.emp <- cov(X)
## parameters in the simulation 
# sample size to be changed 
n.vec.train <- 2^seq(7, 8, by = 1)
n.test <- 1000
draws <- 100
lambda.iter <- seq(0, 10^5, by = 200)
rmse.train <- array(NA, dim = c(length(n.vec.train), draws, 7+6)) # 10 weighting methods 
rmse.test  <- array(NA, dim = c(length(n.vec.train), draws, 7+6)) # 10 weighting methods 
lambda <- array(NA, dim = c(length(n.vec.train), draws, 6))
for(j in 1:length(n.vec.train)){
  n <- n.vec.train[j]
  for(k in 1:draws){
    set.seed(j*k*123)

    # generate the true state
    y.true <- rnorm(100, muy.emp, sdy.emp)
    x.true <- rmvn(100, mu.emp+muy.emp, sigma.emp)
    muy.true <- mean(y.true)
    sdy.true <- sd(y.true)
    mu.true <- colMeans(x.true - y.true)
    sigma.true <- cov(x.true)
    
    # generate a set of judgments 
    Y.train <- rnorm(n, muy.true, sdy.true)
    X.train <- rmvn(n, mu.true+muy.true, sigma.true)
    Y.test <- rnorm(n.test, muy.true, sdy.true)
    X.test <- rmvn(n.test, mu.true+muy.true, sigma.true)
    
    # in-sample prediction and out-of-sample prediction 
    sa.train <- rowMeans(X.train)
    sa.test <- rowMeans(X.test)
    med.train <- apply(X.train, 1, median)
    med.test <- apply(X.test, 1, median)
    weights.optw <- optw(X.train, Y.train)
    wa.optw.train <- X.train%*%weights.optw
    wa.optw.test <- X.test%*%weights.optw
    weights.optwpos <- optw.sum1pos(X.train, Y.train)
    wa.optwpos.train <- X.train%*%weights.optwpos
    wa.optwpos.test <- X.test%*%weights.optwpos
    weights.cwm <- CWM(X.train, Y.train)$weight.CWM
    wa.cwm.train <- X.train%*%weights.cwm
    wa.cwm.test <- X.test%*%weights.cwm
    weights.ssde <- SeqSearch(X.train, Y.train)$weight.SSDE
    wa.ssde.train <- X.train%*%weights.ssde
    wa.ssde.test <- X.test%*%weights.ssde
    # common correlation 
    bias <- colMeans(X.train - Y.train)
    ccorw_exo <- optw_var_1cor_exo(X.train, Y.train)$weight
    ccorw_bme <- optw_var_1cor_BayesMean(X.train, Y.train)$weight
    ccorw_bmap <- optw_var_1cor_BayesMAP(X.train, Y.train)$weight
    ccorw_mean <- optw_var_1cor_mean(X.train, Y.train)$weight
    allweights <- cbind(ccorw_exo, ccorw_bme, ccorw_bmap, ccorw_mean)
    debias.mat <- matrix(rep(bias, each = nrow(X.train)), ncol = ncol(X.train), nrow = nrow(X.train))
    ccorw_error <- colMeans(((X.train - debias.mat) %*% allweights - Y.train)^2)
    ccorw <- allweights[ , which.min(ccorw_error)]
    debias.mat.test <- matrix(rep(bias, each = nrow(X.test)), ncol = ncol(X.test), nrow = nrow(X.test))
    wa.ccor.train <- (X.train - debias.mat) %*% ccorw
    wa.ccor.test <- (X.test - debias.mat.test) %*% ccorw
    
    wa.regequal.lasso.train <- wa.regequal.ridge.train <- wa.regcwm.lasso.train <- wa.regcwm.ridge.train <- wa.regssde.lasso.train <- wa.regssde.ridge.train <- rep(NA, length(sa.train))
    wa.regequal.lasso.test <- wa.regequal.ridge.test <- wa.regcwm.lasso.test <- wa.regcwm.ridge.test <- wa.regssde.lasso.test <- wa.regssde.ridge.test <- rep(NA, length(sa.test))
    lambda.best <- rep(NA, 6)
    lambda.rmse.piece <- matrix(NA, ncol = 6, nrow = length(lambda.iter))
    lambda.rsdse.piece <- matrix(NA, ncol = 6, nrow = length(lambda.iter))
    tryCatch({
      weights.regEqual.lasso <- RegW_EstLambda_h1o_lasso(X.train, Y.train, lambda.iter, wprior = rep(1/ncol(X.train), ncol(X.train)))
      wa.regequal.lasso.train <- X.train%*%weights.regEqual.lasso$reg_weights
      wa.regequal.lasso.test <- X.test%*%weights.regEqual.lasso$reg_weights
      lambda.best[1] <- weights.regEqual.lasso$est_lambda
      lambda.rmse.piece[ , 1] <- weights.regEqual.lasso$rmse 
      lambda.rsdse.piece[ , 1] <- weights.regEqual.lasso$rsdse
    }, error = function(e){cat("No Solution of RegEqual (Lasso)!")})
    tryCatch({
      weights.regEqual.ridge <- RegW_EstLambda_h1o_ridge(X.train, Y.train, lambda.iter, wprior = rep(1/ncol(X.train), ncol(X.train)))
      wa.regequal.ridge.train <- X.train%*%weights.regEqual.ridge$reg_weights
      wa.regequal.ridge.test <- X.test%*%weights.regEqual.ridge$reg_weights
      lambda.best[2] <- weights.regEqual.ridge$est_lambda
      lambda.rmse.piece[ , 2] <- weights.regEqual.ridge$rmse 
      lambda.rsdse.piece[ , 2] <- weights.regEqual.ridge$rsdse
    }, error = function(e){cat("No Solution of RegEqual (Ridge)! ")})
    tryCatch({
      weights.regCWM.lasso <- RegW_EstLambda_h1o_lasso(X.train, Y.train, lambda.iter, wprior = weights.cwm)
      wa.regcwm.lasso.train <- X.train%*%weights.regCWM.lasso$reg_weights
      wa.regcwm.lasso.test <- X.test%*%weights.regCWM.lasso$reg_weights
      lambda.best[3] <- weights.regCWM.lasso$est_lambda
      lambda.rmse.piece[ , 3] <- weights.regCWM.lasso$rmse 
      lambda.rsdse.piece[ , 3] <- weights.regCWM.lasso$rsdse
    }, error = function(e){cat("No Solution of RegCWM (Lasso)! ")})
    tryCatch({
      weights.regCWM.ridge <- RegW_EstLambda_h1o_ridge(X.train, Y.train, lambda.iter, wprior = weights.cwm)
      wa.regcwm.ridge.train <- X.train%*%weights.regCWM.ridge$reg_weights
      wa.regcwm.ridge.test <- X.test%*%weights.regCWM.ridge$reg_weights
      lambda.best[4] <- weights.regCWM.ridge$est_lambda
      lambda.rmse.piece[ , 4] <- weights.regCWM.ridge$rmse 
      lambda.rsdse.piece[ , 4] <- weights.regCWM.ridge$rsdse
    }, error = function(e) cat("No Solution of RegCWM (Ridge)! "))
    tryCatch({
      weights.regSSDE.lasso <- RegW_EstLambda_h1o_lasso(X.train, Y.train, lambda.iter, wprior = weights.ssde)
      wa.regssde.lasso.train <- X.train%*%weights.regSSDE.lasso$reg_weights
      wa.regssde.lasso.test <- X.test%*%weights.regSSDE.lasso$reg_weights
      lambda.best[5] <- weights.regSSDE.lasso$est_lambda
      lambda.rmse.piece[ , 5] <- weights.regSSDE.lasso$rmse 
      lambda.rsdse.piece[ , 5] <- weights.regSSDE.lasso$rsdse
    }, error = function(e) cat("No Solution of RegSSDE (Lasso)! "))
    tryCatch({
      weights.regSSDE.ridge <- RegW_EstLambda_h1o_ridge(X.train, Y.train, lambda.iter, wprior = weights.ssde)
      wa.regssde.ridge.train <- X.train%*%weights.regSSDE.ridge$reg_weights
      wa.regssde.ridge.test <- X.test%*%weights.regSSDE.ridge$reg_weights
      lambda.best[6] <- weights.regSSDE.ridge$est_lambda
      lambda.rmse.piece[ , 6] <- weights.regSSDE.ridge$rmse 
      lambda.rsdse.piece[ , 6] <- weights.regSSDE.ridge$rsdse
    }, error = function(e) cat("No Solution of RegSSDE (Ridge)! "))
    all.pred.train <- cbind(sa.train, med.train, wa.optw.train, wa.optwpos.train, wa.cwm.train, wa.ssde.train, wa.ccor.train,
                            wa.regequal.lasso.train, wa.regequal.ridge.train, wa.regcwm.lasso.train, wa.regcwm.ridge.train, 
                            wa.regssde.lasso.train, wa.regssde.ridge.train)
    all.pred.test <- cbind(sa.test, med.test, wa.optw.test, wa.optwpos.test, wa.cwm.test, wa.ssde.test, wa.ccor.test,
                           wa.regequal.lasso.test, wa.regequal.ridge.test, wa.regcwm.lasso.test, wa.regcwm.ridge.test, 
                           wa.regssde.lasso.test, wa.regssde.ridge.test)
    
    # record all lambda values and prediction errors 
    rmse.train[j, k, ] <- sqrt(colMeans((all.pred.train - Y.train)^2))
    rmse.test[j, k, ]  <- sqrt(colMeans((all.pred.test - Y.test)^2))
    lambda[j, k, ] <- lambda.best
    print(k)
    # recoding data line by line
    temp4 <- cbind(lambda.rmse.piece, rep(n, length(lambda.iter)), rep(k, length(lambda.iter)), lambda.iter)
    write.table(temp4, "inf-lambda_rmse.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
    temp5 <- cbind(lambda.rsdse.piece, rep(n, length(lambda.iter)), rep(k, length(lambda.iter)), lambda.iter)
    write.table(temp5, "inf-lambda_rsdse.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  # write csv file 
  temp1 <- cbind(rmse.train[j, , ], rep(n.vec.train[j], draws), 1:draws)
  write.table(temp1, "inf-rmse_train.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  temp2 <- cbind(rmse.test[j, , ], rep(n.vec.train[j], draws), 1:draws)
  write.table(temp2, "inf-rmse_test.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  temp3 <- cbind(lambda[j, , ], rep(n.vec.train[j], draws), 1:draws)
  write.table(temp3, "inf-lambda.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
}


############################################################
######### Bootstrap based on Unemployment Rate #############
############################################################
## Prepare the data 
unempdata <- read.csv("unemp_filtered.csv", header = T, sep = ",")
unempdata <- unempdata[ , -1]
M <- 8
X <- as.matrix(unempdata[ , 4:11]) # this is for inflation rate data
Y <- unempdata[ , 3]
## estimate the bias and covariance matrix
muy.emp <- mean(Y)
sdy.emp <- sd(Y)
mu.emp <- colMeans(X - Y)
sigma.emp <- cov(X)
## parameters in the simulation 
# sample size to be changed 
n.vec.train <- 2^seq(7, 8, by = 1)
n.test <- 1000
draws <- 100
lambda.iter <- seq(0, 10^5, by = 200)
rmse.train <- array(NA, dim = c(length(n.vec.train), draws, 7+6)) # 10 weighting methods 
rmse.test  <- array(NA, dim = c(length(n.vec.train), draws, 7+6)) # 10 weighting methods 
lambda <- array(NA, dim = c(length(n.vec.train), draws, 6))
for(j in 1:length(n.vec.train)){
  n <- n.vec.train[j]
  for(k in 1:draws){
    set.seed(j*k*123)
    
    # generate the true state
    y.true <- rnorm(100, muy.emp, sdy.emp)
    x.true <- rmvn(100, mu.emp+muy.emp, sigma.emp)
    muy.true <- mean(y.true)
    sdy.true <- sd(y.true)
    mu.true <- colMeans(x.true - y.true)
    sigma.true <- cov(x.true)
    
    # generate a set of judgments 
    Y.train <- rnorm(n, muy.true, sdy.true)
    X.train <- rmvn(n, mu.true+muy.true, sigma.true)
    Y.test <- rnorm(n.test, muy.true, sdy.true)
    X.test <- rmvn(n.test, mu.true+muy.true, sigma.true)
    
    # in-sample prediction and out-of-sample prediction 
    sa.train <- rowMeans(X.train)
    sa.test <- rowMeans(X.test)
    med.train <- apply(X.train, 1, median)
    med.test <- apply(X.test, 1, median)
    weights.optw <- optw(X.train, Y.train)
    wa.optw.train <- X.train%*%weights.optw
    wa.optw.test <- X.test%*%weights.optw
    weights.optwpos <- optw.sum1pos(X.train, Y.train)
    wa.optwpos.train <- X.train%*%weights.optwpos
    wa.optwpos.test <- X.test%*%weights.optwpos
    weights.cwm <- CWM(X.train, Y.train)$weight.CWM
    wa.cwm.train <- X.train%*%weights.cwm
    wa.cwm.test <- X.test%*%weights.cwm
    weights.ssde <- SeqSearch(X.train, Y.train)$weight.SSDE
    wa.ssde.train <- X.train%*%weights.ssde
    wa.ssde.test <- X.test%*%weights.ssde
    # common correlation 
    bias <- colMeans(X.train - Y.train)
    ccorw_exo <- optw_var_1cor_exo(X.train, Y.train)$weight
    ccorw_bme <- optw_var_1cor_BayesMean(X.train, Y.train)$weight
    ccorw_bmap <- optw_var_1cor_BayesMAP(X.train, Y.train)$weight
    ccorw_mean <- optw_var_1cor_mean(X.train, Y.train)$weight
    allweights <- cbind(ccorw_exo, ccorw_bme, ccorw_bmap, ccorw_mean)
    debias.mat <- matrix(rep(bias, each = nrow(X.train)), ncol = ncol(X.train), nrow = nrow(X.train))
    ccorw_error <- colMeans(((X.train - debias.mat) %*% allweights - Y.train)^2)
    ccorw <- allweights[ , which.min(ccorw_error)]
    debias.mat.test <- matrix(rep(bias, each = nrow(X.test)), ncol = ncol(X.test), nrow = nrow(X.test))
    wa.ccor.train <- (X.train - debias.mat) %*% ccorw
    wa.ccor.test <- (X.test - debias.mat.test) %*% ccorw
    
    wa.regequal.lasso.train <- wa.regequal.ridge.train <- wa.regcwm.lasso.train <- wa.regcwm.ridge.train <- wa.regssde.lasso.train <- wa.regssde.ridge.train <- rep(NA, length(sa.train))
    wa.regequal.lasso.test <- wa.regequal.ridge.test <- wa.regcwm.lasso.test <- wa.regcwm.ridge.test <- wa.regssde.lasso.test <- wa.regssde.ridge.test <- rep(NA, length(sa.test))
    lambda.best <- rep(NA, 6)
    lambda.rmse.piece <- matrix(NA, ncol = 6, nrow = length(lambda.iter))
    lambda.rsdse.piece <- matrix(NA, ncol = 6, nrow = length(lambda.iter))
    tryCatch({
      weights.regEqual.lasso <- RegW_EstLambda_h1o_lasso(X.train, Y.train, lambda.iter, wprior = rep(1/ncol(X.train), ncol(X.train)))
      wa.regequal.lasso.train <- X.train%*%weights.regEqual.lasso$reg_weights
      wa.regequal.lasso.test <- X.test%*%weights.regEqual.lasso$reg_weights
      lambda.best[1] <- weights.regEqual.lasso$est_lambda
      lambda.rmse.piece[ , 1] <- weights.regEqual.lasso$rmse 
      lambda.rsdse.piece[ , 1] <- weights.regEqual.lasso$rsdse
    }, error = function(e){cat("No Solution of RegEqual (Lasso)!")})
    tryCatch({
      weights.regEqual.ridge <- RegW_EstLambda_h1o_ridge(X.train, Y.train, lambda.iter, wprior = rep(1/ncol(X.train), ncol(X.train)))
      wa.regequal.ridge.train <- X.train%*%weights.regEqual.ridge$reg_weights
      wa.regequal.ridge.test <- X.test%*%weights.regEqual.ridge$reg_weights
      lambda.best[2] <- weights.regEqual.ridge$est_lambda
      lambda.rmse.piece[ , 2] <- weights.regEqual.ridge$rmse 
      lambda.rsdse.piece[ , 2] <- weights.regEqual.ridge$rsdse
    }, error = function(e){cat("No Solution of RegEqual (Ridge)! ")})
    tryCatch({
      weights.regCWM.lasso <- RegW_EstLambda_h1o_lasso(X.train, Y.train, lambda.iter, wprior = weights.cwm)
      wa.regcwm.lasso.train <- X.train%*%weights.regCWM.lasso$reg_weights
      wa.regcwm.lasso.test <- X.test%*%weights.regCWM.lasso$reg_weights
      lambda.best[3] <- weights.regCWM.lasso$est_lambda
      lambda.rmse.piece[ , 3] <- weights.regCWM.lasso$rmse 
      lambda.rsdse.piece[ , 3] <- weights.regCWM.lasso$rsdse
    }, error = function(e){cat("No Solution of RegCWM (Lasso)! ")})
    tryCatch({
      weights.regCWM.ridge <- RegW_EstLambda_h1o_ridge(X.train, Y.train, lambda.iter, wprior = weights.cwm)
      wa.regcwm.ridge.train <- X.train%*%weights.regCWM.ridge$reg_weights
      wa.regcwm.ridge.test <- X.test%*%weights.regCWM.ridge$reg_weights
      lambda.best[4] <- weights.regCWM.ridge$est_lambda
      lambda.rmse.piece[ , 4] <- weights.regCWM.ridge$rmse 
      lambda.rsdse.piece[ , 4] <- weights.regCWM.ridge$rsdse
    }, error = function(e) cat("No Solution of RegCWM (Ridge)! "))
    tryCatch({
      weights.regSSDE.lasso <- RegW_EstLambda_h1o_lasso(X.train, Y.train, lambda.iter, wprior = weights.ssde)
      wa.regssde.lasso.train <- X.train%*%weights.regSSDE.lasso$reg_weights
      wa.regssde.lasso.test <- X.test%*%weights.regSSDE.lasso$reg_weights
      lambda.best[5] <- weights.regSSDE.lasso$est_lambda
      lambda.rmse.piece[ , 5] <- weights.regSSDE.lasso$rmse 
      lambda.rsdse.piece[ , 5] <- weights.regSSDE.lasso$rsdse
    }, error = function(e) cat("No Solution of RegSSDE (Lasso)! "))
    tryCatch({
      weights.regSSDE.ridge <- RegW_EstLambda_h1o_ridge(X.train, Y.train, lambda.iter, wprior = weights.ssde)
      wa.regssde.ridge.train <- X.train%*%weights.regSSDE.ridge$reg_weights
      wa.regssde.ridge.test <- X.test%*%weights.regSSDE.ridge$reg_weights
      lambda.best[6] <- weights.regSSDE.ridge$est_lambda
      lambda.rmse.piece[ , 6] <- weights.regSSDE.ridge$rmse 
      lambda.rsdse.piece[ , 6] <- weights.regSSDE.ridge$rsdse
    }, error = function(e) cat("No Solution of RegSSDE (Ridge)! "))
    all.pred.train <- cbind(sa.train, med.train, wa.optw.train, wa.optwpos.train, wa.cwm.train, wa.ssde.train, wa.ccor.train,
                            wa.regequal.lasso.train, wa.regequal.ridge.train, wa.regcwm.lasso.train, wa.regcwm.ridge.train, 
                            wa.regssde.lasso.train, wa.regssde.ridge.train)
    all.pred.test <- cbind(sa.test, med.test, wa.optw.test, wa.optwpos.test, wa.cwm.test, wa.ssde.test, wa.ccor.test,
                           wa.regequal.lasso.test, wa.regequal.ridge.test, wa.regcwm.lasso.test, wa.regcwm.ridge.test, 
                           wa.regssde.lasso.test, wa.regssde.ridge.test)
    
    # record all lambda values and prediction errors 
    rmse.train[j, k, ] <- sqrt(colMeans((all.pred.train - Y.train)^2))
    rmse.test[j, k, ]  <- sqrt(colMeans((all.pred.test - Y.test)^2))
    lambda[j, k, ] <- lambda.best
    print(k)
    # recoding data line by line
    temp4 <- cbind(lambda.rmse.piece, rep(n, length(lambda.iter)), rep(k, length(lambda.iter)), lambda.iter)
    write.table(temp4, "unemp-lambda_rmse.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
    temp5 <- cbind(lambda.rsdse.piece, rep(n, length(lambda.iter)), rep(k, length(lambda.iter)), lambda.iter)
    write.table(temp5, "unemp-lambda_rsdse.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  # write csv file 
  temp1 <- cbind(rmse.train[j, , ], rep(n.vec.train[j], draws), 1:draws)
  write.table(temp1, "unemp-rmse_train.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  temp2 <- cbind(rmse.test[j, , ], rep(n.vec.train[j], draws), 1:draws)
  write.table(temp2, "unemp-rmse_test.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  temp3 <- cbind(lambda[j, , ], rep(n.vec.train[j], draws), 1:draws)
  write.table(temp3, "unemp-lambda.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
}

