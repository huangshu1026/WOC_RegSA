#### code for the data application section in the regularized weighting method paper ####
source("WeightingModels.R")

#### SPF FRB ####
item.name <- c(rep("NGDP", 5), rep("PGDP", 4), rep("RCONSUM", 7), rep("RGDP", 7), rep("RNRESIN", 7), rep("RRESINV", 7))
time.name <- c("currentQ", "futureQ1", "futureQ2", "futureQ3", "futureQ4", 
               "futureQ1", "futureQ2", "futureQ3", "futureQ4",
               rep(c("currentQ", "currentY", "futureQ1", "futureQ2", "futureQ3", "futureQ4", "futureY1"), 4))
M.name <- c(rep(4, 5), rep(5, 4), rep(4, 28))
# pre-settings
lambda.vec <- seq(0, 10^5, by = 100)
crit <- 0.8 # CHANGE crit FOR OTHER DATA SPLIT RATIOS
#
rmse.train <- rmse.test <- matrix(NA, nrow = length(item.name), ncol = 12)
best.lambda <- matrix(NA, nrow = length(item.name), ncol = 3)
lambda.rmse <- array(NA, dim = c(length(lambda.vec), 3, length(item.name)))
# for each subset of data 
for(i in c(1:27, 29:length(item.name))){
  filename <- paste("spf_usa/", item.name[i], "_", time.name[i], "_M", M.name[i], ".csv", sep = "")
  data <- read.csv(filename, stringsAsFactors = F)
  data <- data[ , -1]
  training.index <- 1:floor(nrow(data)*crit)
  testing.index <- (floor(nrow(data)*crit)+1):nrow(data)
  training <- data[training.index, ]
  testing <- data[testing.index, ]
  X.train <- as.matrix(training[ , 5:(ncol(training)-1)])
  Y.train <- as.vector(training$actual)
  X.test <- as.matrix(testing[ , 5:(ncol(testing)-1)])
  Y.test <- as.vector(testing$actual)
  
  # estimating all the weights 
  weights.optw <- optw(X.train, Y.train)
  weights.optwpos <- optw.sum1pos(X.train, Y.train)
  weights.cwm <- CWM(X.train, Y.train)$weight.CWM
  weights.cewm <- CWM(X.train, Y.train)$weight.CEWM
  weights.ssin <- SeqSearch(X.train, Y.train)$weight.SSIN
  weights.ssde <- SeqSearch(X.train, Y.train)$weight.SSDE
  bias <- colMeans(X.train - Y.train)
  ccorw_exo <- optw_var_1cor_exo(X.train, Y.train)$weight
  ccorw_bme <- optw_var_1cor_BayesMean(X.train, Y.train)$weight
  ccorw_bmap <- optw_var_1cor_BayesMAP(X.train, Y.train)$weight
  ccorw_mean <- optw_var_1cor_mean(X.train, Y.train)$weight
  allweights <- cbind(ccorw_exo, ccorw_bme, ccorw_bmap, ccorw_mean)
  debias.mat <- matrix(rep(bias, each = nrow(X.train)), ncol = ncol(X.train), nrow = nrow(X.train))
  ccorw_error <- colMeans(((X.train - debias.mat) %*% allweights - Y.train)^2)
  weights.ccorw <- allweights[ , which.min(ccorw_error)]
  weights.regEqual.lasso <- RegW_EstLambda_h1o_lasso(X.train, Y.train, lambda.vec, wprior = rep(1/ncol(X.train), ncol(X.train)))
  weights.regCWM.lasso <- RegW_EstLambda_h1o_lasso(X.train, Y.train, lambda.vec, wprior = weights.cwm)
  weights.regSSDE.lasso <- RegW_EstLambda_h1o_lasso(X.train, Y.train, lambda.vec, wprior = weights.ssde)
  # in-sample prediction 
  sa.train <- rowMeans(X.train)
  med.train <- apply(X.train, 1, median)
  wa.optw.train <- X.train%*%weights.optw
  wa.optwpos.train <- X.train%*%weights.optwpos
  wa.cwm.train <- X.train%*%weights.cwm
  wa.cewm.train <- X.train%*%weights.cewm
  wa.ssin.train <- X.train%*%weights.ssin
  wa.ssde.train <- X.train%*%weights.ssde
  wa.ccor.train <- (X.train - debias.mat) %*% weights.ccorw
  wa.regequal.lasso.train <- X.train%*%weights.regEqual.lasso$reg_weights
  wa.regcwm.lasso.train <- X.train%*%weights.regCWM.lasso$reg_weights
  wa.regssde.lasso.train <- X.train%*%weights.regSSDE.lasso$reg_weights
  all.pred.train <- cbind(sa.train, med.train, wa.optw.train, wa.optwpos.train, 
                          wa.cwm.train, wa.cewm.train, wa.ssin.train, wa.ssde.train, wa.ccor.train,
                          wa.regequal.lasso.train, wa.regcwm.lasso.train, wa.regssde.lasso.train)
  sa.test <- rowMeans(X.test)
  med.test <- apply(X.test, 1, median)
  wa.optw.test <- X.test%*%weights.optw
  wa.optwpos.test <- X.test%*%weights.optwpos
  wa.cwm.test <- X.test%*%weights.cwm
  wa.cewm.test <- X.test%*%weights.cewm
  wa.ssin.test <- X.test%*%weights.ssin
  wa.ssde.test <- X.test%*%weights.ssde
  debias.mat.test <- matrix(rep(bias, each = nrow(X.test)), ncol = ncol(X.test), nrow = nrow(X.test))
  wa.ccor.test <- (X.test - debias.mat.test) %*% weights.ccorw
  wa.regequal.lasso.test <- X.test%*%weights.regEqual.lasso$reg_weights
  wa.regcwm.lasso.test <- X.test%*%weights.regCWM.lasso$reg_weights
  wa.regssde.lasso.test <- X.test%*%weights.regSSDE.lasso$reg_weights
  all.pred.test <- cbind(sa.test, med.test, wa.optw.test, wa.optwpos.test, 
                         wa.cwm.test, wa.cewm.test, wa.ssin.test, wa.ssde.test, wa.ccor.test,
                         wa.regequal.lasso.test, wa.regcwm.lasso.test, wa.regssde.lasso.test)
  # record all lambda values and prediction errors 
  rmse.train[i, ] <- sqrt(colMeans((all.pred.train - Y.train)^2))
  rmse.test[i, ]  <- sqrt(colMeans((all.pred.test - Y.test)^2))
  best.lambda[i, ] <- c(weights.regEqual.lasso$est_lambda[1], weights.regCWM.lasso$est_lambda[1], weights.regSSDE.lasso$est_lambda[1])
  lambda.rmse[ , , i] <- cbind(weights.regEqual.lasso$rmse, weights.regCWM.lasso$rmse, weights.regSSDE.lasso$rmse)
  print(i)
  # write all the results into csv files 
  write.table(t(c(i, rmse.train[i, ])), "rmse_train_spf_usa.csv", row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
  write.table(t(c(i, rmse.test[i, ])), "rmse_test_spf_usa.csv", row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
  write.table(t(c(i, best.lambda[i, ])), "best_lambda_spf_usa.csv", row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
  write.table(t(cbind(rep(i, length(lambda.vec)), lambda.rmse[ , , i])), "lambda_rmse_spf_usa.csv", row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
} 
#### SPF ECB ####
# inflation - 1 year
infdata <- read.csv("~/Documents/WOC_TEST_2021/2empirical_analysis/inf_filtered.csv", header = T, sep = ",")
infdata <- infdata[ , -1]
M <- 9 # inf
infdata1 <- infdata[which(infdata$Year == 1), ]
X <- as.matrix(infdata1[ , 4:12]) # this is for inflation rate data
Y <- infdata1[ , 3]
# inflation - 2 years
infdata <- read.csv("~/Documents/WOC_TEST_2021/2empirical_analysis/inf_filtered.csv", header = T, sep = ",")
infdata <- infdata[ , -1]
M <- 9 # inf
infdata2 <- infdata[which(infdata$Year == 2), ]
X <- as.matrix(infdata2[ , 4:12]) # this is for inflation rate data
Y <- infdata2[ , 3]
# unemployment - 1 year 
infdata <- read.csv("~/Documents/WOC_TEST_2021/2empirical_analysis/unemp_filtered.csv", header = T, sep = ",")
infdata <- infdata[ , -1]
M <- 8 # unemp
infdata1 <- infdata[which(infdata$Year == 1), ]
X <- as.matrix(infdata1[ , 4:11]) # this is for unemployment rate data
Y <- infdata1[ , 3]
# unemployment - 2 years 
infdata <- read.csv("~/Documents/WOC_TEST_2021/2empirical_analysis/unemp_filtered.csv", header = T, sep = ",")
infdata <- infdata[ , -1]
M <- 8 # unemp
infdata2 <- infdata[which(infdata$Year == 2), ]
X <- as.matrix(infdata2[ , 4:11]) # this is for unemployment rate data
Y <- infdata2[ , 3]

# pre-settings
lambda.vec <- c(seq(0, 1, by = 0.01), seq(2, 100, by = 1), seq(200, 10^5, by = 100))
crit <- seq(0.1, 0.9, by = 0.1)
# for each subset of data 
for(i in 1:length(crit)){
  training.index <- 1:floor(nrow(X)*crit[i])
  testing.index <- (floor(nrow(X)*crit[i])+1):nrow(X)
  X.train <- as.matrix(X[training.index, ])
  Y.train <- as.vector(Y[training.index])
  X.test <- as.matrix(X[testing.index, ])
  Y.test <- as.vector(Y[testing.index])
  
  # estimating all the weights 
  weights.optw <- optw(X.train, Y.train)
  weights.optwpos <- optw.sum1pos(X.train, Y.train)
  weights.cwm <- CWM(X.train, Y.train)$weight.CWM
  weights.cewm <- CWM(X.train, Y.train)$weight.CEWM
  weights.ssin <- SeqSearch(X.train, Y.train)$weight.SSIN
  weights.ssde <- SeqSearch(X.train, Y.train)$weight.SSDE
  bias <- colMeans(X.train - Y.train)
  ccorw_exo <- optw_var_1cor_exo(X.train, Y.train)$weight
  ccorw_bme <- optw_var_1cor_BayesMean(X.train, Y.train)$weight
  ccorw_bmap <- optw_var_1cor_BayesMAP(X.train, Y.train)$weight
  ccorw_mean <- optw_var_1cor_mean(X.train, Y.train)$weight
  allweights <- cbind(ccorw_exo, ccorw_bme, ccorw_bmap, ccorw_mean)
  debias.mat <- matrix(rep(bias, each = nrow(X.train)), ncol = ncol(X.train), nrow = nrow(X.train))
  ccorw_error <- colMeans(((X.train - debias.mat) %*% allweights - Y.train)^2)
  weights.ccor <- allweights[ , which.min(ccorw_error)]
  weights.regEqual.lasso <- RegW_EstLambda_h1o_lasso(X.train, Y.train, lambda.vec, wprior = rep(1/ncol(X.train), ncol(X.train)))
  weights.regCWM.lasso <- RegW_EstLambda_h1o_lasso(X.train, Y.train, lambda.vec, wprior = weights.cwm)
  weights.regSSDE.lasso <- RegW_EstLambda_h1o_lasso(X.train, Y.train, lambda.vec, wprior = weights.ssde)
  # in-sample prediction 
  sa.train <- rowMeans(X.train)
  med.train <- apply(X.train, 1, median)
  wa.optw.train <- X.train%*%weights.optw
  wa.optwpos.train <- X.train%*%weights.optwpos
  wa.cwm.train <- X.train%*%weights.cwm
  wa.cewm.train <- X.train%*%weights.cewm
  wa.ssin.train <- X.train%*%weights.ssin
  wa.ssde.train <- X.train%*%weights.ssde
  wa.ccor.train <- (X.train - debias.mat) %*% weights.ccor
  wa.regequal.lasso.train <- X.train%*%weights.regEqual.lasso$reg_weights
  wa.regcwm.lasso.train <- X.train%*%weights.regCWM.lasso$reg_weights
  wa.regssde.lasso.train <- X.train%*%weights.regSSDE.lasso$reg_weights
  all.pred.train <- cbind(sa.train, med.train, wa.optw.train, wa.optwpos.train, 
                          wa.cwm.train, wa.cewm.train, wa.ssin.train, wa.ssde.train, wa.ccor.train,
                          wa.regequal.lasso.train, wa.regcwm.lasso.train, wa.regssde.lasso.train)
  sa.test <- rowMeans(X.test)
  med.test <- apply(X.test, 1, median)
  wa.optw.test <- X.test%*%weights.optw
  wa.optwpos.test <- X.test%*%weights.optwpos
  wa.cwm.test <- X.test%*%weights.cwm
  wa.cewm.test <- X.test%*%weights.cewm
  wa.ssin.test <- X.test%*%weights.ssin
  wa.ssde.test <- X.test%*%weights.ssde
  debias.mat.test <- matrix(rep(bias, each = nrow(X.test)), ncol = ncol(X.test), nrow = nrow(X.test))
  wa.ccor.test <- (X.test - debias.mat.test) %*% weights.ccor
  wa.regequal.lasso.test <- X.test%*%weights.regEqual.lasso$reg_weights
  wa.regcwm.lasso.test <- X.test%*%weights.regCWM.lasso$reg_weights
  wa.regssde.lasso.test <- X.test%*%weights.regSSDE.lasso$reg_weights
  all.pred.test <- cbind(sa.test, med.test, wa.optw.test, wa.optwpos.test, 
                         wa.cwm.test, wa.cewm.test, wa.ssin.test, wa.ssde.test, wa.ccor.test, 
                         wa.regequal.lasso.test, wa.regcwm.lasso.test, wa.regssde.lasso.test)
  # record all lambda values and prediction errors 
  temp1 <- sqrt(colMeans((all.pred.train - Y.train)^2))
  temp2 <- sqrt(colMeans((all.pred.test - Y.test)^2))
  temp3 <- c(weights.regEqual.lasso$est_lambda[1], weights.regCWM.lasso$est_lambda[1], weights.regSSDE.lasso$est_lambda[1])
  temp4 <- cbind(weights.regEqual.lasso$rmse, weights.regCWM.lasso$rmse, weights.regSSDE.lasso$rmse)
  # write all the results into csv files 
  write.table(t(c(i, temp1)), "rmse_train_inf.csv", row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
  write.table(t(c(i, temp2)), "rmse_test_inf.csv", row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
  write.table(t(c(i, temp3)), "best_lambda_inf.csv", row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
  write.table(cbind(rep(i, 3), t(temp4)), "lambda_rmse_inf.csv", row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
  print(i)
} 

#### M4 data ####
item.name <- c(rep("demo", 5), rep("finance", 5), rep("industry", 5), 
               rep("macro", 5), rep("micro", 5), rep("other", 6))
time.name <- c(rep(c("daily", "monthly", "quarterly", "weekly", "yearly"), 5), 
               c("daily", "hourly", "monthly", "quarterly", "weekly", "yearly"))
training.index.end <- c(rep(c(11, 14, 6, 10, 5), 5), c(11, 38, 14, 6, 10, 5))
testing.index.end <- c(rep(c(14, 18, 8, 13, 6), 5), c(14, 48, 18, 8, 13, 6))
# pre-settings
lambda.vec <- seq(0, 10^5, by = 100)
# for each subset of data 
for(i in 1:length(item.name)){
  # prepare the data 
  filename <- paste("2M4comp2018/", item.name[i],"_",time.name[i],".csv", sep = "")
  filetag <- paste("_", item.name[i],"_",time.name[i],".csv", sep = "")
  data <- read_delim(filename, ",", col_names = TRUE)
  data <- data[ , 1:28]
  training.index <- 1:training.index.end[i]
  testing.index <- (training.index.end[i] + 1):testing.index.end[i]
  item <- unique(data$label)
  item.num <- length(unique(data$label))
  series.num <- nrow(data)/item.num # series.num - 2 > 25
  for(j in 1:item.num){
    subdat <- data[which(data$label ==  item[j]), ]
    training <- subdat[training.index, ]
    testing <- subdat[testing.index, ]
    X.train <- as.matrix(training[ , 4:ncol(training)])
    Y.train <- as.vector(training$target)
    X.test <- as.matrix(testing[ , 4:ncol(testing)])
    Y.test <- as.vector(testing$target)
    # estimating all the weights 
    weights.optw <- optw(X.train, Y.train)
    weights.optwpos <- optw.sum1pos(X.train, Y.train)
    weights.cwm <- CWM(X.train, Y.train)$weight.CWM
    weights.cewm <- CWM(X.train, Y.train)$weight.CEWM
    weights.ssin <- SeqSearch(X.train, Y.train)$weight.SSIN
    weights.ssde <- SeqSearch(X.train, Y.train)$weight.SSDE
    bias <- colMeans(X.train - Y.train)
    ccorw_exo <- optw_var_1cor_exo(X.train, Y.train)$weight
    ccorw_bme <- optw_var_1cor_BayesMean(X.train, Y.train)$weight
    ccorw_bmap <- optw_var_1cor_BayesMAP(X.train, Y.train)$weight
    ccorw_mean <- optw_var_1cor_mean(X.train, Y.train)$weight
    allweights <- cbind(ccorw_exo, ccorw_bme, ccorw_bmap, ccorw_mean)
    debias.mat <- matrix(rep(bias, each = nrow(X.train)), ncol = ncol(X.train), nrow = nrow(X.train))
    ccorw_error <- colMeans(((X.train - debias.mat) %*% allweights - Y.train)^2)
    weights.ccor <- allweights[ , which.min(ccorw_error)]
    weights.regEqual.lasso <- RegW_EstLambda_h1o_lasso(X.train, Y.train, lambda.vec, wprior = rep(1/ncol(X.train), ncol(X.train)))
    weights.regCWM.lasso <- RegW_EstLambda_h1o_lasso(X.train, Y.train, lambda.vec, wprior = weights.cwm)
    weights.regSSDE.lasso <- RegW_EstLambda_h1o_lasso(X.train, Y.train, lambda.vec, wprior = weights.ssde)
    # in-sample prediction 
    sa.train <- rowMeans(X.train)
    med.train <- apply(X.train, 1, median)
    wa.optw.train <- X.train%*%weights.optw
    wa.optwpos.train <- X.train%*%weights.optwpos
    wa.cwm.train <- X.train%*%weights.cwm
    wa.cewm.train <- X.train%*%weights.cewm
    wa.ssin.train <- X.train%*%weights.ssin
    wa.ssde.train <- X.train%*%weights.ssde
    wa.ccor.train <- (X.train - debias.mat) %*% weights.ccor
    wa.regequal.lasso.train <- X.train%*%weights.regEqual.lasso$reg_weights
    wa.regcwm.lasso.train <- X.train%*%weights.regCWM.lasso$reg_weights
    wa.regssde.lasso.train <- X.train%*%weights.regSSDE.lasso$reg_weights
    all.pred.train <- cbind(sa.train, med.train, wa.optw.train, wa.optwpos.train, 
                            wa.cwm.train, wa.cewm.train, wa.ssin.train, wa.ssde.train, wa.ccor.train, 
                            wa.regequal.lasso.train, wa.regcwm.lasso.train, wa.regssde.lasso.train)
    sa.test <- rowMeans(X.test)
    med.test <- apply(X.test, 1, median)
    wa.optw.test <- X.test%*%weights.optw
    wa.optwpos.test <- X.test%*%weights.optwpos
    wa.cwm.test <- X.test%*%weights.cwm
    wa.cewm.test <- X.test%*%weights.cewm
    wa.ssin.test <- X.test%*%weights.ssin
    wa.ssde.test <- X.test%*%weights.ssde
    debias.mat.test <- matrix(rep(bias, each = nrow(X.test)), ncol = ncol(X.test), nrow = nrow(X.test))
    wa.ccor.test <- (X.test - debias.mat.test) %*% weights.ccor
    wa.regequal.lasso.test <- X.test%*%weights.regEqual.lasso$reg_weights
    wa.regcwm.lasso.test <- X.test%*%weights.regCWM.lasso$reg_weights
    wa.regssde.lasso.test <- X.test%*%weights.regSSDE.lasso$reg_weights
    all.pred.test <- cbind(sa.test, med.test, wa.optw.test, wa.optwpos.test, 
                           wa.cwm.test, wa.cewm.test, wa.ssin.test, wa.ssde.test, wa.ccor.test, 
                           wa.regequal.lasso.test, wa.regcwm.lasso.test, wa.regssde.lasso.test)
    # record all lambda values and prediction errors 
    temp1 <- sqrt(colMeans((all.pred.train - Y.train)^2))
    temp2 <- sqrt(colMeans((all.pred.test - Y.test)^2))
    temp3 <- c(weights.regEqual.lasso$est_lambda[1], weights.regCWM.lasso$est_lambda[1], weights.regSSDE.lasso$est_lambda[1])
    temp4 <- cbind(weights.regEqual.lasso$rmse, weights.regCWM.lasso$rmse, weights.regSSDE.lasso$rmse)
    # write all the results into csv files 
    write.table(t(c(i, j, temp1)), "rmse_train_m4.csv", row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
    write.table(t(c(i, j, temp2)), "rmse_test_m4.csv", row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
    write.table(t(c(i, j, temp3)), "best_lambda_m4.csv", row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
    write.table(cbind(rep(i, 3), rep(j, 3), t(temp4)), "lambda_rmse_m4.csv", row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
    print(paste(i, "-", j, sep = ""))
  }
} 
#### Epidemic Forecasting ####
item <- c("wILI1", "wILI2", "wILI3", "wILI4", "peak")
location <- c("loc1", "loc2", "loc3", "loc4", "loc5", "loc6", "loc7", "loc8", "loc9", "loc10", "loc11")
# pre-settings
lambda.vec <- seq(0, 10^5, by = 100)
crit <- seq(0.1, 0.9, by = 0.1)
# for each subset of data 
for(i in 1:length(item)){
  for(j in 1:length(location)){
    for(k in 1:length(crit)){
      # prepare huthe data 
      filename <- paste("flu/", item[i], "_", location[j], ".csv", sep = "")
      data <- read.csv(filename, header = T, stringsAsFactors = F)
      training.index <- 1:floor(nrow(data)*crit[k])
      testing.index <- (floor(nrow(data)*crit[k])+1):nrow(data)
      training <- data[training.index, ]
      testing <- data[testing.index, ]
      X.train <- as.matrix(training[ , c(9:23, 25:30)])
      Y.train <- as.vector(training$obs_value)
      X.test <- as.matrix(testing[ , c(9:23, 25:30)])
      Y.test <- as.vector(testing$obs_value)
      # estimating all the weights 
      weights.optw <- optw(X.train, Y.train)
      weights.optwpos <- optw.sum1pos(X.train, Y.train)
      weights.cwm <- CWM(X.train, Y.train)$weight.CWM
      weights.cewm <- CWM(X.train, Y.train)$weight.CEWM
      weights.ssin <- SeqSearch(X.train, Y.train)$weight.SSIN
      weights.ssde <- SeqSearch(X.train, Y.train)$weight.SSDE
      bias <- colMeans(X.train - Y.train)
      ccorw_exo <- optw_var_1cor_exo(X.train, Y.train)$weight
      ccorw_bme <- optw_var_1cor_BayesMean(X.train, Y.train)$weight
      ccorw_bmap <- optw_var_1cor_BayesMAP(X.train, Y.train)$weight
      ccorw_mean <- optw_var_1cor_mean(X.train, Y.train)$weight
      allweights <- cbind(ccorw_exo, ccorw_bme, ccorw_bmap, ccorw_mean)
      debias.mat <- matrix(rep(bias, each = nrow(X.train)), ncol = ncol(X.train), nrow = nrow(X.train))
      ccorw_error <- colMeans(((X.train - debias.mat) %*% allweights - Y.train)^2)
      weights.ccor <- allweights[ , which.min(ccorw_error)]
      weights.regEqual.lasso <- RegW_EstLambda_h1o_lasso(X.train, Y.train, lambda.vec, wprior = rep(1/ncol(X.train), ncol(X.train)))
      weights.regCWM.lasso <- RegW_EstLambda_h1o_lasso(X.train, Y.train, lambda.vec, wprior = weights.cwm)
      weights.regSSDE.lasso <- RegW_EstLambda_h1o_lasso(X.train, Y.train, lambda.vec, wprior = weights.ssde)
      # in-sample prediction 
      sa.train <- rowMeans(X.train)
      med.train <- apply(X.train, 1, median)
      wa.optw.train <- X.train%*%weights.optw
      wa.optwpos.train <- X.train%*%weights.optwpos
      wa.cwm.train <- X.train%*%weights.cwm
      wa.cewm.train <- X.train%*%weights.cewm
      wa.ssin.train <- X.train%*%weights.ssin
      wa.ssde.train <- X.train%*%weights.ssde
      wa.ccor.train <- (X.train - debias.mat) %*% weights.ccor
      wa.regequal.lasso.train <- X.train%*%weights.regEqual.lasso$reg_weights
      wa.regcwm.lasso.train <- X.train%*%weights.regCWM.lasso$reg_weights
      wa.regssde.lasso.train <- X.train%*%weights.regSSDE.lasso$reg_weights
      all.pred.train <- cbind(sa.train, med.train, wa.optw.train, wa.optwpos.train, 
                              wa.cwm.train, wa.cewm.train, wa.ssin.train, wa.ssde.train, wa.ccor.train, 
                              wa.regequal.lasso.train, wa.regcwm.lasso.train, wa.regssde.lasso.train)
      sa.test <- rowMeans(X.test)
      med.test <- apply(X.test, 1, median)
      wa.optw.test <- X.test%*%weights.optw
      wa.optwpos.test <- X.test%*%weights.optwpos
      wa.cwm.test <- X.test%*%weights.cwm
      wa.cewm.test <- X.test%*%weights.cewm
      wa.ssin.test <- X.test%*%weights.ssin
      wa.ssde.test <- X.test%*%weights.ssde
      debias.mat.test <- matrix(rep(bias, each = nrow(X.test)), ncol = ncol(X.test), nrow = nrow(X.test))
      wa.ccor.test <- (X.test - debias.mat.test) %*% weights.ccor
      wa.regequal.lasso.test <- X.test%*%weights.regEqual.lasso$reg_weights
      wa.regcwm.lasso.test <- X.test%*%weights.regCWM.lasso$reg_weights
      wa.regssde.lasso.test <- X.test%*%weights.regSSDE.lasso$reg_weights
      all.pred.test <- cbind(sa.test, med.test, wa.optw.test, wa.optwpos.test, 
                             wa.cwm.test, wa.cewm.test, wa.ssin.test, wa.ssde.test, wa.ccor.test, 
                             wa.regequal.lasso.test, wa.regcwm.lasso.test, wa.regssde.lasso.test)
      # record all lambda values and prediction errors 
      temp1 <- sqrt(colMeans((all.pred.train - Y.train)^2))
      temp2 <- sqrt(colMeans((all.pred.test - Y.test)^2))
      temp3 <- c(weights.regEqual.lasso$est_lambda[1], weights.regCWM.lasso$est_lambda[1], weights.regSSDE.lasso$est_lambda[1])
      temp4 <- cbind(weights.regEqual.lasso$rmse, weights.regCWM.lasso$rmse, weights.regSSDE.lasso$rmse)
      # write all the results into csv files 
      write.table(t(c(i, j, k, temp1)), "rmse_train_epi.csv", row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
      write.table(t(c(i, j, k, temp2)), "rmse_test_epi.csv", row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
      write.table(t(c(i, j, k, temp3)), "best_lambda_epi.csv", row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
      write.table(cbind(rep(i, 3), rep(j, 3), rep(k, 3), t(temp4)), "lambda_rmse_epi.csv", row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
      print(paste(i, "-", j, "-", k, sep = ""))
    }
  }
} 


