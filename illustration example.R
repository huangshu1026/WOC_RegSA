#### code for the illustration example section in the regularized weighting method paper ####
source("WeightingModels.R")

# 1. sample size: small, moderate and large 
# 2. parameters: 
## 2.1. target value: y ~ N(0, 1)
## 2.2. forecasters: bias, variance, correlation, predictive validity
##### define all the parameters #####
M <- 2
n <- c(10, 50, 100)
n.test <- 100
# parameters for all three: two forecasters + target value 
mu <- matrix(c(0,0,0, -1,1,0, -2,2,0, 1,2,0, 2,1,0), ncol = 3, nrow = 5, byrow = T)
sd <- matrix(c(1,1,1, 0.5,2,1, 2,0.5,1), ncol = 3, nrow = 3, byrow = T)
cor1 <- matrix(c(1,0,0, 0,1,0, 0,0,1), ncol = 3, nrow = 3)
cor2 <- matrix(c(1,0.2,0.2, 0.2,1,0.2, 0.2,0.2,1), ncol = 3, nrow = 3)
cor3 <- matrix(c(1,0.5,0.2, 0.5,1,0.2, 0.2,0.2,1), ncol = 3, nrow = 3)
cor4 <- matrix(c(1,0.8,0.2, 0.8,1,0.2, 0.2,0.2,1), ncol = 3, nrow = 3)
cor5 <- matrix(c(1,0.2,0.5, 0.2,1,0.5, 0.5,0.5,1), ncol = 3, nrow = 3)
cor6 <- matrix(c(1,0.5,0.5, 0.5,1,0.5, 0.5,0.5,1), ncol = 3, nrow = 3)
cor7 <- matrix(c(1,0.8,0.5, 0.8,1,0.5, 0.5,0.5,1), ncol = 3, nrow = 3)
cor8 <- matrix(c(1,0.2,0.8, 0.2,1,0.8, 0.8,0.8,1), ncol = 3, nrow = 3)
cor9 <- matrix(c(1,0.5,0.8, 0.5,1,0.8, 0.8,0.8,1), ncol = 3, nrow = 3)
cor10 <- matrix(c(1,0.8,0.8, 0.8,1,0.8, 0.8,0.8,1), ncol = 3, nrow = 3)
cor.array <- array(NA, dim = c(3, 3, 10))
cor.array[ , , 1] <- cor1
cor.array[ , , 2] <- cor2
cor.array[ , , 3] <- cor3
cor.array[ , , 4] <- cor4
cor.array[ , , 5] <- cor5
cor.array[ , , 6] <- cor6
cor.array[ , , 7] <- cor7
cor.array[ , , 8] <- cor8
cor.array[ , , 9] <- cor9
cor.array[ , , 10] <- cor10
cov.array <- array(NA, dim = c(3, 3, 3, 10))
for(i in 1:3){
  for(j in 1:10){
    cov.array[ , , i, j] <- cor2cov(cor.array[ , , j], sd[i, ])
  }
}
##### generate the judgments, calculate weights and get RMSE #####
lambda.vec <- seq(0, 10^5, by = 100)
allweights <- array(NA, dim = c(M, 3, 3, 5, 3, 10)) # 3 weights, 3 sample sizes, 5 biases, 3 variances and 10 correlations
rmse.test <- array(NA, dim = c(3, 3, 5, 3, 10))
for(i in 1:length(n)){
  for(j in 1:nrow(mu)){
    for(k in 1:nrow(sd)){
      for(r in 1:10){
        print(paste(i, j, k, r, sep = "-"))
        # generate judgments
        Xy.train <- rmvn(n[i], mu[j, ], cov.array[ , , k, r])
        X.train <- Xy.train[ , 1:2]
        y.train <- Xy.train[ , 3]
        Xy.test <- rmvn(n.test, mu[j, ], cov.array[ , , k, r])
        X.test <- Xy.test[ , 1:2]
        y.test <- Xy.test[ , 3]
        # calculate weights and rmse 
        weights.optw <- optw(X.train, y.train)
        weights.eq <- rep(1/M, M)
        weights.regEqual.lasso <- RegW_EstLambda_h1o_lasso(X.train, y.train, lambda.vec, wprior = weights.eq)$reg_weights
        allweights[ , , i, j, k, r] <- cbind(weights.optw, weights.eq, weights.regEqual.lasso)
        # prediction 
        wa.optw.test <- X.test%*%weights.optw
        sa.test <- rowMeans(X.test)
        wa.regequal.lasso.test <- X.test%*%weights.regEqual.lasso
        all.pred.test <- cbind(wa.optw.test, sa.test, wa.regequal.lasso.test)
        # calculate rmse
        rmse.test[ , i, j, k, r] <- sqrt(colMeans((all.pred.test - y.test)^2))
      }
    }
  }
}
# record the weights and rmse into csv file 
output <- c()
for(i in 1:length(n)){
  for(j in 1:nrow(mu)){
    for(k in 1:nrow(sd)){
      for(r in 1:10){
        temp1 <- c(i, j, k, r, as.vector(allweights[ , , i, j, k, r]), as.vector(rmse.test[ , i, j, k, r]))
        output <- rbind(output, temp1)
      }
    }
  }
}
output <- as.data.frame(output)
colnames(output) <- c("ni", "muj", "sdk", "corr", "optw1", "optw2", "eqw1", "eqw2", "regw1", "regw2", "rmse_optw", "rmse_eq", "rmse_reg")
# calculate the true weights and find the best performing weights 
truew <- c()
for(i in 1:length(n)){
  for(j in 1:nrow(mu)){
    for(k in 1:nrow(sd)){
      for(r in 1:10){
        Xy.train <- rmvn(100000, mu[j, ], cov.array[ , , k, r])
        X.train <- Xy.train[ , 1:2]
        y.train <- Xy.train[ , 3]
        truew <- rbind(truew, optw(X.train, y.train))
      }
    }
  }
}
output <- cbind(output, truew)
colnames(output) <- c("ni", "muj", "sdk", "corr", "optw1", "optw2", "eqw1", "eqw2", "regw1", "regw2", "rmse_optw", "rmse_eq", "rmse_reg", "truew1", "truew2")
# analyze the output results 
# output_0 <- output[which(abs(output$optw1 - output$eqw1) <= 0.01), ] # optimal weights are equal weights 
# output_diff <- output[which(abs(output$optw1 - output$eqw1) > 0.01), ]
# output_1 <- output_diff[which(abs(output_diff$regw1 - 0.5) <= 0.01), ] # close to equal weights
# output_11 <- output_1[which(abs(output_1$truew1 - 0.5) <= 0.01), ] # equal weights are true weights 
# output_12 <- output_1[which(abs(output_1$truew1 - 0.5) > 0.01), ] # equal weights are not true weights 
# output_2 <- output_diff[which(abs(output_diff$regw1 - output_diff$optw1) <= 0.01), ] # close to the optimal weights
# output_21 <- output_2[which(abs(output_2$truew1 - output_2$optw1) <= 0.01), ] # optimal weights are true weights 
# output_22 <- output_2[which(abs(output_2$truew1 - output_2$optw1) > 0.01), ] # optimal weights are not true weights 
# output_3 <- output_diff[which(abs(output_diff$regw1 - 0.5) > 0.01 & abs(output_diff$regw1 - output_diff$optw1) > 0.01), ]
index1 <- which(abs(output$regw1-0.5)+abs(output$regw2-0.5) <= 0.01)
index2 <- which(abs(output$regw1-output$optw1)+abs(output$regw2-output$optw2) <= 0.01)
index3 <- c(1:nrow(output))[-unique(c(index1, index2))]
output1 <- output[index1, ]
output2 <- output[index2, ]
output3 <- output[index3, ]
index11 <- which(output1$rmse_reg < output1$rmse_eq & output1$rmse_reg < output1$rmse_optw)
index12 <- which(output1$rmse_reg > output1$rmse_eq & output1$rmse_reg > output1$rmse_optw)
index13 <- c(1:nrow(output1))[-unique(c(index11, index12))]
output11 <- output1[index11, ]
output12 <- output1[index12, ]
output13 <- output1[index13, ]
index21 <- which(output2$rmse_reg < output2$rmse_eq & output2$rmse_reg < output2$rmse_optw)
index22 <- which(output2$rmse_reg > output2$rmse_eq & output2$rmse_reg > output2$rmse_optw)
index23 <- c(1:nrow(output2))[-unique(c(index21, index22))]
output21 <- output2[index21, ]
output22 <- output2[index22, ]
output23 <- output2[index23, ]
index31 <- which(output3$rmse_reg < output3$rmse_eq & output3$rmse_reg < output3$rmse_optw)
index32 <- which(output3$rmse_reg > output3$rmse_eq & output3$rmse_reg > output3$rmse_optw)
index33 <- c(1:nrow(output3))[-unique(c(index31, index32))]
output31 <- output3[index31, ]
output32 <- output3[index32, ]
output33 <- output3[index33, ]

minrmse1 <- apply(output1[ , 11:13], 1, min)
minrmse2 <- apply(output2[ , 11:13], 1, min)
minrmse3 <- apply(output3[ , 11:13], 1, min)
bestw1 <- unlist(sapply(1:nrow(output1), function(x) which(output1[x, 11:13] == minrmse1[x])))
bestw2 <- unlist(sapply(1:nrow(output2), function(x) which(output2[x, 11:13] == minrmse2[x])))
bestw3 <- unlist(sapply(1:nrow(output3), function(x) which(output3[x, 11:13] == minrmse3[x])))
table(bestw1)
table(bestw2)
table(bestw3)


size1 <- which(output1$ni == 1)
size2 <- which(output1$ni == 2)
size3 <- which(output1$ni == 3)
length(which(size1 %in% index11))
length(which(size1 %in% index12))
length(which(size1 %in% index13))
length(which(size2 %in% index11))
length(which(size2 %in% index12))
length(which(size2 %in% index13))
length(which(size3 %in% index11))
length(which(size3 %in% index12))
length(which(size3 %in% index13))

size1 <- which(output2$ni == 1)
size2 <- which(output2$ni == 2)
size3 <- which(output2$ni == 3)
length(which(size1 %in% index21))
length(which(size1 %in% index22))
length(which(size1 %in% index23))
length(which(size2 %in% index21))
length(which(size2 %in% index22))
length(which(size2 %in% index23))
length(which(size3 %in% index21))
length(which(size3 %in% index22))
length(which(size3 %in% index23))

size1 <- which(output3$ni == 1)
size2 <- which(output3$ni == 2)
size3 <- which(output3$ni == 3)
length(which(size1 %in% index31))
length(which(size1 %in% index32))
length(which(size1 %in% index33))
length(which(size2 %in% index31))
length(which(size2 %in% index32))
length(which(size2 %in% index33))
length(which(size3 %in% index31))
length(which(size3 %in% index32))
length(which(size3 %in% index33))

# visualization 
View(output3)
setEPS()
postscript("illustration.eps", width = 8, height = 4)
plot(-10, 1, ylim = c(-2, 2), xlim = c(1, nrow(output3)), xlab = "Scenario Index", ylab = "w1")
abline(h = 0.5, lwd = 3)
for(s in 1:nrow(output3)){
  if(s %in% index31){
    color_spec <- "red"
  }
  if(s %in% index32){
    color_spec <- "blue"
  }
  if(s %in% index33){
    color_spec <- "black"
  }
  points(s, output3$optw1[s], pch = 0, col = color_spec)
  points(s, output3$regw1[s], pch = 1, col = color_spec)
  points(s, output3$truew1[s], pch = 3, col = color_spec)
  #segments(s, output3$optw1[s], s, output3$regw1[s], col = color_spec)
}
legend("bottomleft", c("TrueWeights", "OW", "RegSA"), pch = c(3, 0, 1))
dev.off()
