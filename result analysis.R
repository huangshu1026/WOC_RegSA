#### code for the result analysis in the regularized weighting method paper ####
source("WeightingModels.R")
# additional functions for results analysis 
SignTest <- function(ratios_vec){
  total <- length(ratios_vec)
  success <- length(which(ratios_vec > 0))
  failure <- length(which(ratios_vec < 0))
  #signtest <- binom.test(x = success, n = total, p = 0.5, alternative = "two.sided", conf.level = 0.95)
  signtest <- binom.test(x = c(success, failure), p = 0.5, alternative = "two.sided", conf.level = 0.95)
  pvalue <- signtest$p.value
  estimate <- signtest$estimate
  return(list(pvalue = pvalue, estimate = estimate))
}
SignTest_paired <- function(ratios_vec1, ratios_vec2){
  total <- length(ratios_vec1)
  success <- length(which(ratios_vec1 > ratios_vec2))
  failure <- length(which(ratios_vec1 < ratios_vec2))
  #signtest <- binom.test(x = success, n = total, p = 0.5, alternative = "two.sided", conf.level = 0.95)
  signtest <- binom.test(x = c(success, failure), p = 0.5, alternative = "two.sided", conf.level = 0.95)
  pvalue <- signtest$p.value
  estimate <- signtest$estimate
  return(list(pvalue = pvalue, estimate = estimate))
}

##########################################################
##### Results analysis for common correlation method #####
##########################################################
##### SPF FRB #####
rmse.test.sa <- read.csv("results/rmse_test_spfus.csv", header = F)
rmse.test.ccor <- read.csv("results/frb_ccor_rmse_test.csv", header = F)
res_sa <- rmse.test.sa[ , 3]
res_ccor <- rmse.test.ccor[ , 2]
ratio_ccor <- (res_sa - res_ccor)/res_sa
rmse_ratio_ccor <- cbind(rmse.test.ccor$V1, ratio_ccor)
rmse.ratio.spf <- read.csv("results/rmse_ratio_spf.csv", header = T, stringsAsFactors = F)
rmse.ratio.spf.frb <- rmse.ratio.spf[which(rmse.ratio.spf$Data %in% c(1:37)), ]
rmse.ratio.spf.frb.test <- rmse.ratio.spf.frb[which(rmse.ratio.spf.frb$TrainOrTest == "Test"), ]
rmse.ratio.spf.frb.test <- cbind(rmse.ratio.spf.frb.test[ , c(2:3, 5, 7, 10:13)], rmse_ratio_ccor[ , 2])
colnames(rmse.ratio.spf.frb.test) <- c("data", "samplesize", "OW", "CWM", "SSDE", "RegSA", "RegCWM", "RegSSDE", "CCor")
rmse.ratio.spf.frb.test <- rmse.ratio.spf.frb.test[which(rmse.ratio.spf.frb.test$data %in% c(1:5, 10:37)), ]
# data split ratio = 8
rmse.ratio.spf.frb.test.sub8 <- rmse.ratio.spf.frb.test[which(rmse.ratio.spf.frb.test$samplesize == 8), ]
ratio_frb8 <- rmse.ratio.spf.frb.test.sub8[ , c(3:9)]
med_ratio_frb8 <- apply(ratio_frb8, 2, median)
pvalue_frb8 <- apply(ratio_frb8, 2, function(x) SignTest(x)$pvalue)
round(med_ratio_frb8, 4)
round(pvalue_frb8, 4)
# sign test
names(med_ratio_frb8) <- c("OW", "CWM", "SSDE", "RegSA", "RegCWM", "RegSSDE", "CCor")
SignTest_paired(ratio_frb8$RegE, ratio_frb8$CWM)
SignTest_paired(ratio_frb8$RegE, ratio_frb8$SSDE)
SignTest_paired(ratio_frb8$RegE, ratio_frb8$CCor)
# visualization 
setEPS()
postscript("barplot_frb_median.eps", height = 3.5, width = 5)
plot_frb <- barplot(med_ratio_frb8[c(1:3, 7, 4)]*100, ylab = "%RMSE Reduction to SA", ylim = c(-1, 5))
text(plot_frb, c(-0.5, 0.5, 0.5, 0.5, 0.5), paste(round(med_ratio_frb8[c(1:3, 7, 4)]*100, 2), "%", sep = ""), cex = 1)
text(plot_frb, c(NA, med_ratio_frb8[2:3]*100+0.5, NA, NA), c("", "*", "*", ""), cex = 2)
dev.off()

##### SPF ECB #####
rmse.ecb_inf1 <- read.csv("results/rmse_test_inf1.csv", header = F, stringsAsFactors = F)
rmse.ecb_inf2 <- read.csv("results/rmse_test_inf2.csv", header = F, stringsAsFactors = F)
rmse.ecb_emp1 <- read.csv("results/rmse_test_unemp1.csv", header = F, stringsAsFactors = F)
rmse.ecb_emp2 <- read.csv("results/rmse_test_unemp2.csv", header = F, stringsAsFactors = F)
rmse.ccor <- read.csv("results/ecb_ccor_rmse_test.csv", header = F, stringsAsFactors = F)
# given data split ratio = 8
k = 8
rmse.ecb8 <- rbind(c(rmse.ecb_inf1[which(rmse.ecb_inf1$V1 == k), ], rmse.ccor$V3[which(rmse.ccor$V1 == "inf1" & rmse.ccor$V2 == k)]),
                   c(rmse.ecb_inf2[which(rmse.ecb_inf2$V1 == k), ], rmse.ccor$V3[which(rmse.ccor$V1 == "inf2" & rmse.ccor$V2 == k)]),
                   c(rmse.ecb_emp1[which(rmse.ecb_emp1$V1 == k), ], rmse.ccor$V3[which(rmse.ccor$V1 == "unemp1" & rmse.ccor$V2 == k)]),
                   c(rmse.ecb_emp2[which(rmse.ecb_emp2$V1 == k), ], rmse.ccor$V3[which(rmse.ccor$V1 == "unemp2" & rmse.ccor$V2 == k)]))
rmse.ecb8 <- as.data.frame(matrix(unlist(rmse.ecb8), nrow = 4, ncol = 13))
colnames(rmse.ecb8) <- c("splitratio", "SA", "Med", "OW", "OWpos", "CWM", "CEWM", "SSIN", "SSDE", "RegSA", "RegCWM", "RegSSDE", "CCor")
ratio_ecb8 <- (rmse.ecb8[ , 2] - rmse.ecb8[ , c(4, 6, 9:13)])/rmse.ecb8[ , 2]
med_ratio_ecb8 <- apply(ratio_ecb8, 2, median)
pvalue_ecb8 <- apply(ratio_ecb8, 2, function(x) SignTest(x)$pvalue)
round(med_ratio_ecb8, 4)
round(pvalue_ecb8, 4)
# sign test
SignTest_paired(ratio_ecb8$RegSA, ratio_ecb8$CWM)
SignTest_paired(ratio_ecb8$RegSA, ratio_ecb8$SSDE)
SignTest_paired(ratio_ecb8$RegSA, ratio_ecb8$CCor)
# visualization 
setEPS()
postscript("barplot_ecb_median.eps", height = 3.5, width = 5)
plot_ecb <- barplot(med_ratio_ecb8[c(1:3, 7, 4)]*100, ylab = "%RMSE Reduction to SA", ylim = c(-12, 10))
text(plot_ecb, c(-2, 2, 2, -2, 2), paste(round(med_ratio_ecb8[c(1:3, 7, 4)]*100, 2), "%", sep = ""), cex = 1)
#text(plot_ecb, c(NA, med_ratio_frb8[2:3]*100+0.5, NA), c("", "*", "*", ""), cex = 2)
dev.off()

##### Epidemic Forecasting #####
rmse.ccor <- read.csv("results/epi_ccor_rmse_test.csv", header = F, stringsAsFactors = F)
rmse.sa <- read.csv("results/rmse_test_epi.csv", header = F)
rmse.sa <- rmse.sa[-c(397:399), ]
ratio_ccor <- (rmse.sa$V4 - rmse.ccor$V4)/rmse.sa$V4
ratio.flu <- read.csv("results/rmse_ratio_flu.csv", header = T, stringsAsFactors = F)
ratio.flu <- ratio.flu[which(ratio.flu$TrainOrTest == "Test"), ]
ratio.flu <- cbind(ratio.flu[ , c(2:4, 6, 8, 11:14)], ratio_ccor)
# given data split ratio = 8 
ratio.flu8 <- ratio.flu[which(ratio.flu$Data.split.ratio == 8), ]
ratio.flu8 <- ratio.flu8[ , c(4:10)]
med_ratio_flu8 <- apply(ratio.flu8, 2, median)
pvalue_flu8 <- apply(ratio.flu8, 2, function(x) SignTest(x)$pvalue)
round(med_ratio_flu8, 4)
round(pvalue_flu8, 4)
# sign test
names(med_ratio_flu8) <- c("OW", "CWM", "SSDE", "RegSA", "RegCWM", "RegSSDE", "CCor")
SignTest_paired(ratio.flu8$RegE, ratio.flu8$CWM)
SignTest_paired(ratio.flu8$RegE, ratio.flu8$SSDE)
SignTest_paired(ratio.flu8$RegE, ratio.flu8$ratio_ccor)
SignTest_paired(ratio.flu8$RegE, ratio.flu8$OW)
# visualization 
setEPS()
postscript("barplot_fluall_median.eps", height = 3.5, width = 5)
plot_flu <- barplot(med_ratio_flu8[c(1:3, 7, 4)]*100, ylab = "%RMSE Reduction to SA", ylim = c(-17, 3))
text(plot_flu, c(-1.5, -1.5, 1.5, -1.5, -1.5), paste(round(med_ratio_flu8[c(1:3, 7, 4)]*100, 2), "%", sep = ""), cex = 1)
text(plot_flu, c(med_ratio_flu8[c(1:3, 7, 4)]*100-1), c("**", "", "", "***", ""), cex = 2)
dev.off()

# for each target wILI value 
ratio.flu8 <- ratio.flu[which(ratio.flu$Data.split.ratio == 8 & ratio.flu$Target == 4), ]
ratio.flu8 <- ratio.flu8[ , c(4:10)]
med_ratio_flu8 <- apply(ratio.flu8, 2, median)
pvalue_flu8 <- apply(ratio.flu8, 2, function(x) SignTest(x)$pvalue)
round(med_ratio_flu8, 4)
round(pvalue_flu8, 4)
names(med_ratio_flu8) <- c("OW", "CWM", "SSDE", "RegSA", "RegCWM", "RegSSDE")

##### M4 data #####
rmse.ccor <- read.csv("results/m4/m4_ccor_rmse_test.csv", header = F)
colnames(rmse.ccor) <- c("dataset", "target", "CCor")

rmse.demo <- read.csv("results/m4/part1/rmse_test_m4_demo.csv", header = F)
colnames(rmse.demo) <- c("dataset", "target", "SA", "Med", "OW", "OWpos", "CWM", "CEWM", "SSIN", "SSDE", "RegSA", "RegCWM", "RegSSDE")
rmse.demo <- rmse.demo[ , -c(4, 6, 8, 9)]
rmse.demo.monthly <- read.csv("results/m4/part1/rmse_test_m4_demo_monthly_sub.csv", header = F)
colnames(rmse.demo.monthly) <- c("dataset", "target", "SA", "OW", "CWM", "SSDE", "RegSA", "RegCWM", "RegSSDE")
rmse.demo <- rbind(rmse.demo, rmse.demo.monthly)

rmse.finance <- read.csv("results/m4/part1/rmse_test_m4_finance.csv", header = F)
colnames(rmse.finance) <- c("dataset", "target", "SA", "Med", "OW", "OWpos", "CWM", "CEWM", "SSIN", "SSDE", "RegSA", "RegCWM", "RegSSDE")
rmse.finance <- rmse.finance[ , -c(4, 6, 8, 9)]
rmse.finance.new <- read.csv("results/m4/part2/rmse_test_m4_finance_new.csv", header = F)
colnames(rmse.finance.new) <- c("dataset", "target", "SA", "Med", "OW", "OWpos", "CWM", "CEWM", "SSIN", "SSDE", "RegSA", "RegCWM", "RegSSDE")
rmse.finance.new <- rmse.finance.new[ , -c(4, 6, 8, 9)]
rmse.finance.monthly1 <- read.csv("results/m4/part1/rmse_test_m4_finance_monthly_sub.csv", header = F)
colnames(rmse.finance.monthly1) <- c("dataset", "target", "SA", "OW", "CWM", "SSDE", "RegSA", "RegCWM", "RegSSDE")
rmse.finance.monthly2 <- read.csv("results/m4/part2/rmse_test_m4_finance_monthly_sub.csv", header = F)
colnames(rmse.finance.monthly2) <- c("dataset", "target", "SA", "OW", "CWM", "SSDE", "RegSA", "RegCWM", "RegSSDE")
rmse.finance <- rbind(rmse.finance, rmse.finance.new, rmse.finance.monthly1, rmse.finance.monthly2)

rmse.industry <- read.csv("results/m4/part1/rmse_test_m4_industry.csv", header = F)
colnames(rmse.industry) <- c("dataset", "target", "SA", "Med", "OW", "OWpos", "CWM", "CEWM", "SSIN", "SSDE", "RegSA", "RegCWM", "RegSSDE")
rmse.industry <- rmse.industry[ , -c(4, 6, 8, 9)]
rmse.industry.new <- read.csv("results/m4/part2/rmse_test_m4_industry_new.csv", header = F)
colnames(rmse.industry.new) <- c("dataset", "target", "SA", "Med", "OW", "OWpos", "CWM", "CEWM", "SSIN", "SSDE", "RegSA", "RegCWM", "RegSSDE")
rmse.industry.new <- rmse.industry.new[ , -c(4, 6, 8, 9)]
rmse.industry.monthly1 <- read.csv("results/m4/part1/rmse_test_m4_industry_monthly_sub.csv", header = F)
colnames(rmse.industry.monthly1) <- c("dataset", "target", "SA", "OW", "CWM", "SSDE", "RegSA", "RegCWM", "RegSSDE")
rmse.industry.monthly2 <- read.csv("results/m4/part2/rmse_test_m4_industry_monthly_sub.csv", header = F)
colnames(rmse.industry.monthly2) <- c("dataset", "target", "SA", "OW", "CWM", "SSDE", "RegSA", "RegCWM", "RegSSDE")
rmse.industry <- rbind(rmse.industry, rmse.industry.new, rmse.industry.monthly1, rmse.industry.monthly2)

rmse.macro <- read.csv("results/m4/part1/rmse_test_m4_macro.csv", header = F)
colnames(rmse.macro) <- c("dataset", "target", "SA", "Med", "OW", "OWpos", "CWM", "CEWM", "SSIN", "SSDE", "RegSA", "RegCWM", "RegSSDE")
rmse.macro <- rmse.macro[ , -c(4, 6, 8, 9)]
rmse.macro.new <- read.csv("results/m4/part2/rmse_test_m4_macro_new.csv", header = F)
colnames(rmse.macro.new) <- c("dataset", "target", "SA", "Med", "OW", "OWpos", "CWM", "CEWM", "SSIN", "SSDE", "RegSA", "RegCWM", "RegSSDE")
rmse.macro.new <- rmse.macro.new[ , -c(4, 6, 8, 9)]
rmse.macro.monthly1 <- read.csv("results/m4/part1/rmse_test_m4_macro_monthly_sub.csv", header = F)
colnames(rmse.macro.monthly1) <- c("dataset", "target", "SA", "OW", "CWM", "SSDE", "RegSA", "RegCWM", "RegSSDE")
rmse.macro.monthly2 <- read.csv("results/m4/part2/rmse_test_m4_macro_monthly_sub.csv", header = F)
colnames(rmse.macro.monthly2) <- c("dataset", "target", "SA", "OW", "CWM", "SSDE", "RegSA", "RegCWM", "RegSSDE")
rmse.macro <- rbind(rmse.macro, rmse.macro.new, rmse.macro.monthly1, rmse.macro.monthly2)

rmse.micro <- read.csv("results/m4/part1/rmse_test_m4_micro.csv", header = F)
colnames(rmse.micro) <- c("dataset", "target", "SA", "Med", "OW", "OWpos", "CWM", "CEWM", "SSIN", "SSDE", "RegSA", "RegCWM", "RegSSDE")
rmse.micro <- rmse.micro[ , -c(4, 6, 8, 9)]
rmse.micro.new <- read.csv("results/m4/part2/rmse_test_m4_micro_new.csv", header = F)
colnames(rmse.micro.new) <- c("dataset", "target", "SA", "Med", "OW", "OWpos", "CWM", "CEWM", "SSIN", "SSDE", "RegSA", "RegCWM", "RegSSDE")
rmse.micro.new <- rmse.micro.new[ , -c(4, 6, 8, 9)]
rmse.micro.monthly1 <- read.csv("results/m4/part1/rmse_test_m4_micro_monthly_sub.csv", header = F)
colnames(rmse.micro.monthly1) <- c("dataset", "target", "SA", "OW", "CWM", "SSDE", "RegSA", "RegCWM", "RegSSDE")
rmse.micro.monthly2 <- read.csv("results/m4/part2/rmse_test_m4_micro_monthly_sub.csv", header = F)
colnames(rmse.micro.monthly2) <- c("dataset", "target", "SA", "OW", "CWM", "SSDE", "RegSA", "RegCWM", "RegSSDE")
rmse.micro <- rbind(rmse.micro, rmse.micro.new, rmse.micro.monthly1, rmse.micro.monthly2)

rmse.other <- read.csv("results/m4/part1/rmse_test_m4_other.csv", header = F)
colnames(rmse.other) <- c("dataset", "target", "SA", "Med", "OW", "OWpos", "CWM", "CEWM", "SSIN", "SSDE", "RegSA", "RegCWM", "RegSSDE")
rmse.other <- rmse.other[ , -c(4, 6, 8, 9)]
rmse.other.hourly1 <- read.csv("results/m4/part1/rmse_test_m4_other_hourly_sub.csv", header = F)
colnames(rmse.other.hourly1) <- c("dataset", "target", "SA", "OW", "CWM", "SSDE", "RegSA", "RegCWM", "RegSSDE")
rmse.other.hourly2 <- read.csv("results/m4/part2/rmse_test_m4_other_hourly_sub.csv", header = F)
colnames(rmse.other.hourly2) <- c("dataset", "target", "SA", "OW", "CWM", "SSDE", "RegSA", "RegCWM", "RegSSDE")
rmse.other <- rbind(rmse.other, rmse.other.hourly1, rmse.other.hourly2)

rmse.test <- rbind(rmse.demo, rmse.finance, rmse.industry, rmse.macro, rmse.micro, rmse.other)
rmse.test <- merge(rmse.test, rmse.ccor, by = c("dataset", "target"), all.x = TRUE)
rmse.test <- rmse.test[which(rmse.test$dataset %in% c(1:30)), ]
ratio.m4 <- (rmse.test[ , 3] - rmse.test[ , 4:ncol(rmse.test)])/rmse.test[ , 3]
med_ratio_m4 <- apply(ratio.m4, 2, function(x) median(x, na.rm = T))
pvalue_m4 <- apply(ratio.m4, 2, function(x) SignTest(x)$pvalue)
round(med_ratio_m4, 4)
round(pvalue_m4, 4)
# sign test
SignTest_paired(ratio.m4$RegSA, ratio.m4$CWM)
SignTest_paired(ratio.m4$RegSA, ratio.m4$SSDE)
SignTest_paired(ratio.m4$RegSA, ratio.m4$CCor)
# visualization 
setEPS()
postscript("barplot_m4_median.eps", height = 3.5, width = 5)
plot_m4 <- barplot(med_ratio_m4[c(1:3, 7, 4)]*100, ylab = "%RMSE Reduction to SA", ylim = c(-120, 100))
text(plot_m4, c(-10, 10, 10, 10, 10), paste(round(med_ratio_m4[c(1:3, 7, 4)]*100, 2), "%", sep = ""), cex = 1)
text(plot_m4, c(med_ratio_m4[1]*100-10, med_ratio_m4[c(2, 3, 7, 4)]*100+10), c("***", "***", "***", "***", "***"), cex = 2)
dev.off()

# for each domain and frequency data 
data_index <- sort(unique(rmse.test$dataset))
data_index <- c(1, 4, 2, 3, 6, 9, 7, 8, 11, 14, 12, 13, 16, 19, 17, 18, 21, 24, 22, 23, 27, 26, 30, 28, 29)
med_ratio_m4_each <- c()
pvalue_m4_each <- c()
for(i in 1:length(data_index)){
  index <- which(rmse.test$dataset == data_index[i])
  sub_ratio_m4 <- ratio.m4[index, ]
  med_ratio_m4_each <- rbind(med_ratio_m4_each, apply(sub_ratio_m4, 2, function(x) median(x, na.rm = T)))
  pvalue_m4_each <- rbind(pvalue_m4_each, apply(sub_ratio_m4, 2, function(x) SignTest(x)$pvalue))
}
round(med_ratio_m4_each, 4)
round(pvalue_m4_each, 4)
