#------------------------------------------------------------------------------#
#                    Simulation study Scenario 4                               #
#         Copyright, Oswaldo Gressani. All rights reserved.                    #
#------------------------------------------------------------------------------#

library("EpiLPS")
library("survival")

rm(list = ls())
source("lpsfrail.R")
source("simfrail2.R")

set.seed(1992)

# Scenario specification
S <- 300               # number of replications
ng <- 50               # number of clusters
nindiv <- rep(20, ng)  # number of individuals per cluster
censrate <- 0.10
routinefail <- 0       # counter for datasets for which the routine failed

# For LPS
betafits <- matrix(0, nrow = S, ncol = 2)
fpfits <- c()
S0fits <- matrix(0, nrow = S, ncol = 1000)
CI90inc <- matrix(0, nrow = S, ncol = 3)
CI95inc <- matrix(0, nrow = S, ncol = 3)
colnames(CI90inc) <- c("beta1","beta2","frail")
colnames(CI95inc) <- c("beta1","beta2","frail")

# For coxph
betafitscoxph <- matrix(0, nrow = S, ncol = 2)
fpfitscoxph <- c()
CI90inc_coxph <- matrix(0, nrow = S, ncol = 3)
CI95inc_coxph <- matrix(0, nrow = S, ncol = 3)
colnames(CI90inc_coxph) <- c("beta1","beta2","frail")
colnames(CI95inc_coxph) <- c("beta1","beta2","frail")

tic <- proc.time()

progress <- utils::txtProgressBar(min = 1, max = S, initial = 1,
                                  style = 3, char =">")

for(s in 1:S){
  # Data generation
  data <- simfrail2(ngroups = ng, ni = nindiv, censrate = censrate)
  
  #---------------------- Fit with lpsfrail ----------------------------------
  lpsfit <- try(lpsfrail(x = data$data, K = 15), silent = TRUE)
  while(substr(lpsfit[1],1,5)=="Error" || anyNA(lpsfit$CI90)){
    routinefail <- routinefail + 1
    data <- simfrail2(ngroups = ng, ni = nindiv, censrate = censrate)
    lpsfit <- try(lpsfrail(x = data$data, K = 15), silent = TRUE)
  }

  # Estimation of baseline survival
  tt <- seq(0, lpsfit$tmax, length = 1000)
  dtt <- tt[2] - tt[1]
  Bf <- EpiLPS:::Rcpp_KercubicBspline(tt, lower = 0, upper = lpsfit$tmax, K = lpsfit$K)
  h0fit <- as.numeric(exp(Bf%*%lpsfit$thetahat))
  S0fit <- exp(-cumsum(h0fit)*dtt)
  S0fits[s,] <- S0fit
  
  # Point estimation of betas and frailty param.
  betafits[s,] <- lpsfit$betahat
  fpfits[s] <- lpsfit$frailhat
  
  # Set estimation of betas and frailty param.
  tarcoeffs <- c(data$betas,data$frailpar)
  CI90inc[s, ] <- lpsfit$CI90[,1] <= tarcoeffs  & tarcoeffs <= lpsfit$CI90[,2] 
  CI95inc[s, ] <- lpsfit$CI95[,1] <= tarcoeffs  & tarcoeffs <= lpsfit$CI95[,2] 
  
  #---------------------- Fit with survival package  --------------------------
  datasurv <- do.call(rbind, data$data)
  fitcoxph <- coxph(Surv(Time, delta) ~ z1 + z2 + frailty.gamma(Cluster), 
                    data = datasurv)
  
  # Point estimation of betas and frailty param.
  betafitscoxph[s, ] <- fitcoxph$coefficients
  fpfitscoxph[s] <- 1/(tail(fitcoxph$history$`frailty.gamma(Cluster)`$history,1)[1])
  
  coxCI90low <- fitcoxph$coefficients-qnorm(0.95) * sqrt(diag(fitcoxph$var))
  coxCI90up <- fitcoxph$coefficients+qnorm(0.95) * sqrt(diag(fitcoxph$var))
  coxCI95low <- fitcoxph$coefficients-qnorm(0.975) * sqrt(diag(fitcoxph$var))
  coxCI95up <- fitcoxph$coefficients+qnorm(0.975) * sqrt(diag(fitcoxph$var))
  
  # Set estimation of betas and frailty param.
  CI90inc_coxph[s,1:2] <- coxCI90low <= data$betas  & data$betas <= coxCI90up 
  CI95inc_coxph[s,1:2] <- coxCI95low <= data$betas  & data$betas <= coxCI95up 

  # print(s)
  # print("LPS")
  # print(round(apply(CI90inc,2,sum)*100/s),2)
  # print("coxph")
  # print(round(apply(CI90inc_coxph,2,sum) * 100 /s),2)
  utils::setTxtProgressBar(progress, s)
  
}
close(progress)

#------------ Output results for LPS
outreslps <- matrix(0, nrow = 3, ncol = 7)
colnames(outreslps) <- c("True", "Mean", "Bias", "ESE", "RMSE", "CP90%", "CP95%")
rownames(outreslps) <- c("beta1","beta2","frailty")
outreslps[,1] <- c(data$betas, data$frailpar)
parestim_mean <- colMeans(cbind(betafits, fpfits))
parestim_mean_mat <- matrix(rep(parestim_mean, S), nrow = S, byrow = TRUE)
partrue_mat   <- cbind(matrix(rep(data$betas,S),nrow = S, byrow = TRUE), 
                     rep(data$frailpar,S))
Bias <- colMeans(cbind(betafits,fpfits) - partrue_mat)
ESE <- sqrt((1/(S-1)) * colSums((cbind(betafits,fpfits) - parestim_mean_mat)^2))
RMSE <- sqrt(colMeans((cbind(betafits,fpfits) - partrue_mat)^2))
outreslps[,2] <- parestim_mean
outreslps[,3] <- Bias
outreslps[,4] <- ESE
outreslps[,5] <- RMSE
outreslps[,6] <- colMeans(CI90inc) * 100
outreslps[,7] <- colMeans(CI95inc) * 100

#------------ Output results for coxph
outrescoxph <- matrix(0, nrow = 3, ncol = 7)
colnames(outrescoxph) <-  c("True", "Mean", "Bias", "ESE", "RMSE", "CP90%", "CP95%")
rownames(outrescoxph) <- c("beta1","beta2","frailty")
outrescoxph[,1] <- c(data$betas, data$frailpar)
parestim_mean_cox <- colMeans(cbind(betafitscoxph, fpfitscoxph))
parestim_mean_matcox <- matrix(rep(parestim_mean_cox, S), nrow = S, byrow = TRUE)
Biascoxph <- colMeans(cbind(betafitscoxph,fpfitscoxph) - partrue_mat)
ESEcoxph <- sqrt((1/(S-1)) * colSums((cbind(betafitscoxph,fpfitscoxph) -
                                 parestim_mean_matcox)^2))
RMSEcoxph <- sqrt(colMeans((cbind(betafitscoxph,fpfitscoxph) - partrue_mat)^2))
outrescoxph[,2] <- parestim_mean_cox
outrescoxph[,3] <- Biascoxph
outrescoxph[,4] <- ESEcoxph
outrescoxph[,5] <- RMSEcoxph
outrescoxph[,6] <- colMeans(CI90inc_coxph) * 100
outrescoxph[,7] <- colMeans(CI95inc_coxph) * 100
outrescoxph[3,6:7] <- c(NA, NA)


png(filename = "lpsS0-S4.png", width = 830, height = 600)
par(mfrow = c(1,1))
# Plot estimates of baseline survival
plot(tt, data$S0(tt), type = "l", xlab = "Survival time", ylab = "Baseline survival",
     xlim = c(0,data$ttrunk))
for(s in 1:S){
  lines(tt, S0fits[s,], type = "l", col = "gray65", lty = 1, lwd = 1)
}
lines(tt, data$S0(tt), type = "l", lwd = 2)
lines(tt, apply(S0fits, 2, "median"), col = "red", lty = 2, lwd = 2)
dev.off()


toc <- proc.time() - tic
toc

write.table(round(outreslps,3),"lpsS4.txt")
write.table(round(outrescoxph,3),"coxphS4.txt")


cat("--------- Simulation results for LPS ---------------------\n")
print(round(outreslps,3))
cat("--------- Simulation results for coxph ---------------------\n")
print(round(outrescoxph,3))
cat("----- Proportion of datasets for which routine failed ------\n")
print(paste0(round(routinefail/S * 100, 2),"%"))













