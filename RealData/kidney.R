#------------------------------------------------------------------------------#
#   Application of lpsfrail to real data available in 'collett' package CRAN   #
#            Copyright, Oswaldo Gressani. All rights reserved.                 #
#------------------------------------------------------------------------------#


library("collett")   # Load package to access datasets
library("survival")  # Includes coxph routine
source("lpsfrail.R") # Load function that fits a frailty model with LPS method

# Dataset on survival following kidney transplantation
data <- tplant # 434 observations of 7 variables (see ?tplant for details)


# Fit model with coxph function

fitcoxph <- coxph(Surv(time, status) ~ age + diabetes +
                    frailty.gamma(donor), ties = "breslow", data = tplant)

fitcoxphfrail <- 1/(tail(fitcoxph$history$`frailty.gamma(donor)`$history,1)[1])

# Fit model with lpsfrail function

# Organize data in a list of clusters
x <- list()
nclust <- length(unique(tplant$donor))
for(j in 1:nclust){
  dat_temp <- tplant[which(tplant$donor == j), ]
  x[[j]] <- data.frame(Time = dat_temp$time, Censor = NA, delta = dat_temp$status,
             z1 = dat_temp$age, z2 = dat_temp$diabetes, Cluster = dat_temp$donor)
}

tmax <- max(tplant$time) + 1

fitlps <- lpsfrail(x = x, K = 20, tmax = tmax)

outres <- matrix(0, nrow = 3, ncol = 3)
colnames(outres) <- c("coef", "95%CrI_low", "95%CrI_up")
rownames(outres) <- c("age", "diabetes", "frailty.gamma")
outres[1,] <- c(fitlps$betahat[1], fitlps$CI95[1,])
outres[2,] <- c(fitlps$betahat[2], fitlps$CI95[2,])
outres[3,] <- c(fitlps$frailhat, fitlps$CI95[3,])
outres <- round(outres, 4)


# Compare fitcoxph and fitlps
fitcoxph
outres

