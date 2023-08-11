#------------------------------------------------------------------------------#
#          Simulating survival times from a frailty model                      #
#         Copyright, Oswaldo Gressani. All rights reserved.                    #
#------------------------------------------------------------------------------#


simfrail <- function(a = 2.4, b = 4, frailpar = 1.7, betas = c(1.5,-0.3,0.2),
                     ngroups = 3, ni = rep(150, ngroups), censrate = 0.15){
  
  data = list()
  p <- length(betas)
  ainv <- 1 / a
  Cr <- data.frame(matrix(nrow = ngroups, ncol = 2))
  colnames(Cr) <- c("Group", "Censoring rate")
  
  for (i in 1:ngroups) {
    # Draw a frailty term shared across a group
    u <- rgamma(n = 1, shape = frailpar, rate = frailpar)
    data[[i]] <- data.frame(matrix(nrow = ni[i], ncol = p + 3))
    covarnames <- paste0(rep("z", p), seq(p))
    colnames(data[[i]]) <- c("Time", "Censor", "delta", covarnames)
    # Draw covariates and survival time within group
    for (j in 1:ni[i]) {
      for (l in 1:p) data[[i]][j, 3 + l] <- rnorm(1)
      expbetaz <- exp(sum(betas * data[[i]][j, 4:(3 + p)]))
      data[[i]][j, 1] <- b * (-log(1 - runif(1)) / (u * expbetaz)) ^ (ainv)
    }
    # Fill censoring and event indicator
    Tij <- data[[i]]$Time
    rootfun <- function(rho){
      val <- (1-censrate) - mean(exp(-rho * Tij))
      return(val)
    }
    
    rho_lb <- 0.1
    rho_ub <- 0.1
    while(rootfun(rho_lb)>=0){
      rho_lb <- rho_lb/10
    }
    while(rootfun(rho_ub)<=0){
      rho_ub <- rho_ub * 10
    }
  
    rhostar <- uniroot(rootfun, lower = rho_lb, upper = rho_ub)$root
    C <- rexp(n = ni[i], rate = rhostar)
    delta <- as.numeric(data[[i]]$Time <= C)
    data[[i]][,2] <- C
    data[[i]][,3] <- delta
    Cr[i,1] <- i
    Cr[i,2] <-  (1/ni[i])*sum(data[[i]][,3]==0) 
  }
  
  outlist <- list(data = data, Cr = Cr)
  return(outlist)
}
  
# Example
simfrail(ngroups = 5, ni = rep(100, 5), censrate = 0.2)





















