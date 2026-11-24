#------------------------------------------------------------------------------#
#         Simulating survival times from a frailty model                       #
#         simfrail2 routine                                                    #
#         Copyright, Oswaldo Gressani. All rights reserved.                    #
#------------------------------------------------------------------------------#


simfrail2 <- function(a = 5, b = 70, frailpar = 1.5, betas = c(log(2),-0.15),
                     ngroups = 5, ni = rep(40, 5), censrate = 0.10){
  
  data = list()
  p <- length(betas)
  ainv <- 1 / a
  Cr <- data.frame(matrix(nrow = ngroups, ncol = 2))
  colnames(Cr) <- c("Group", "Censoring rate")
  n <- sum(ni)
  
  # Mean and variance of the baseline distribution
  basemean <- b * gamma(1+1/a)
  basevar <- b^2 * (gamma(1+2/a) - (gamma(1+1/a))^2) 
  ttrunk <- 120
  
  for (i in 1:ngroups) {
    # Draw a frailty term shared across a group
    u <- rgamma(n = 1, shape = frailpar, rate = frailpar)
    data[[i]] <- data.frame(matrix(nrow = ni[i], ncol = p + 4))
    covarnames <- paste0(rep("z", p), seq(p))
    colnames(data[[i]]) <- c("Time", "Censor", "delta", covarnames, "Cluster")
    # Draw covariates and survival time within group
    for (j in 1:ni[i]) {
      data[[i]][j, 4:(4+p-1)] <- c(rbinom(1,1,prob = 0.5) , rnorm(1))
      expbetaz <- exp(sum(betas * data[[i]][j, 4:(3 + p)]))
      data[[i]][j, 1] <- b * (-log(1 - runif(1)) / (u * expbetaz)) ^ (ainv)
      #To be adjusted!! truncation
      while(data[[i]][j, 1]>ttrunk){
        data[[i]][j, 1] <- b * (-log(1 - runif(1)) / (u * expbetaz)) ^ (ainv)
      }
      data[[i]][j,6] <- i
    }
  }
  
  # Censoring scheme
  Tij <- unlist(lapply(data, "[[", "Time"))
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
  for (i in 1:ngroups) {
      C <- rexp(n = ni[i], rate = rhostar)
      delta <- as.numeric(data[[i]]$Time <= C)
      data[[i]][,2] <- C
      data[[i]][,3] <- delta
      Cr[i,1] <- i
      Cr[i,2] <- (1/ni[i])*sum(data[[i]][,3]==0)
  }
  
  S0 <- function(t){
    fS0 <- exp(-(t/b)^a)
    return(fS0)
  }
  
  h0 <- function(t){
    fh0 <- (a/(b^a)) * t^(a-1)
    return(fh0)
  }
  
  outlist <- list(data = data, Cr = Cr, S0 = S0, h0 = h0, betas = betas,
                  censrate = censrate, frailpar = frailpar, ttrunk = ttrunk,
                  n = n, ngroups = ngroups, ni = ni)
  return(outlist)
}
  
# Example
# x <- simfrail2(censrate = 0.2, ni = rep(1000, 5))





















