#------------------------------------------------------------------------------#
#                 Fitting a frailty model with LPS                             #
#         Copyright, Oswaldo Gressani. All rights reserved.                    #
#------------------------------------------------------------------------------#


source("simfrail.R")
library("blapsr")
library("numDeriv")
library("microbenchmark")

# Data
x <- simfrail(ngroups = 3, ni = c(100,150,110), censrate = 0.2)
D <- x$data

tr <- max(unlist(lapply(D, '[[', 1))) # Largest observed follow-up time
K <- 20                               # Number of B-splines
I <- length(D)                        # Number of groups
ni <- lengths(lapply(D, '[[', 1))     # Number of subjects in each group
p <- 3                              
delta <- lapply(D, '[[', 3)           # Indicator function
B <- lapply(lapply(lapply(D, '[[', 1), cubicbs, lower = 0, upper = tr, K = K), '[[', 8)
Z <- list(as.matrix(D[[1]][, 4:6]), as.matrix(D[[2]][, 4:6]), as.matrix(D[[3]][, 4:6]))

#-- Bspline basis 
J <- 300                                                           # Number of bins
phij <- seq(0, tr, length = J + 1)                                 # Grid points
Deltaj <- phij[2] - phij[1]                                        # Interval width
sj <- phij[1:J] + Deltaj * 0.5                                     # Interval midpoint
Bsj <- blapsr::cubicbs(sj, lower = 0, upper = tr, K = K)$Bmatrix   # B-spline basis
upto <- list()
for(j in 1:I){
  upto_temp <-as.integer(D[[j]]$Time/Deltaj) + 1
  upto_temp[which(upto_temp == (J + 1))] <- J
  upto[[j]] <- upto_temp
}
  
# Penalty 
Dr <- diag(K)
penorder <- 2
for(k in 1:penorder){Dr <- diff(Dr)}
P <- t(Dr) %*% Dr          # Penalty matrix of dimension c(K,K)
P <- P + diag(1e-06, K)    # Perturbation to ensure P is full rank

#Log-likelihood
loglik <- function(latent, fp){
  theta <- latent[1:K]
  beta <- tail(latent, p)
  sumd <- sapply(delta, "sum")
  cfr <- fp * log(fp) + lgamma(sumd + fp) - lgamma(fp)
  ll1 <- sapply(Map('+',Map('*',lapply(B,'%*%',theta), delta), 
             Map('*',lapply(Z, '%*%', beta), delta)), "sum")
  log_g1 <- c()
  for(j in 1:I){
  log_g1[j] <- log(sum((cumsum(exp(Bsj%*%theta) * Deltaj)[upto[[j]]]) * exp(Z[[j]]%*%beta)) + fp)
  }
  ll2 <- (sumd + fp) * log_g1
  ll <- sum(cfr + ll1 - ll2)
  return(ll)
}

#Gradient of log-likelihood
Dloglik <- function(latent, fp){
  theta <- latent[1:K]
  beta <- tail(latent, p)
  sumd <- sapply(delta, "sum")
  g1 <- c()
  for(j in 1:I){
    g1[j] <- sum((cumsum(exp(Bsj%*%theta) * Deltaj)[upto[[j]]]) * exp(Z[[j]]%*%beta)) + fp
  }
  # Gradient with respect to theta
  Dtheta1 <- sapply(Map('*',B,delta),"colSums")
  D2mat <- matrix(0, nrow = K, ncol = I)
  for (j in 1:I) {
    for (k in 1:K) {
      D2mat[k, j] <- (sum((cumsum(exp(Bsj %*% theta) * Bsj[, k] * Deltaj)[upto[[j]]]) * 
                            exp(Z[[j]] %*% beta)))
    }
  }
    
  Dtheta2 <- t((sumd + fp) * (1 / g1) * t(D2mat))
  grad_theta <- rowSums(Dtheta1 - Dtheta2)
  # Gradient with respect to beta
  Dbeta1 <- sapply(Map('*',Z,delta),"colSums")
  D3mat <- matrix(0, nrow = p, ncol = I)
  for (j in 1:I) {
    for (m in 1:p) {
      D3mat[m, j] <- (sum((cumsum(exp(Bsj %*% theta) * Deltaj)[upto[[j]]]) * 
                            exp(Z[[j]] %*% beta) * Z[[j]][,m]))
    }
  }
  Dbeta2 <- t((sumd + fp) * (1 / g1) * t(D3mat))
  grad_beta <- as.numeric(rowSums(Dbeta1-Dbeta2))
  gg <- c(grad_theta, grad_beta)
  return(gg)
}

#Hessian of log-likelihood
D2loglik <- function(latent, fp){
  theta <- latent[1:K]
  beta <- tail(latent, p)
  sumd <- sapply(delta, "sum")
  
  g1 <- c()
  for(j in 1:I){
    g1[j] <- sum((cumsum(exp(Bsj%*%theta) * Deltaj)[upto[[j]]]) * exp(Z[[j]]%*%beta)) + fp
  }
  # Gradient with respect to theta
  g2 <- matrix(0, nrow = K, ncol = I)
  for (j in 1:I) {
    for (k in 1:K) {
      g2[k, j] <- (sum((cumsum(exp(Bsj %*% theta) * Bsj[, k] * Deltaj)[upto[[j]]]) * 
                            exp(Z[[j]] %*% beta)))
    }
  }
  
  dg1 <- g2
  dg2 <- list()
  for (j in 1:I) {
    dg2[[j]] <- matrix(0, nrow = K, ncol = K)
    for (k in 1:K) {
      for(l in 1:K) {
        dg2[[j]][k, l] <- (sum((cumsum(exp(Bsj %*% theta) * Bsj[, k] * Bsj[,l] * Deltaj)[upto[[j]]]) * 
                         exp(Z[[j]] %*% beta)))
    }
   }
  }
  
  outermats <- list()
  for(j in 1:I){
    outermats[[j]] <- outer(g2[,j],dg1[,j])
  }
  #Upper-left block
  D2theta <- (-1) * Reduce('+',Map('*',as.list((sumd + fp) * (1 / g1^2)),
                 Map('-',Map('*', dg2, as.list(g1)), outermats)))


  #Upper-right block
  g3 <- matrix(0, nrow = p, ncol = I)
  for (j in 1:I) {
    for (m in 1:p) {
      g3[m, j] <- (sum((cumsum(exp(Bsj %*% theta) * Deltaj)[upto[[j]]]) * 
                            exp(Z[[j]] %*% beta) * Z[[j]][,m]))
    }
  }
  
  dg2b <- list()
  for (j in 1:I) {
    dg2b[[j]] <- matrix(0, nrow = K, ncol = p)
    for (k in 1:K) {
      for(m in 1:p)
      dg2b[[j]][k, m] <- (sum((cumsum(exp(Bsj %*% theta) * Bsj[, k] * Deltaj)[upto[[j]]]) * 
                         exp(Z[[j]] %*% beta) * Z[[j]][,m]))
    }
  }
  
  outermats2 <- list()
  for(j in 1:I){
    outermats2[[j]] <- outer(g2[,j],g3[,j])
  }
  
  Dcross12 <- (-1) * Reduce('+',Map('*',as.list((sumd + fp) * (1 / g1^2)),
                                    Map('-',Map('*', dg2b, as.list(g1)), outermats2)))
  
  Dcross21 <- t(Dcross12)
  
  #Lower-right block
  dg1b <- g3
  
  dg3b <- list()
  for (j in 1:I) {
    dg3b[[j]] <- matrix(0, nrow = p, ncol = p)
    for (l in 1:p) {
      for(m in 1:p)
        dg3b[[j]][l, m] <- (sum((cumsum(exp(Bsj %*% theta) * Deltaj)[upto[[j]]]) * 
                                  exp(Z[[j]] %*% beta) * Z[[j]][,l] * Z[[j]][,m]))
    }
  }
  
  outermats3 <- list()
  for(j in 1:I){
    outermats3[[j]] <- outer(g3[,j],g3[,j])
  }
  
  D2beta <- (-1) * Reduce('+',Map('*',as.list((sumd + fp) * (1 / g1^2)),
                                    Map('-',Map('*', dg3b, as.list(g1)), outermats3)))
  
  Hessmat <- matrix(0, nrow = (K+p), ncol = (K+p))
  Hessmat[1:K,1:K] <- D2theta
  Hessmat[1:K, (K+1):(K+p)] <- Dcross12
  Hessmat[(K+1):(K+p), 1:K] <- Dcross21
  Hessmat[(K+1):(K+p),(K+1):(K+p)] <- D2beta
  
  return(Hessmat)
  
}

# Checks------------------------------------------------------------
set.seed(123)

#Check: analytical gradient vs numerical gradient with numDeriv
latent <- MASS::mvrnorm(n=1, mu = rep(0, (K+p)), Sigma = diag(K+p))
grad(func = loglik, x = latent, fp = 5) # Numerical derivative
Dloglik(latent, fp = 5) # Analytical derivative

round(grad(func = loglik, x = latent, fp = 5)-Dloglik(latent, fp = 5),5)

#Check: analytical Hessian vs numerical Hessian with numDeriv
latent <- MASS::mvrnorm(n=1, mu = rep(0, (K+p)), Sigma = diag(K+p))
hessian(func = loglik, x = latent, fp = 5) # Numerical Hessian
D2loglik(latent, fp = 5) # Analytical Hessian
round(hessian(func = loglik, x = latent, fp = 5)-
        D2loglik(latent, fp = 5),5)


microbenchmark(Dloglik(latent, fp = 5))

microbenchmark(D2loglik(latent, fp = 5))








