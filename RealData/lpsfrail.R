#------------------------------------------------------------------------------#
#               Fitting a Gamma frailty model with LPS                         #
#         Copyright, Oswaldo Gressani. All rights reserved.                    #
#------------------------------------------------------------------------------#

lpsfrail <- function(x, K = 15, tmax = 120){
  
  D <- x
  tr <- tmax                                # Largest observed follow-up time
  tvec <- unlist(lapply(D, '[[', 1))        # Vector of follow-up times
  deltalist <- lapply(D, '[[', 3)           # Event status as list
  sumdeltac <- sapply(deltalist, "sum")     # Sum of indicator functions per cluster
  deltavec <- unlist(deltalist)             # Event status as vector
  BB <- EpiLPS:::Rcpp_KercubicBspline(tvec, lower = 0, upper = tr, K = K)
  Zmat <- as.matrix(do.call(rbind, lapply(D,"[",4:5)))
  
  I <- length(x)                        # Number of groups
  ni <- lengths(lapply(x, '[[', 1))     # Number of subjects in each group
  n <- sum(ni)                          # Total number of subjects
  p <- 2                                # Number of covariates                              
  delta <- lapply(D, '[[', 3)           # Indicator function
  sumdeltac <- sapply(delta, "sum")     # Sum of indicator functions per cluster
  
  gindicmat <- matrix(0, nrow = n, ncol = I)
  rowstart <- 1
  for (i in 1:I) {
    gindicmat[rowstart:(rowstart + ni[i] - 1),i] <- rep(1, ni[i])
    rowstart <- rowstart + ni[i]
  }
    
  #-- Bspline basis 
  J <- 200                                                  # Number of bins
  phij <- seq(0, tr, length = J + 1)                        # Grid points
  Deltaj <- phij[2] - phij[1]                               # Interval width
  sj <- phij[1:J] + Deltaj * 0.5                            # Interval midpoint
  Bsj <- EpiLPS:::Rcpp_KercubicBspline(sj, lower = 0, upper = tr, K = K)  # B-spline basis
  upto <- list()
  for(j in 1:I){
    upto_temp <-as.integer(D[[j]]$Time/Deltaj) + 1
    upto_temp[which(upto_temp == (J + 1))] <- J
    upto[[j]] <- upto_temp
  }
  uptovec <- unlist(upto)
  
  #Hyperparameters for priors
  nu <- 3
  a_kappa <- 1e-04
  b_kappa <- 1e-04
  gamma_prec <- 1e-06
  
  # Penalty 
  Dr <- diag(K)
  penorder <- 2
  for(k in 1:penorder){Dr <- diff(Dr)}
  P <- t(Dr) %*% Dr          # Penalty matrix of dimension c(K,K)
  P <- P + diag(1e-06,K)     # Diagonal perturbation to make P full rank
  Q <- function(lambda){
    Qmat <- matrix(0, nrow = (K+p+1), ncol = (K+p+1))
    Qmat[1:K,1:K] <- lambda * P
    Qmat[(K+1):(K+p+1),(K+1):(K+p+1)] <- gamma_prec * diag(p+1)
    return(Qmat)
  }

  logtar <- function(latent, lambda, grad = FALSE, hess = FALSE){
    # Define quantities
    theta <- latent[1:K]
    beta <- latent[(K+1):(K+p)]
    efp <- exp(tail(latent, 1))
    Zb <- as.numeric(Zmat%*%beta)
    eBsj <- exp(Bsj%*%theta)
    cmsum <- cumsum(eBsj * Deltaj)[uptovec]
    cms <- cmsum * exp(Zb)
    Qlamb <- Q(lambda)
    
    #log-likelihood
    ll1 <- sum(deltavec * ((BB%*%theta) + Zb))
    ll2 <- sum((sumdeltac + efp) * log(colSums(gindicmat * cms) + efp))
    loglik <- I * (efp * log(efp) - lgamma(efp)) + sum(lgamma(sumdeltac + efp)) +
      ll1 - ll2
    
    #log-target (log conditional posterior of latent vector)
    logtarget <-  loglik - 0.5 * sum((latent*Qlamb)%*%latent)
    outlist <- list(loglik = loglik, logtarget = logtarget)
    
    if(isTRUE(grad)){
      cmskmat <- apply((Bsj * as.numeric(eBsj) * Deltaj), 2, "cumsum")[uptovec, ] * exp(Zb)
      g1 <- colSums(gindicmat * cms) + efp
      g1inv <- 1 / g1
      fact <- colSums(t(gindicmat) * (sumdeltac + efp) * g1inv)
      grad_theta <- colSums(BB * deltavec - cmskmat * fact)
      
      # Gradient with respect to beta (only for balanced design)
      grad_beta <- as.numeric(colSums(Zmat * (deltavec - cms * fact)))
      
      # Gradient with respect to gamma_tilde (log of frailty parameter)
      grad_logfp <- efp * sum((1+log(efp)) + digamma(sumdeltac+efp) - digamma(efp) -
                                log(g1) - (sumdeltac+efp) * g1inv)
      
      Dloglik <- c(grad_theta, grad_beta, grad_logfp)
      Dlogplat <- Dloglik - colSums(Qlamb*latent)
      outlist <- list(loglik = loglik, logtarget = logtarget, Dlogtarget = Dlogplat)
    }
    
    if(isTRUE(hess)){
      g1 <- colSums(gindicmat * cms) + efp
      f1 <- (sumdeltac + efp)/(g1^2)
      d2theta <- matrix(0, nrow = K, ncol = K)
      H13 <- c()
      H12 <- matrix(0, nrow = K, ncol = p)  # Block 12 (only for balanced design)
      g2mat <- apply((Bsj * as.numeric(eBsj) * Deltaj), 2, "cumsum")[uptovec, ] * exp(Zb)
      dg1mat <- g2mat
      
      colidx <- 1
      for(k in 1:K){
        for(l in colidx:K){
          g2 <-  colSums(g2mat[,k] * gindicmat)
          dg1 <- colSums(dg1mat[,l] * gindicmat)
          dg2 <- colSums(cumsum(as.numeric(eBsj) * Bsj[,k] * Bsj[,l] * Deltaj)[uptovec] * exp(Zb) * gindicmat)
          d2theta[k,l] <- (-1) * sum(f1 * (dg2 * g1 - g2 * dg1))
        }
        colidx <- colidx + 1
        H13[k] <- efp * sum((f1 - (1/g1)) * g2)
        for(l in 1:p){
          g3l <- colSums(Zmat[,l] * cms * gindicmat)
          dg2_beta <- colSums(g2mat[,k] * Zmat[,l] * gindicmat)
          H12[k,l] <- (-1) * sum(f1 * (dg2_beta * g1 - g2 * g3l))
        }
      }
      H11 <- d2theta + (t(d2theta) - t(d2theta) * diag(1,K))
      
      # Block 22 (only for balanced design)
      H22 <- matrix(0, nrow = p, ncol = p)
      H23 <- c()
      for(l in 1:p){
        g3l <- colSums(Zmat[,l] * cms * gindicmat)
        for(m in 1:p){
          dg3_beta <- colSums(cms * Zmat[,m] * Zmat[,l] * gindicmat)
          g3m <- colSums(Zmat[,m] * cms * gindicmat)
          H22[l,m] <- (-1) * sum(f1 * (dg3_beta * g1 - g3l * g3m))
        }
        H23[l] <- efp * sum((f1-(1/g1))*g3l)
      }
      
      g4 <- (1+log(efp)) + digamma(sumdeltac + efp) - digamma(efp) -
        log(g1) - (sumdeltac + efp) * (1/g1)
      H33 <- efp * (sum(g4) + sum(1+efp*(trigamma(sumdeltac + efp)-trigamma(efp)-2*(1/g1) + f1)))
      
      D2loglik <- matrix(0, nrow = (K + p + 1), ncol = (K+ p + 1))
      D2loglik[1:K,1:K] <- H11
      D2loglik[1:K, (K+1):(K+p)] <- H12
      D2loglik[1:K, (K+p+1)] <- H13
      D2loglik[(K+1):(K+p), 1:K] <- t(H12)
      D2loglik[(K+1):(K+p),(K+1):(K+p)] <- H22
      D2loglik[(K+1):(K+p), (K+p+1)] <- H23
      D2loglik[(K+p+1),] <- c(H13, H23, H33)
      D2logplat <- D2loglik - Qlamb
      if(isTRUE(grad)){
        outlist <- list(loglik = loglik, logtarget = logtarget, Dlogtarget = Dlogplat,
                        D2logtarget = D2logplat)}
      else{
        outlist <- list(loglik = loglik, logtarget = logtarget, D2logtarget = D2logplat)
      }
    }
    return(outlist)
  }
  
  LM <- function(x, lambda=3, maxiter = 100, tol = 1e-04, damp = 1e-03){
    dist <- 2 * tol
    iter <- 0
    dampitermax <- 5
    conv <- 0
    while(dist > tol & iter < maxiter){
      logtarcall <- logtar(x, lambda = lambda, grad = TRUE, hess = TRUE)
      f <- logtarcall$logtarget
      df <- logtarcall$Dlogtarget
      H <- logtarcall$D2logtarget
      
      #Eigen decomposition
      Hdiaginv <- diag((1 / diag(H)))
      Omega <- Hdiaginv %*% H
      Omegaeig <- eigen(Omega)
      U <- Omegaeig$vectors
      Uinv <- solve(U)
      EV <- Omegaeig$values
      UinvHdiaginv <- Uinv %*% Hdiaginv
      
      LMHinv <- U %*% diag(1 / (EV + damp)) %*% UinvHdiaginv
      xn <- as.numeric(x - LMHinv%*%df)
      
      diff <- logtar(xn, lambda)$logtarget - f
      if(diff > 0){
        damp <- damp * (1/10)
      } else{
        dampiter <- 0
        while(diff <=0 & dampiter < dampitermax){
          damp <- 10 * damp
          LMHinv <- U %*% diag(1 / (EV + damp)) %*% UinvHdiaginv
          xn <- as.numeric(x - LMHinv%*%df)
          diff <- logtar(xn, lambda)$logtarget - f
          dampiter <- dampiter + 1
        }
      }
      dist <- sqrt(sum((xn-x)^2))
      if(dist < tol) conv <- 1
      x <- xn
      iter <- iter + 1
    }
    outlist <- list("par" = x, "conv" = conv)
    return(outlist)
  }
  
  latinit <- rep(0, (K+p+1))
  
  lamblow <- 1
  trylat <- try(LM(x = latinit, lambda = lamblow)$par,silent = TRUE)
  while(substr(trylat[1],1,5)=="Error"){
    lamblow <- lamblow + 1
    trylat <- try(LM(x = latinit, lambda = lamblow)$par,silent = TRUE)
  }
  latinit <- LM(x = latinit, lambda = lamblow)$par
  
  logpv <- function(v){
    lambda <- exp(v)
    latopt <- LM(x = latinit, lambda = lambda)$par
    logtarcall <- logtar(latopt, lambda = lambda, hess = TRUE)
    Sigopt <- solve((-1) * logtarcall$D2logtarget)
    val <- 0.5 * log(det(Sigopt)) + logtarcall$loglik -
      0.5 * sum((latopt * Q(lambda))%*%latopt) +
      0.5 * (K+nu) * v -
      (0.5 * nu + a_kappa) * log(exp(v) * 0.5 * nu + b_kappa) 
    return(val)
  }
  
  #----- Exploration of penalty over grid
  loglambdas <- seq(0, 5, length = 100)
  lamb_grid <- round(10 ^ loglambdas)
  logplamb_val <- c()
  logplamb_val[1] <- logpv(lamb_grid[1])
  for(j in 2:length(lamb_grid)){
    logplamb_val[j] <- logpv(lamb_grid[j])
    if(logplamb_val[j]<logplamb_val[j-1] | is.na(logplamb_val[j])) 
      break
  }
  jidx <- length(logplamb_val) - 1
  vopt <- lamb_grid[jidx]
  
  lambmap <- exp(vopt)
  Lap <- LM(x = latinit, lambda = lambmap)
  latenthat <- as.numeric(Lap$par)
  Sigmap <-  diag(diag(solve((-1) * logtar(latenthat, lambda = lambmap, hess = TRUE)$D2logtarget)))
  Ptilde <- matrix(0,nrow = (K+p+1), ncol = (K+p+1))
  Ptilde[1:K,1:K] <- P
  thetahat <- latenthat[1:K]
  betahat <- latenthat[(K+1):(K+p)]
  frailhat <- exp(latenthat[(K+p+1)])
  lambdahat <- lambmap

  # 90% and 95% CrI
  CI90low <- c(betahat - qnorm(0.95) * sqrt(diag(Sigmap)[(K+1):(K+p)]),
               exp(latenthat[(K+p+1)] - qnorm(0.95) * sqrt(diag(Sigmap)[(K+p+1)])))
  CI90up <- c(betahat + qnorm(0.95) * sqrt(diag(Sigmap)[(K+1):(K+p)]),
              exp(latenthat[(K+p+1)] + qnorm(0.95) * sqrt(diag(Sigmap)[(K+p+1)])))
  CI90 <- cbind(CI90low, CI90up)
  rownames(CI90) <- c("beta0", "beta1", "frail")
  
  CI95low <- c(betahat - qnorm(0.975) * sqrt(diag(Sigmap)[(K+1):(K+p)]),
               exp(latenthat[(K+p+1)] - qnorm(0.975) * sqrt(diag(Sigmap)[(K+p+1)])))
  CI95up <- c(betahat + qnorm(0.975) * sqrt(diag(Sigmap)[(K+1):(K+p)]),
              exp(latenthat[(K+p+1)] + qnorm(0.975) * sqrt(diag(Sigmap)[(K+p+1)])))
  CI95 <- cbind(CI95low, CI95up)
  rownames(CI95) <- c("beta0", "beta1", "frail")
  
  
  outlist <- list(thetahat = thetahat, betahat = betahat,
                  frailhat = frailhat, lambdahat = lambdahat,
                  CI90 = CI90, CI95 = CI95,
                  tmax = tmax, K = K)
  
  return(outlist)
}
  
  
  