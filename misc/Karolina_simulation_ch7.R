# Karolina_simulation_ch7.R
# Code from Karolina's thesis appendix F used in chapter 7.

para <- data.frame( a = rep(c(0.3,1,3,2,4), 8),
                    b = rep(c(2.7,2,6,2,0.5), 8),
                    p0 = c(rep(1/3, 20), rep(1/2, 20)),
                    n = c(rep( c(rep(12, 5), rep(18, 5), rep(30, 5), rep(100, 5)), 2)),
                    k = c(rep( c(rep(4, 5), rep(12, 5), rep(10, 5), rep(3, 5)), 2)) )
#load sensR package
library(sensR)
#number of simulated points for each parameter set
N <- 10000
for (j in 1:40){ #for each parameter set
  #set random seed
  set.seed(123)
  #get parameters
  a <- para$a[j]
  b <- para$b[j]
  p0 <- para$p0[j]
  n <- para$n[j]
  k <- para$k[j]
  mu <- a/(a+b)
  gamma <- 1/(a+b+1)
  #file name
  if (p0 == 1/2){
    filename <- paste(j, "_sim_", "a", a, "b", b, "n", n, "k", k, "N", N, "duotrio.csv", sep = "")
  } else{
    filename <- paste(j, "_sim_", "a", a, "b", b, "n", n, "k", k, "N", N, "triangle.csv", sep = "")
  }
  #create the data frame
  simdat <- matrix(NA, nrow = N, ncol = n)
  simdat <- data.frame(simdat)
  simdat$conv <- NA
  simdat$lik1.5 <- NA
  simdat$lik2 <- NA
  simdat$loglik.mle <- NA
  simdat$loglik.true <- NA
  simdat$loglik.true.mu <- NA
  simdat$loglik.true.gamma <- NA
  simdat$prof.mu.bool <- NA
  simdat$prof.gamma.bool <- NA
  simdat$wald.error <- NA
  simdat$wald.na <- NA
  simdat$wald.mu.95.bool <- NA
  simdat$wald.gamma.95.bool <- NA
  simdat$wald.mugamma.bool <- NA
  for (i in 1:N){
    #keep track of the progress
    print(i)
    #simulate probabilities
    sim_pd <- rbeta(n, a, b)
    sim_pc <- p0 + (1-p0) * sim_pd
    #simulate responses
    y <- rbinom(n = rep(1, n), size = rep(k, n), prob = sim_pc)
    simdat[i, 1:n] <- y
    #optimize
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # This line....
    optim_result <- optimcbb(y, k, p0)
    # ...replaces the following line (missing function)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    #save if converged
    simdat$conv[i] <- optim_result$conv
    if (optim_result$conv){
      #log-likelihood for true, underlying parameters
      loglik_true <- loglikcbb(c(mu, gamma), y, k, p0)
      #threshold values of log-likelihood
      loglik_bound15 <- optim_result$neg.loglik + qchisq(0.95, 1.5)/2 #1.5 d.f.
      loglik_bound2 <- optim_result$neg.loglik + qchisq(0.95, 2)/2 #2 d.f.
      #check if the true point is in the confidence region
      simdat$lik1.5[i] <- ( loglik_true < loglik_bound15 )
      simdat$lik2[i] <- ( loglik_true < loglik_bound2 )
      #get values of log-likelihood
      simdat$loglik.mle[i] <- optim_result$neg.loglik #in MLE
      simdat$loglik.true[i] <- loglikcbb(c(mu, gamma), y, k, p0) #in true parameters
      simdat$loglik.true.mu[i] <- optimize(f = loglikcbb.mu, interval = c(0,1), mu = mu,
                                           y = y, k = k, p0 = p0)$objective #for true mu
      simdat$loglik.true.gamma[i] <- optimize(f = loglikcbb.gamma, interval = c(0,1), gamma = gamma,
                                              y = y, k = k, p0 = p0)$objective #for true gamma
      #threshold value for the profile likelihood
      loglik_bound1 <- optim_result$lik + qchisq(0.95, 1)/2
      #profile mu and gamma
      lik_mu <- optimize(f = loglikcbb.mu, interval = c(0,1), mu = mu,
                         y = y, k = k, p0 = p0)$objective
      lik_gamma <- optimize(f = loglikcbb.gamma, interval = c(0,1), gamma = gamma,
                            y = y, k = k, p0 = p0)$objective
      #check if true points are inside confidence intervals
      simdat$prof.mu.bool[i] <- (lik_mu < loglik_bound1)
      simdat$prof.gamma.bool[i] <- (lik_gamma < loglik_bound1)
    }
    #fit the betabin from sensR
    responses <- cbind(y, rep(k,n))
    fit1 <- NA
    #try to use betabin function
    if (p0 == 1/2){
      try( fit1 <- betabin(responses, method = "duotrio"), silent = T )
    } else{
      try( fit1 <- betabin(responses, method = "triangle"), silent = T )
    }
    #record errors
    simdat$wald.error[i] <- (length(fit1) == 1) #check if it is NA, because it failed
    if (!is.na(fit1)){
      if ( is.na(summary(fit1)$coefficients[1,2]) ){
        simdat$wald.na[i] <- TRUE #record non-estimable SEs
      } else{
        simdat$wald.na[i] <- FALSE
        #check if parameters are within confidence intervals
        matr <- summary(fit1, level = 0.95)$coefficients
        simdat$wald.mu.95.bool[i] <- ((mu > matr[1,3]) && (mu < matr[1,4]))
        simdat$wald.gamma.95.bool[i] <- ((gamma > matr[2,3]) && (gamma < matr[2,4]))
        matr2 <- summary(fit1, level = 0.975)$coefficients
        simdat$wald.mugamma.bool[i] <- ((mu > matr2[1,3]) && (mu < matr2[1,4]) &&
                                          (gamma > matr2[2,3]) && (gamma < matr2[2,4]))
      }
    }
  } #end of the loop for N = 1:10,000
  #save the data frame
  write.csv(simdat, file = filename, row.names = FALSE)
} #end of the parameter loop 1:40
