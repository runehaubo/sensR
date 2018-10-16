# karolina_thesis_functions.R - Functions from Karolina Stachlewska's master 
# thesis, appendix E.

#' Likelihood function for the chance-corrected beta-binomial model
#'
#' @param a parameter
#' @param b parameter
#' @param x number of correct responses
#' @param k number of replications
#' @param p0 guessing probability
#'
#' @rdname cbb
#' @export
cbb <- function(a, b, x, k, p0){
  ((1 - p0)^k) / beta(a,b) * choose(k, x) *
    sum( sapply( c(0:x), function(i){ choose(x,i) * (p0/(1 - p0))^(x-i) * beta(a+i, k+b-x) } ) )
}

#' Likelihood function for the chance-corrected beta-binomial model
#'
#' @param parameters mu and gamma parameter vector
#' @param y vector of correct answers
#' @param k number of replications
#' @param p0 guessing probability
#'
#' @export
loglikcbb <- function(parameters, y, k, p0){
  mu <- parameters[1]
  gamma <- parameters[2]
  #points outside the parameter space
  if (mu <= 0 || mu >= 1 || gamma <= 0 || gamma >= 1) return(Inf)
  #reparameterize
  a <- mu * (1/gamma - 1)
  b <- (1 - mu) * (1/gamma - 1)
  #calculate negative log-likelihood
  - sum( sapply(y, function(j){ log( cbb(a, b, j, k, p0) ) } ) )
}

#' Optimize the chance-corrected beta-binomial likelihood
#'
#' @param y vector of correct answers
#' @param k number of replications
#' @param p0 guessing probability
#' @param init vector of initial parameters
#' @param max.try maximum number of re-tries
#'
#' @export
optimcbb <- function(y, k, p0, init = c(0.5,0.2), max.try = 2){
  #check the arguments
  stopifnot(all(y%%1 == 0), all(y >=0), all(y <= k), p0 < 1, p0 > 0, max.try >= 0)
  #check initial values passesd to the optimizer
  if(any(init <= 0) | any(init >= 1)){
    init <- c(0.5,0.2)
    message("Both initial values have to satisfy 0 < initial_value < 1.
            \nChanging to default: (0.5,0.2)")
  }
  #optimize
  bbRho <- bbEnvir(parent.frame(), X = as.matrix(cbind(y, k)), 
                   corrected = TRUE, pGuess = p0, start = init)
  result <- NA
  i <- 0
  while ( i <= max.try && (is.na(result[1]) || result$convergence != 0) ){
    if (i != 0){
      init <- runif(n = 2, min = 1e-6, max = 1 - 1e-6)
      message(paste("A convergence problem or an error occurred.
                    \nChanging initial values to randomly generated: (",
                    init[1], ", ", init[2], ")", sep = ""))
    }
    suppressWarnings(try(
      result <- optim(par = init, fn = loglikcbb, 
                      y = y, k = k, p0 = p0, method = "L-BFGS-B", hessian = FALSE, 
                      lower = bbRho$lbounds, upper = bbRho$ubounds, 
                      control = list(parscale = c(0.01, 0.01))), silent = TRUE))
    i <- i + 1
  }
  #display results
  if (is.na(result[1])){
    stop("Optimization could not be performed.")
  } else {
    if (result$convergence != 0){
      warning("Algorithm did not converge. Estimates are not reliable.")
    }
    thelist <- list("mu" = result$par[1], "gamma" = result$par[2],
                    "neg.loglik" = result$value, "conv" = (result$convergence == 0))
    return(thelist)
  }
  }


#' Plot likelihood confidence region for mu and gamma
#'
#' @param y vector of correct answers
#' @param k number of replications
#' @param p0 guessing probability
#' @param init vector of initial parameters
#' @param max.try maximum number of re-tries
#' @param level confidence level
#' @param length.out number of grid points at each axis 
#' @param p df for the chi-square approximation (between 1 and 2)
#'
#' @export
plotcbb <- function(y, k, p0, level = 0.95, length.out = 101, p = 2, init = c(0.5,0.2), max.try = 2){
  #check the arguments
  stopifnot(level < 1, level > 0, length.out > 2, (p >= 1 & p <= 2))
  #optimize
  optimization <- optimcbb(y = y, k = k, p0 = p0, init = init, max.try = max.try)
  if (!optimization$conv){
    stop("Algorithm did not converge. Results can not be used to find confidence region.")
  }
  thelist <- optimization[1:3] #list of results
  loglikbound <- thelist$neg.loglik + qchisq(level, p)/2 #threshold value
  thelist$loglik.bound <- loglikbound #add to the results
  #construct the grid
  mu_vals <- seq(from = 0, to = 1, length.out = length.out)
  gamma_vals <- seq(from = 0, to = 1, length.out = length.out)
  loglik_vals <- matrix(data = NA, length(mu_vals), length(gamma_vals))
  for (i in 2:(length.out - 1)){
    for (j in 2:(length.out - 1)){
      loglik_vals[i,j] <- loglikcbb(c(mu_vals[i], gamma_vals[j]), y, k, p0)
    }
  }
  #choose the points within the likelihood region
  good_points <- (loglik_vals < loglikbound)
  #plot the points that are within the region and the ones that give NAs
  par(mfrow = c(1,1), mar = c(2.5,2.5,1.5,0.5), mgp = c(1.2,0.2,0), cex.lab = 1.2)
  plot(NA, NA, type = 'n', xlim = c(0,1), ylim = c(0,1.08),
       xlab = expression(mu), ylab = expression(gamma),
       main = paste(level*100, "% likelihood region", sep = ""))
  for (i in 1:length(mu_vals)){
    for (j in 1:length(gamma_vals)){
      if (!is.na(good_points[i,j]) & good_points[i,j]){
        points(mu_vals[i], gamma_vals[j], pch = 20, col = "green", cex = 0.3) #within the region
      }
      if (is.na(good_points[i,j])){
        points(mu_vals[i], gamma_vals[j], pch = 20, col = "red", cex = 0.3) #giving NA
      }
    }
  }
  points(thelist$mu, thelist$gamma, pch = 20, col = "forestgreen") #MLE
  legend("topleft", legend = c("point estimate", "within the region", "not computable"),
         xpd = TRUE, bty = 'n', horiz = TRUE, cex = 0.8, pt.cex = c(2,5,5),
         fill = c("forestgreen", "green", "red"), border = c("forestgreen", "green", "red"))
  #give results
  return(thelist)
}

#' @rdname cbb
#' @param mu mu parameter
#' @param gamma gamma parameter
#' @param y vector of correct answers
loglikcbb.gamma <- function(mu, gamma, y, k, p0){ #to optimize with respect to mu
  #points outside the parameter space
  if (mu >= 1 || mu <= 0 || gamma >= 1 || gamma <= 0) return(Inf)
  #reparameterize
  a <- mu * (1/gamma - 1)
  b <- (1 - mu) * (1/gamma - 1)
  #calculate negative log-likelihood
  - sum( sapply(y, function(j){ log( cbb(a, b, j, k, p0) ) } ) )
}

#' @rdname cbb
loglikcbb.mu <- function(gamma, mu, y, k, p0){ #to optimize with respect to gamma
  #points outside the parameter space
  if (gamma >= 1 || gamma <= 0 || mu >= 1 || mu <= 0) return(Inf)
  #reparameterize
  a <- mu * (1/gamma - 1)
  b <- (1 - mu) * (1/gamma - 1)
  #calculate negative log-likelihood
  - sum( sapply(y, function(j){ log( cbb(a, b, j, k, p0) ) } ) )
}

#' Profile CBB likelihood function
#'
#' @param y vector of correct answers
#' @param k number of replications
#' @param p0 guessing probability
#' @param init vector of initial parameters
#' @param max.try maximum number of re-tries
#' @param level confidence level
#' @param length.out number of grid points at each axis 
#' @param plot.profile logical; plot the profile likelihoods?
#'
#' @export
profilecbb <- function(y, k, p0, plot.profile = TRUE, level = 0.95,
                       length.out = 101, init = c(0.5,0.2), max.try = 2){
  #check the arguments
  stopifnot(level < 1, level > 0, length.out > 2, is.logical(plot.profile))
  #optimize
  optimization <- optimcbb(y = y, k = k, p0 = p0, init = init, max.try = max.try)
  if (!optimization$conv){
    stop("Algorithm did not converge. Results can not be used to find confidence region.")
  }
  loglik_bound <- optimization$neg.loglik + qchisq(level, 1)/2 #threshold value
  #construct the sequence for mu
  mu_vals <- sort( c(seq(from = 0, to = 1, length.out = length.out), optimization$mu) )
  mu_profile <- sapply( X = mu_vals, FUN = function(x){
    suppressWarnings(optimize(f = loglikcbb.mu,
                              interval = c(0,1), mu = x, y = y, k = k, p0 = p0)$objective)} )
  #construct the sequence for gamma
  gamma_vals <- sort( c(seq(from = 0, to = 1, length.out = length.out), optimization$gamma) )
  gamma_profile <- sapply( X = gamma_vals, FUN = function(x){
    suppressWarnings(optimize(f = loglikcbb.gamma,
                              interval = c(0,1), gamma = x, y = y, k = k, p0 = p0)$objective)} )
  #get approximated break points
  #mu
  mu_profile[mu_profile == Inf | mu_profile == -Inf] <- NaN
  spl_mu <- spline(mu_vals, mu_profile, xmin = 0, xmax = 1) #spline interpolation
  mu_lower <- mu_upper <- NA
  try( mu_lower <- approx(spl_mu$y[spl_mu$x < optimization$mu],
                          spl_mu$x[spl_mu$x < optimization$mu],
                          xout = loglik_bound)$y, silent = T ) #linear interpolation, lower bound
  try( mu_upper <- approx(spl_mu$y[spl_mu$x > optimization$mu],
                          spl_mu$x[spl_mu$x > optimization$mu],
                          xout = loglik_bound)$y, silent = T ) #linear interpolation, upper bound
  if (is.na(mu_lower)) mu_lower <- 0
  if (is.na(mu_upper)) mu_upper <- 1
  #gamma
  gamma_profile[gamma_profile == Inf | gamma_profile == -Inf] <- NaN
  spl_gamma <- spline(gamma_vals, gamma_profile, xmin = 0, xmax = 1) #spline interpolation
  gamma_lower <- gamma_upper <- NA
  try( gamma_lower <- approx(spl_gamma$y[spl_gamma$x < optimization$gamma],
                             spl_gamma$x[spl_gamma$x < optimization$gamma],
                             xout = loglik_bound)$y, silent = T ) #linear interpolation, lower bound
  try( gamma_upper <- approx(spl_gamma$y[spl_gamma$x > optimization$gamma],
                             spl_gamma$x[spl_gamma$x > optimization$gamma],
                             xout = loglik_bound)$y, silent = T ) #linear interpolation, upper bound
  if (is.na(gamma_lower)) gamma_lower <- 0
  if (is.na(gamma_upper)) gamma_upper <- 1
  #put the values in the list to return
  thelist <- list(mu_lower = NA, mu_upper = NA, gamma_lower = NA, gamma_upper = NA)
  thelist$mu_lower <- mu_lower
  thelist$mu_upper <- mu_upper
  thelist$gamma_lower <- gamma_lower
  thelist$gamma_upper <- gamma_upper
  #plot the profiles
  if (isTRUE(plot.profile)) {
    par(mfrow = c(2,1), mar = c(2.5,2.5,1,0.5), mgp = c(1.2,0.2,0), cex.lab = 1.2)
    layout(matrix(c(1,1,1,2,2,2,3), 7, 1, byrow = TRUE))
    #plot profile likelihood for mu
    plot(NA, NA, type = "n", xlab = expression(mu), ylab = "neg. log-lik",
         ylim = c(min(c(mu_profile, loglik_bound), na.rm = T),
                  max(c(mu_profile, loglik_bound), na.rm = T)),
         xlim = c(0,1), main = bquote("Profile likelihood for "* ~ mu))
    bool_mu <- (spl_mu$x <= mu_upper & spl_mu$x >= mu_lower)
    lines(spl_mu$x, spl_mu$y, col = "red", lwd = 2)
    lines(spl_mu$x[bool_mu], spl_mu$y[bool_mu], col = "green", lwd = 2)
    abline(h = loglik_bound, col = "gray60", lty = 2, lwd = 2)
    #plot profile likelihood for gamma
    plot(NA, NA, type = "n", xlab = expression(gamma), ylab = "neg. log-lik",
         ylim = c(min(c(gamma_profile, loglik_bound), na.rm = T),
                  max(c(gamma_profile, loglik_bound), na.rm = T)),
         xlim = c(0,1), main = bquote("Profile likelihood for "* ~ gamma))
    bool_gamma <- (spl_gamma$x <= gamma_upper & spl_gamma$x >= gamma_lower)
    lines(spl_gamma$x, spl_gamma$y, col = "red", lwd = 2)
    lines(spl_gamma$x[bool_gamma], spl_gamma$y[bool_gamma], col = "green", lwd = 2)
    abline(h = loglik_bound, col = "gray60", lty = 2, lwd = 2)
    #add legend
    plot(NA, NA, type = "n", xlim = c(0,1), ylim = c(0,1), xaxt = 'n', yaxt = 'n',
         ann = FALSE, bty = 'n')
    legend("center", legend = c("within", "outside", "threshold"), xpd = TRUE, bty = 'n',
           lty = c(1,1,2), lwd = c(2,2,2), col = c("green", "red", "gray60"), cex = 1.2,
           horiz = TRUE, title = paste(level*100, "% likelihood interval", sep = ""))
  }
  return(thelist)
}
