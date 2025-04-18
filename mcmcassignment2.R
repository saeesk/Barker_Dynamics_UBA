#######################################################
#### Here is a completely chatgpt solution : GIves NaN
#######################################################
poissonHorse <- function(y, X, n_iter, beta_start, lambda_start, tau_start) {
  # Dimensions
  n <- length(y)
  p <- ncol(X)
  
  # Storage
  beta <- matrix(0, nrow = n_iter, ncol = p)
  lambda <- matrix(0, nrow = n_iter, ncol = p)
  tau <- numeric(n_iter)
  
  # Initial values
  beta_curr <- beta_start
  lambda_curr <- lambda_start
  tau_curr <- tau_start
  
  # Prior parameters for inverse-gamma augmentation
  # Half-Cauchy(0, 1) ~ Inv-Gamma(1/2, 1/(2 * z)) where z ~ Inv-Gamma(1/2, 1)
  a0 <- 0.5
  b0 <- 1.0
  
  # Proposal standard deviation for Metropolis (tune this)
  proposal_sd <- 0.1
  
  for (iter in 1:n_iter) {
    ### --- 1. Update beta using Metropolis-Hastings ---
    beta_prop <- beta_curr + rnorm(p, 0, proposal_sd)
    
    loglik_curr <- sum(y * (X %*% beta_curr) - exp(X %*% beta_curr))
    logprior_curr <- sum(dnorm(beta_curr, 0, sqrt((tau_curr^2) * lambda_curr^2), log = TRUE))
    
    loglik_prop <- sum(y * (X %*% beta_prop) - exp(X %*% beta_prop))
    logprior_prop <- sum(dnorm(beta_prop, 0, sqrt((tau_curr^2) * lambda_curr^2), log = TRUE))
    
    log_accept_ratio <- (loglik_prop + logprior_prop) - (loglik_curr + logprior_curr)
    
    if (log(runif(1)) < log_accept_ratio) {
      beta_curr <- beta_prop
    }
    
    ### --- 2. Update lambda_j^2 (local scales) using Gibbs ---
    lambda2 <- lambda_curr^2
    for (j in 1:p) {
      shape <- a0 + 0.5
      rate <- (beta_curr[j]^2) / (2 * tau_curr^2) + b0
      lambda2[j] <- 1 / rgamma(1, shape = shape, rate = rate)
    }
    lambda_curr <- sqrt(lambda2)
    
    ### --- 3. Update tau^2 (global scale) using Gibbs ---
    shape <- a0 + p / 2
    rate <- sum(beta_curr^2 / lambda2) / 2 + b0
    tau2 <- 1 / rgamma(1, shape = shape, rate = rate)
    tau_curr <- sqrt(tau2)
    
    ### --- Store samples ---
    beta[iter, ] <- beta_curr
    lambda[iter, ] <- lambda_curr
    tau[iter] <- tau_curr
  }
  
  samples <- list("beta" = beta, "lambda" = lambda, "tau" = tau)
  return(samples)
}

################################################################################
#####Testing if function works correctly 
###Data 
dat = read.csv('betting.csv')
X = as.matrix(dat[ ,-1], ncol = p)
y = dat[,1]
###Starting values 
fit_glm <- glm(y ~ X - 1, family = poisson)
beta_start <- coef(fit_glm)
lambda_start = abs(rt(p, df = 1))
tau_start=1 
time = system.time({
  out <- poissonHorse(y, X, n_iter = 1e4, beta_start, lambda_start, tau_start)
})

dim(out[[1]])
dim(out[[2]])
dim(out[[3]])
library(mcmcse)
ess_beta = min(ess(out[[1]]))
ess_lambda = min(ess(out[[2]]))
ess_tau = min(ess(out[[3]]))
minEss = min(ess_beta, ess_lambda , ess_tau)
eval = minEss/time[3]
eval

