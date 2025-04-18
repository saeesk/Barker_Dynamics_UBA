#####################################################
#### Componentwise MH sampler with:
#### - Gaussian proposal for beta 
#### - Lomax proposals for lambda_j (updated separately) 
#### - Lomax proposal for tau
#####################################################

# Log-posterior function 
log_f <- function(y, X, beta, lambda, tau) {
  n <- nrow(X)
  p <- ncol(X)
  
  if (all(lambda > 0) && tau > 0 && all(is.finite(lambda))) {
    eta <- X %*% beta
    log_lik <- sum(y * eta - exp(eta))  # Poisson log-likelihood
    
    log_prior_beta <- sum(dnorm(beta, mean = 0, sd = tau * lambda, log = TRUE))
    log_prior_lambda <- sum(-log(1 + lambda^2))  # Half-Cauchy
    log_prior_tau <- -log(1 + tau^2)             # Half-Cauchy
    
    return(log_lik + log_prior_beta + log_prior_lambda + log_prior_tau)
  } else {
    return(-Inf)
  }
}

# Lomax (Pareto type II) functions
r_lomax <- function(alpha, lambda, s = 1) {
  lambda * ((1 - runif(s))^(-1/alpha) - 1)
}

d_lomax <- function(x, alpha, lam) {
  ifelse(x > 0, 
         log(alpha) - log(lam) - (alpha + 1) * log(1 + x/lam),
         -Inf)
}

poissonHorse <- function(y, X, n_iter, beta_start, lambda_start, tau_start) {
  p <- ncol(X)
  beta <- matrix(0, ncol = p, nrow = n_iter)
  lambda <- matrix(0, ncol = p, nrow = n_iter)
  tau <- numeric(n_iter)
  
  beta[1,] <- beta_start
  lambda[1,] <- lambda_start
  tau[1] <- tau_start
  
  proposal_sd_beta <- 0.02   # Tune this
  alpha_lomax <- 0.5         # Shape parameter for Lomax proposals
  lambda_lomax <- 1.0        # Scale parameter for Lomax proposals
  
  for (i in 2:n_iter) {
    # Current values
    beta_curr <- beta[i-1,]
    lambda_curr <- lambda[i-1,]
    tau_curr <- tau[i-1]
    
    ### 1. Update beta (vectorized)
    beta_prop <- beta_curr + rnorm(p, 0, proposal_sd_beta)
    log_ratio <- log_f(y, X, beta_prop, lambda_curr, tau_curr) - 
      log_f(y, X, beta_curr, lambda_curr, tau_curr)
    
    if (is.finite(log_ratio) && log(runif(1)) < log_ratio) {
      beta[i,] <- beta_prop
    } else {
      beta[i,] <- beta_curr
    }
    
    ### 2. Component-wise update for lambda (vectorized implementation)
    # Propose new values for all lambdas
    lambda_prop <- r_lomax(alpha_lomax, lambda_lomax, p)
    log_target_prop <- log_f(y, X, beta[i,], lambda_prop, tau_curr)
    log_target_curr <- log_f(y, X, beta[i,], lambda_curr, tau_curr)
      
    # Proposal densities (only for the j-th component)
    log_q_curr <- d_lomax(lambda_curr, alpha_lomax, lambda_lomax)
    log_q_prop <- d_lomax(lambda_prop, alpha_lomax, lambda_lomax)
    log_ratio <- (log_target_prop - log_target_curr) + (log_q_curr - log_q_prop)
    for(j in 1:p)
    {
      lambda[i,j ] = ifelse( log(runif(1)) < log_ratio[j] , lambda_prop[j] , lambda_curr[j])
    }
    
    ### 3. Update tau
    tau_prop <- r_lomax(alpha_lomax, lambda_lomax)
    log_target_prop <- log_f(y, X, beta[i,], lambda[i,], tau_prop)
    log_target_curr <- log_f(y, X, beta[i,], lambda[i,], tau_curr)
    
    log_q_curr <- d_lomax(tau_curr, alpha_lomax, lambda_lomax)
    log_q_prop <- d_lomax(tau_prop, alpha_lomax, lambda_lomax)
    
    log_ratio <- (log_target_prop - log_target_curr) + (log_q_curr - log_q_prop)
    
    if (is.finite(log_ratio) && log(runif(1)) < log_ratio) {
      tau[i] <- tau_prop
    } else {
      tau[i] <- tau_curr
    }
    
    if (i %% 1000 == 0) cat("Iteration:", i, "\n")
  }
  
  list(beta = beta, lambda = lambda, tau = tau)
}

# Run the analysis
dat <- read.csv('betting.csv')
X <- as.matrix(dat[,-1])
y <- dat[,1]
p <- ncol(X)

fit_glm <- glm(y ~ X - 1, family = poisson)
beta_start <- coef(fit_glm)
lambda_start <- abs(rt(p, df = 1))
tau_start <- 1 

time <- system.time({
  out <- poissonHorse(y, X, n_iter = 1e3, beta_start, lambda_start, tau_start)
})

# Compute ESS and efficiency
library(mcmcse)
ess_beta <- min(ess(out$beta))
ess_lambda <- min(ess(out$lambda))
ess_tau <- ess(out$tau)
minEss <- min(ess_beta, ess_lambda, ess_tau)
eval <- minEss/time[3]
eval