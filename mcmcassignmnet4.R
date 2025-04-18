###########################################################
####Beta updated by MH, Lambda and Tau using augmentation
###########################################################
log_f = function(y, X, beta, lambda, tau) {
  n = nrow(X)
  p = ncol(X)
  beta = as.numeric(beta)
  
  l = sum(lambda > 0)
  t = tau > 0 
  flag = l + t 
  
  if (flag == p + 1) {
    # Log-likelihood
    power = X %*% beta
    log.lik = sum(y * power - exp(power))
    
    # Prior for beta
    log.beta = sum(dnorm(beta, mean = 0, sd = tau * lambda, log = TRUE))
    
    # Half-Cauchy priors for lambda and tau
    log.lambda = sum(-log(1 + lambda^2))
    log.tau = -log(1 + tau^2)
    
    log.post = log.lik + log.beta + log.lambda + log.tau
  } else {
    log.post = -Inf
  }
  
  return(log.post)
}

poissonHorse <- function(y, X, n_iter, beta_start, lambda_start, tau_start) {
  n <- length(y)
  p <- ncol(X)
  
  # Initialize storage
  beta <- matrix(0, nrow = n_iter, ncol = p)
  lambda <- matrix(0, nrow = n_iter, ncol = p)
  tau <- numeric(n_iter)
  
  # Augmented variables
  nu <- rep(1, p)
  xi <- 1
  
  # Initial values
  beta_current <- beta_start
  lambda_current <- lambda_start
  tau_current <- tau_start
  h <- 1e-3
  
  for (i in 1:n_iter) {
    lambda2 <- lambda_current^2
    tau2 <- tau_current^2
    
    # 1. Update beta via MH
    beta_prop <- beta_current + rnorm(p, mean = 0, sd = h)
    
    log_ratio <- log_f(y, X, beta_prop, lambda_current, tau_current) - 
      log_f(y, X, beta_current, lambda_current, tau_current)
    
    if (log(runif(1)) < log_ratio) {
      beta_current <- beta_prop
    }
    
    # 2. Update lambda_j^2
    lambda2_new <- numeric(p)
    for (j in 1:p) {
      shape <- 1 + 0.5
      rate <- beta_current[j]^2 / (2 * tau2) + 1 / nu[j]
      lambda2_new[j] <- 1 / rgamma(1, shape = shape, rate = rate)
    }
    lambda_current <- sqrt(lambda2_new)
    
    # 3. Sample nu_j
    for (j in 1:p) {
      nu[j] <- 1 / rgamma(1, shape = 2, rate = 1 + 1 / lambda2_new[j])
    }
    
    # 4. Sample tau^2
    shape_tau <- (p + 1) / 2
    rate_tau <- 1 / xi + sum(beta_current^2 / lambda2_new) / 2
    tau2 <- 1 / rgamma(1, shape = shape_tau, rate = rate_tau)
    tau_current <- sqrt(tau2)
    
    # 5. Sample xi
    xi <- 1 / rgamma(1, shape = 2, rate = 1 + 1 / tau2)
    
    # Store samples
    beta[i, ] <- beta_current
    lambda[i, ] <- lambda_current
    tau[i] <- tau_current
  }
  
  return(list(beta = beta, lambda = lambda, tau = tau))
}

####Testing if function works correctly
###Data
dat = read.csv('betting.csv')
X = as.matrix(dat[ ,-1])
y = dat[,1]
p = dim(X)[2]
###Starting values
fit_glm = glm(y ~ X - 1, family = poisson)
beta_start = coef(fit_glm)
lambda_start = abs(rt(p, df = 1))
tau_start=1
time = system.time({
  out = poissonHorse(y, X, n_iter = 1e5, beta_start, lambda_start, tau_start)
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

