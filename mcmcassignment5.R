log_f <- function(y, X, beta, lambda, tau) {
  beta <- as.numeric(beta)
  l <- sum(lambda > 0)
  t <- tau > 0 
  flag <- l + t 
  
  if (flag == length(beta) + 1) {
    power <- X %*% beta
    log.lik <- sum(y * power - exp(power))
    log.beta <- sum(dnorm(beta, mean = 0, sd = tau * lambda, log = TRUE))
    log.lambda <- sum(-log(1 + lambda^2))
    log.tau <- -log(1 + tau^2)
    return(log.lik + log.beta + log.lambda + log.tau)
  } else {
    return(-Inf)
  }
}

poissonHorse <- function(y, X, n_iter, beta_start, lambda_start, tau_start) {
  n <- length(y)
  p <- ncol(X)
  
  beta <- matrix(0, nrow = n_iter, ncol = p)
  lambda <- matrix(0, nrow = n_iter, ncol = p)
  tau <- numeric(n_iter)
  
  # Initialization
  beta_current <- beta_start
  lambda_current <- lambda_start
  tau_current <- tau_start
  nu <- rep(1, p)
  xi <- 1
  h <- 1e-3
  
  lambda2 <- lambda_current^2
  tau2 <- tau_current^2
  
  for (t in 1:n_iter) {
    ### 1. Metropolis step for beta
    beta_prop <- beta_current + rnorm(p, mean = 0, sd = h)
    
    log_ratio <- log_f(y, X, beta_prop, sqrt(lambda2), sqrt(tau2)) -
      log_f(y, X, beta_current, sqrt(lambda2), sqrt(tau2))
    
    if (log(runif(1)) < log_ratio) {
      beta_current <- beta_prop
    }
    
    ### 2. Vectorized update of lambda^2
    shape_lambda <- 1.5
    rate_lambda <- beta_current^2 / (2 * tau2) + 1 / nu
    lambda2 <- 1 / rgamma(p, shape = shape_lambda, rate = rate_lambda)
    lambda_current <- sqrt(lambda2)
    
    ### 3. Vectorized update of nu
    nu <- 1 / rgamma(p, shape = 2, rate = 1 + 1 / lambda2)
    
    ### 4. Update tau^2 (scalar)
    shape_tau <- (p + 1) / 2
    rate_tau <- 1 / xi + sum(beta_current^2 / lambda2) / 2
    tau2 <- 1 / rgamma(1, shape = shape_tau, rate = rate_tau)
    tau_current <- sqrt(tau2)
    
    ### 5. Update xi (scalar)
    xi <- 1 / rgamma(1, shape = 2, rate = 1 + 1 / tau2)
    
    ### 6. Store current values
    beta[t, ] <- beta_current
    lambda[t, ] <- lambda_current
    tau[t] <- tau_current
    if(t %% 1000 == 0)
    {
      print(i)
    }
    
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
  out = poissonHorse(y, X, n_iter = 1e3, beta_start, lambda_start, tau_start)
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


