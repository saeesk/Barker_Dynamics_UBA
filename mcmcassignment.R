#######################################
####This code implements RWM (not tuned)
####################################### 

log_f = function(y, X , beta, lambda, tau) {
  # Log-likelihood Poisson regression
  n = dim(X)[1]
  p = dim(X)[2]
  beta = as.matrix(beta, nrow = 1, ncol = p)
  
  l = sum(lambda>0)
  t = s>0 
  flag = l+t 
  
  if(flag == p+1)
  {
    #Log(poisson likelihood)
    power = X %*% beta
    log.lik = sum(y * power - exp(power))  
    
    #Log-prior for beta 
    log.beta = sum(dnorm(beta, mean = 0, sd = tau * lambda, log = TRUE))
    
    # Log-prior for lambda (half-Cauchy(0,1))
    log.lambda = sum(-log(1 + lambda^2))
    
    # Log-prior for tau (half-Cauchy(0,1))
    log.tau = -log(1 + tau^2)
    
    # Sum all components (unnormalized log-f)
    log.post = log.lik + log.beta + log.lambda + log.tau
    
  }
  else
  {
    log.post = -Inf
  }
  
  return(log.post)
}


poissonHorse = function(y, X, n_iter, beta_start, lambda_start, tau_start)
{
  p = dim(X)[2]
  beta = matrix(0, ncol = p, nrow = n_iter)
  lambda = matrix(0, ncol = p, nrow = n_iter)
  tau = numeric(length = n_iter)
  n = dim(X)[1] 
  beta[1,] = beta_start
  lambda[1, ] = lambda_start
  tau[1] = tau_start
  h = 0.1
  
  ##Random Walk MH iterates
  for(i in 2:n_iter) {
    # Current state
    beta_curr = beta[i-1,]
    lambda_curr = lambda[i-1,]
    tau_curr = tau[i-1]
    
    #Draw proposals
    beta_prop = beta_curr + rnorm(p, mean = 0, sd = h)
    lambda_prop = lambda_curr + rnorm(p, mean = 0, sd = h)
    tau_prop = tau_curr + rnorm(1, mean = 0, sd = h)
   
    #log acceptance ratio
    log_ratio = log_f(y, X, beta_prop, lambda_prop, tau_prop) - 
      log_f(y, X, beta_curr, lambda_curr, tau_curr)
    
    #MH Filter
    if(log(runif(1))< log_ratio)
    {
      beta[i,] = beta_prop
      lambda[i,] = lambda_prop
      tau[i] = tau_prop
    }else {
      
      beta[i,] = beta_curr
      lambda[i,] = lambda_curr
      tau[i] = tau_curr
    }
    
    # Optional: Print progress
    if(i %% 100 == 0)
    {
      print(i) 
    }
  }
  
  
  samples = list("beta" = beta, "lambda" = lambda, "tau" = tau)
  return(samples)

}

##COMMENT OUT THE CODE LINES BELOW!!!
################################################################################
################################################################################
#####Testing if function works correctly 
# ###Data 
# dat = read.csv('betting.csv')
# X = as.matrix(dat[ ,-1])
# y = dat[,1]
# p = dim(X)[2]
# ###Starting values 
# fit_glm = glm(y ~ X - 1, family = poisson)
# beta_start = coef(fit_glm)
# lambda_start = abs(rt(p, df = 1))
# tau_start=1 
# time = system.time({
#   out = poissonHorse(y, X, n_iter = 1e5, beta_start, lambda_start, tau_start)
# })
# 
# dim(out[[1]])
# dim(out[[2]])
# dim(out[[3]])
# library(mcmcse)
# ess_beta = min(ess(out[[1]]))
# ess_lambda = min(ess(out[[2]]))
# ess_tau = min(ess(out[[3]]))
# minEss = min(ess_beta, ess_lambda , ess_tau)
# eval = minEss/time[3]
# eval
