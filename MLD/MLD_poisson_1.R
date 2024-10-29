#########################
##PRE positive parameter
#########################
set.seed(0809204)
#Data Generation 
lambdastar = 1e-4
I =50
J = 5
sigma_l = 1e-4
eta = rnorm(n = I , mean = lambdastar , sd = 1)
y = matrix(nrow = I , ncol = J)
for(i in 1:I)
{
  y[i ,] = rpois(n = J , lambda = exp(eta[i]))
}



###################
###MLD for Poisson 
####################
library(mvtnorm)

dU = function(x)
{
  lambda = x[1]
  dU1 = -sum(eta - lambda)
  eta = x[-1]
  dU2 = J*exp(eta) - rowSums(y)
  dU = c(-dU1,-dU2)
  return(dU)
}
MLD = function(x, dt, N)
{
  Y = matrix(ncol = I+1, nrow =  N)
  X = matrix(ncol = I+1, nrow =  N)
  X[1 , ] = x 
  
  for( i in 1: (N-1))
  {
    nu = rmvnorm(I+1 , mean = rep(0 , nrow = I+1))
    Y[i , ] = 1 + log(X[i, ])
    nabU = dU(X[i, ])
    Y[(i+1), ] = Y[i, ]- dt %*% X[i, ]%*% (nabU - 1/(X[i, ]))
    X[(i+1) , ] = exp(Y[(i+1), ]- 1)
  }
  
  rtn = X[-1, ]
  return(rtn)
} 


x0 = c(rexp(1, rate =1/ lambdastar),eta)
r1= MLD(x0, dt = 1e-3 , 1e3)

