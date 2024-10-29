set.seed(123)
####################
###Data Generation
#################### 
library(mvtnorm)
N = 1e4
K= 5
h = 1

X = matrix(rnorm(n = N*K, mean = 0, sd = 3), nrow =N , ncol = K)
beta = rmvnorm(n = 1 , mean = rep(0,K))
p = 1/(1 + exp(-X%*% t(beta)))
Y = rbinom(n =N , size = 1,p )

data = cbind(Y,X)
#####################
###MALA Updates
#####################

##Posterior density calculated at u
f = function(u)
{
  g = Y*log(p)+(1-Y)*log(1-p)
  g1 = exp(sum(g))
  g2 = exp(-u%*%t(u)/2)
  return(g1*g2)
}

## Kernel density calculated at y given x 
Q = function(y,x)
{
  rtn = dmvnorm( y, mean = (x - h*x/2) , sigma = diag(K))
  return(rtn)
}

ratio = function(x,y)
{
  num = f(y)*Q(y,x)
  den = f(x)*Q(x,y)
  r = min( 1, num/den)
  return(r)
}

MALA = function(x0,n)
{
  x = matrix(nrow = n , ncol = K)
  x[1, ] = x0
  for( i in 1:(n-1))
  {
    y = rmvnorm(n=1 , mean = (x[i, ] - h*x[i, ]/2) , sigma = diag(K))
    U = runif(1) 
    v = as.matrix(x[i,], nrow = 1 , ncol =5)
    num = f(y)*Q(y,v)
    den = f(t(v))*Q(t(v),y)
    alpha = min(1, num/den)
    x[i+1,] = ifelse(U<alpha , y , x[i, ])
  }
  return(x)
}
Sys.time()
beta_pst = MALA(beta, 1e5)
colMeans(beta_pst)
beta
Sys.time()



