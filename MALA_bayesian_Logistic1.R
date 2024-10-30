set.seed(123)
##################
##Data Generation
##################
library(mvtnorm)
n = 1e3
k= 5
h = 1

X = matrix(rnorm(n*K, mean = 0, sd = 3), nrow =n , ncol = k)
betastar = matrix(rmvnorm(n = 1 , mean = rep(0,k)), nrow = k , ncol =1 )
p = 1/(1 + exp(- X%*% betastar))
Y = matrix(rbinom(n , size = 1,p ), nrow = 1000, ncol = 1)

data = cbind(Y,X)
#####################
###MALA Updates
#####################


##Posterior density calculated at u
logf = function(beta)
{
  g1 = sum(Y*log(p)+(1-Y)*log(1-p))
  g2 = -t(beta)%*%beta/2
  return(g1+g2)
}

grad_log_f = function(beta)
{
  
  den = 1 + exp(- X %*% beta)
  num = Y %*% t(X)%*% exp( - X%*% beta) - (1-Y)%*%X
  rtn = -beta + colSums(num/den)
  return(rtn)
}

## Kernel density calculated at y given x 
logQ = function(x,y)
{
  rtn = dmvnorm( y, mean = (x - h*grad_log_f(x)/2) , sigma = diag(K), log = TRUE)
  return(rtn)
}

ratio = function(x,y)
{
  num = logf(y)+logQ(y,x)
  den = logf(x)+logQ(x,y)
  logr = min( 0, num-den)
  return(logr)
}

MALA = function(x0,n)
{
  x = matrix(nrow = n , ncol = K)
  x[1, ] = x0
  for( i in 1:(n-1))
  {
    y = rmvnorm(n=1 , mean = (x[i, ] - h*grad_log_f(x[i,])/2) , sigma = diag(K))
    U = runif(1) 
    alpha = ratio(x[i,],y)
    x[i+1,] = ifelse( log(U)<alpha , y , x[i, ])
  }
  return(x)
}
Sys.time()
beta_pst = MALA(beta, 1e5)
colMeans(beta_pst)
beta
Sys.time()


