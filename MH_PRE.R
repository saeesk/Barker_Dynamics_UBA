set.seed(123)
#######################
#### Data Generation
#######################
#Data Generation 
lambdastar = 1e-2
I =50
J = 5
sigma_l = 1e-4
eta = rnorm(n = I , mean = lambdastar , sd = 1)
y = matrix(nrow = I , ncol = J)
for(i in 1:I)
{
  y[i ,] = rpois(n = J , lambda = exp(eta[i]))
}

#####################
###MALA Updates
#####################
library(mvtnorm)
##Posterior density calculated at u
logf = function(x)
{
  lambda = x[1]
  eta = x[-1]
  
  t1 = sum(J*exp(eta)) 
  t2 = 0.5*sum((eta-lambda)^2) + sigma_l*lambda
  eta = as.matrix(eta ,nrow =I,ncol=1 )
  t3 = sum(t(y)%*%eta) 
  
  return(t1+t2-t3)
}

## Kernel density calculated at y given x 
logQ = function(x,y)
{
  rtn = dmvnorm( y, mean = x , sigma = diag(I+1), log = TRUE)
  return(rtn)
}

ratio = function(x,y)
{
  num = logf(y)+logQ(y,x)
  den = logf(x)+logQ(x,y)
  logr = min( 0, num-den)
  return(logr)
}

MH = function(x0,n)
{
  x = matrix(nrow = n , ncol = (I+1))
  x[1, ] = x0
  for( i in 1:(n-1))
  {
    #current 
    current = x[i, ]
    #proposal
    prop = rmvnorm(n=1 , mean = current , sigma = diag(I+1))
    
    U = runif(1) 
    alpha = ratio(current,prop)
    
    x[(i+1),] = ifelse( log(U)<alpha , prop , current)
    
    print(paste('step ' ,i ,'done!'))
  }
  return(x)
}
Sys.time()
temp = c(lambdastar,eta)
lamba_pst = MALA(temp, 1e4)[ ,1]
plot(density(lamba_pst))
Sys.time()

