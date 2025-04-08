##This coe implements MLD for simplex d dimensions
################################
###Mirrored Langevin Algorithm 
#################################
set.seed(123)
###Parmeters of dirichlet posterior
d = 2
nl = rep(3,d)
al= rep(0, d)
pars = nl+al

grad.h = function(x)
{
  rtn = log(x)- log( 1 - sum(x))
  return(rtn)
}

grad.w = function(y)
{
  a0 = sum(pars)
  pars.sub = pars[1:(d-1)]
  rtn = a0*(exp(y)/(1 + sum(exp(y)))) - pars.sub
  return(rtn)
}

grad.hstar = function(y)
{
  rtn = exp(y)/( 1 + sum(exp(y)))
  return(rtn)
}


MLD = function(x0, dt, N)
{
  x = matrix(0 , nrow = N, ncol = d-1) 
  y = matrix(0 , nrow = N , ncol = d-1)
  x[1, ] = x0[-d]
  for( i in 2: N) 
  {
    v = rnorm(d-1)
    y[(i-1), ] = grad.h(x[(i-1), ]) 
    y[i, ] = y[(i-1),  ] - dt*grad.w(y[(i-1), ]) + v*sqrt(2 *dt)
    x[i,] = grad.hstar(y[i, ])
    
  }
  x = cbind( x, 1-x)
  return(x)
}







x0 = rep(1/d, d)
#x0 = c(0.1,0.9)
dt = 1e-3
N <- 1e6
try = MLD(x0, dt= dt,N=N)
ind = 1
foo.x = seq(0,1, length.out = 1e3)
plot(foo.x , dbeta(foo.x , shape1 = pars[ind],shape2 = sum(pars) - pars[ind]), type ='l')
lines(density(try[,ind]), col = 'red')

