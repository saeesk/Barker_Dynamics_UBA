##############################
#Poisson Random Effects Model
##############################
##Cacluate drift 
pois_drift = function(x)
{
  mu = x[1]
  eta = x[2:(I+1)]
  #derivative wrt mu
  dU1 = (mu/(sigma_mu)^2)- sum(eta) + (I+1)*mu 
  
  #derivative wrt eta_i
  dU = J*exp(eta) - rowSums(y) + eta - mu
  
  #derivative of U 
  dU = c(dU1, dU)
  
  return(-1*dU)
}

###EM scheme : ULA
ULA = function(x0,dt ,N)
{
  x = matrix(0 ,nrow = N,ncol= (I + 1))
  x[1 , ] = x0
  for( i in 2:N)
  {
    mu = pois_drift(x[(i-1) , ])
    x[i, ] =  x[(i-1), ]+ dt*mu+ sqrt(2*dt)*rnorm(1)
  }
  return(x)
}

###Skew-Symmetric Scheme : UBA
##Calculate the probability function for jumps 
probf = function(y, dt , v )
{
  mu = pois_drift(y)
  temp = -2*sqrt(dt)*v*mu/sqrt(2)
  return(1/(1 + exp(temp)))
}

##Do the Unadjusted Barker Procedure 
UBA= function(x0,dt,N)
{
  x = matrix(0 ,nrow = N,ncol= (I + 1))
  x[1 , ] = x0
  for( i in 2: N)
  {
    v = rnorm(I+1)
    U = runif(I+1)
    prob = probf(y = x[(i-1), ], dt = dt , v = v)
    b = ifelse(U<= prob , 1, -1)
    x[i , ] = x[(i-1) , ] + sqrt(2*dt)*b*v
    print(i)
  }
  return(x)
}

###############
###############
#Experiments
###Data generation 
mustar = 5
I = 50 
J = 5
sigma_mu = 10 
eta = rnorm(n = I , mean = mustar,sd = 1)
y = matrix(nrow = I , ncol = J) 
for(i in 1 :I)
{
  y[i ,] = rpois(n = J , lambda = eta[i] ) 
}

##Initialize 
x0 = c(mustar,eta) 
foo = 1:1e4
r1= ULA(x0,  dt = 0.025, N = 6e4)
r2 = UBA(x0,  dt = 0.025, N = 6e4)
ULA_mu = r1[-foo ,1]
plot(density(ULA_mu) , main = 'ULA Testing')
UBA_mu = r2[-foo,1]
plot(density(UBA_mu) , main = 'UBA Testing' , col = 'blue')

err = sum((UBA_mu - mustar)^2)/5e4

