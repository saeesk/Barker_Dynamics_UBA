####################################
##Truncated Poisson Random Effects 
#####################################
library(extraDistr)
set.seed(123)
#Data Generation 
lambdastar = 1e5
I =10
J = 5
sigma_0 = 10 #hyperparameter #sigma0=10
eta0 = rexp(n = I , rate = lambdastar)
y = matrix(nrow = I , ncol = J)
for(i in 1:I)
{
  y[i ,] = rpois(n = J , lambda = exp(eta0[i]))
}

###### Projection 
projection = function(x)
{
  x = ifelse(x>0 , x , 0)
  return(x)
}

##AgUBA 
aug_drift = function(x)
{
  u = projection(x)
  lambda = x[1]
  eta = x[-1]
  #Augmented gradient for lambda
  dU1 = ifelse(lambda - u[1]== 0 , I/lambda - sigma_0 - sum(eta),sign(u[1]-lambda)*Inf )
  #Augmented Gradient for eta
  dU2 = ifelse(eta -u[-1]== 0, -lambda - J*exp(eta)+ rowSums(y), sign(u[-1] -eta)*Inf)
  
  dU = c(dU1,dU2)
  return(dU)
}

aug_probf = function(u, dt , v )
{
  mu = aug_drift(u)
  temp = -sqrt(2*dt)*v*mu
  return(1/(1 + exp(temp)))
}

####Calculate preconditioning matrix
C = function(x)
{
  #print('x')
  #print(x)
  u = projection(x)
  ##print('#printing u')
  #print(u)
  s = matrix(u-x , nrow = I+1, ncol = 1)
  #print('#printing s')
  #print(s)
  den = t(s)%*%s
  den = den[1,1]
  
  foo = sum(s==0)
  #print(paste('foo: ', foo))
  
  if(foo>0)
  {
    #print('I am inside')
    I = diag(1,nrow = I+1 , ncol = I+1)
    rtn = I
    #return(I)
  }else
  {
    # print('I am outside the support')
    P = s%*% t(s)
    P = P/den
    rtn = P 
    #return(P)
  }
  #print('#printing rtn')
  #print(rtn)
  return(rtn)
  
}

aug_UBA= function(x0,dt,N)
{
  x = matrix(0 ,nrow = N,ncol= (I + 1))
  #initializing the first row 
  x[1 , ] = x0
  
  for( i in 2: N)
  {
    v = rnorm(I+1)
    U = runif(I+1)
    prob = aug_probf(u = x[(i-1), ], dt = dt , v = v)
    v = ifelse(U<= prob , v, -v)
    mat = C(x[(i-1), ])
    x[i , ] = x[(i-1) , ] + sqrt(2*dt)*mat%*%v
    if(i%% 1e4 ==0)
    {
      print(i)
    }
    
  } 
  #pst_lambda= x[ ,1]
  #return(pst_lambda)
  return(x)
}

####Calculate density
logf = function(x)
{ 
  s = sum(x>0)
  if(s== I+1)
  {
    lambda = x[1]
    ##print('lambda')
    #print(lambda)
    eta = matrix(x[-1],nrow = 1, ncol = I )
    rtn = I*log(lambda) - lambda*sigma_0 - lambda*sum(eta) -J*sum(exp(eta)) + sum(eta%*%y)
    return(rtn)
  }
  else
  {
    return(-Inf)
  }
}


###MH RWM
RWM = function(start, N , h)
{
  out = matrix(nrow = N , ncol = I+1)
  out[1, ] = start
  accept = 0
  for(t in 2:N)
  {
    prop = rnorm(I+1, sd = sqrt(h)) + out[t-1 , ]
    foo1 = logf(prop) 
    foo2 = logf(out[t-1, ]) 
    if(foo1 == -Inf  || foo2 == -Inf)
    {
      log.ratio = -Inf
    }else
    {
      log.ratio <- logf(prop) - logf(out[t-1,])
    }
    
    if(log(runif(1)) < log.ratio)
    {
      out[t, ] <- prop
      accept <- accept + 1
    } else{
      out[t, ] <- out[t-1, ]
    }
    #print('t')
    if(t%%10000 == 0)
    {
      print(paste('iter: ', t, 'accepted values: ', accept))
    }
    
  }
  ##print(paste("Acceptance rate is", accept/N ))
  print(paste("Acceptance probability: ", accept/N))
  return(out)
  
}



x0 = c(rexp(1, rate = lambdastar),eta0)
try1 = aug_UBA(x0, dt = 5e-2, N = 5e5)
l = try1[ ,1]
out = l[l>0]

##Run RWMH
x0 = c(rexp(1, rate = lambdastar),eta0)
out_MH = RWM(x0, N = 5e5, h = 1e-3)
plot(density(out_MH[,1]),
     main = paste('lamdastar = ',lambdastar,'r = ', sum(l<0)/5e5))
lines(density(out), col = 'blue')
legend('topright', legend = c('RWMH','AugUBA'), 
       col = c('black','blue'), lty = 1)
plot.ts(out, main ='AugUBA')
plot.ts(out_MH[,1], main= 'MH')

