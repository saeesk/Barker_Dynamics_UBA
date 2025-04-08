####################################
##Truncated Poisson Random Effects 
#####################################
library(extraDistr)
set.seed(123)
#Data Generation 
lambdastar = 5
I =1e3
J = 10
sigma_0 = 16.5 #hyperparameter
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
  
  lambda = x[1]
  eta = x[-1]
  #Augmented gradient for lambda
  dU1 = ifelse(lambda>0 , I/lambda - sigma_0 - sum(eta), Inf )
  #Augmented Gradient for eta
  dU2 = ifelse(eta>0, -lambda - J*exp(eta)+ rowSums(y), Inf)
  
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
    if(i%% 1e3 ==0)
    {
      print(i)
    }
    
  } 
  #pst_lambda= x[ ,1]
  #return(pst_lambda)
  return(x)
}

#x0 = -c(rexp(1, rate = lambdastar),eta0)
x0 = -c(lambdastar,eta0)
#x0 = rep(10,51)
try1 = aug_UBA(x0, dt = 1e-3, N = 1e4)

l = try1[ ,1]
plot.ts(l, main = 'lambda')
abline(h = lambdastar, col = 'blue', lty=2)
plot(density(l),
     main = paste('Started outside, r = ', sum(l<0)/1e4))
abline(v = lambdastar, col = 'blue', lty = 2)

s = rdunif(n = 10 , min = 1, max = I)
for(i in 1: length(s))
{
  para = try1[ ,s[i]]
  plot.ts(para, main = paste('eta_',s[i]))
  abline(h = eta0[s[i]], col = 'blue', lty=2)
  plot(density(para),
       main = paste('eta_',s[i], ' r = ', sum(para<0)/1e4))
  abline(v = eta0[s[i]], col = 'blue', lty = 2)
}


