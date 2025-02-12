###Augmented drift
drift = function(x)
{
  x1 = x[1]
  x2 = x[2]
  
  if( x1>0 && x2>0)
  {
    return(-x)
  }
  else
  {
    return(c(-Inf, -Inf))
  }
}

#Calculate projetction
projection  = function(u)
{
  if(all(u[1] < 0, na.rm = TRUE) && all(u[2] < 0, na.rm = TRUE))
  {
    return(c(0,0))
  } else if (u[1]>0 & u[2]<0)
  {
    return(c(u[1],0))
  } else if (u[1]<0 & u[2]>0)
  {
    return(c(0,u[2]))
  }
  else{
    return(u)
  }
}

###Calculate Preconditioning matrix 
C = function(x)
{
  u = projection(x)
  s = matrix(u-x , nrow = 2, ncol = 1)
  den = t(s)%*%s
  den = den[1,1]
  if(s[1]== 0 & s[2]==0 )
  {
    I = diag(1,nrow = 2 , ncol = 2)
    return(I)
  }
  else
  {
    P = s%*% t(s)
    P = P/den
    return(P)
  }
  
}

###### Probability function
probf = function(z , v, pmat,dt)
{
  den = e^(-sqrt(2*dt)*v%*%t(C))
  rtn = 1/(1+den)
  return(rtn)
}


PrAugUBA = function(x0, dt, N)
{
  out = matrix(nrow = N , ncol = 2)
  out[1 ,] = x0
  
  for(t in 2:N)
  {
    x = out[t-1, ]
    v = matrix(rnorm(n = 2), ncol = 1)
    U = runif(2)
    
    mu = matrix(drift(x), nrow = 2)
    mat = C(x)
    C_i = mat%*%mu
    
    den = exp(-sqrt(2*dt)*v*C_i)
    p = 1/(1+den)
    
    z = ifelse(U <= p , v , -v)
    out[t , ] = out[t-1 , ] + sqrt(2*dt)*mat%*%z
    print('Current states: ')
    print(out[t, ])
    
  }
  return(out)
}

try1 = PrAugUBA(x0 = c(100,100), dt = 1e-2 , N = 1e3)
plot(try1,pch = 16, xlab = 'x', ylab = 'y')
