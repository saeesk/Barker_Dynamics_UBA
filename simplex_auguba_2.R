###AugUBA for Simplex
set.seed(123)
#Parameters of Dirichlet Distribution
#nl = c(1e5,10,10,rep(0,8))
nl = c(1,1)
al = rep(0.1 , 2) 
eps = 1e-10

##Calculate projection 
projsplx = function(u)
{
  ##Initialize
  y = sort(u)
  n = length(y)
  t = numeric(n-1)
  i = n-1
  
  ##Calculating t_hat
  while(i>= 1)
  {
    #Calculate t_i for i>=1
    t[i] = (sum(y[(i+1):n]) - 1)/(n-i)
    if(t[i] >= y[i])
    {
      t_hat = t[i]
      x = ifelse((u - t_hat) >0 , (u- t_hat) , 0)
      return(x)
    }
    else
    {
      #decrement
      i = i - 1
    }
  }
  
  #Calculate t_hat for i==0
  t_hat = (sum(y)-1)/n
  x = ifelse((u - t_hat) >0 , (u - t_hat) , 0)
  return(x)
  
}


##Calculate augmented drift
aug_drift = function(x)
{
  u = projsplx(x)
  dlogf = (nl+al -1)/x
  rtn = ifelse(u-x < eps , dlogf , sign(u-x)*Inf)
  return(rtn)
}

##Calculate probability function 
probf = function(u , dt, v)
{
  
  mu = aug_drift(u)
  temp = -sqrt(2*dt)*v*mu
  return(1/(1 + exp(temp)))
}


###Preconditioning matrix 
###Calculate preconditioning matrix
C = function(x)
{
  u = projsplx(x)
  s = matrix(u-x , nrow =2 , ncol = 1)
  
  #Projection matrix
  den = t(s)%*%s
  den = den[1,1]
  
  #foo = prod(s==0)
  foo = prod(s < eps)
  count = 0
  if(foo>0)
  {
    I = diag(1,nrow = 2 , ncol = 2)
    rtn = I
    count = count + 1
  }else{
    P = s%*% t(s)
    P = P/den
    rtn = P 
  }
  return(list(rtn ,count))
}

###AugUBA
aug_UBA = function(x0,dt,N)
{
  x = matrix(0 ,nrow = N,ncol=2)
  
  #initializing the first row 
  x[1 , ] = x0
  ind = rep(FALSE, N)
  counter = 0
  for( i in 2: N)
  {
    v = rnorm(2)
    U = runif(2)
    prob = aug_probf(u = x[(i-1), ], dt = dt , v = v)
    v = ifelse(U<= prob , v, -v)
    pack = C(x[(i-1), ])
    mat = pack[[1]]
    temp = pack[[2]]
    if(temp ==1)
    {
      counter = counter + 1
      ind[i] = TRUE
    }
    x[i , ] = x[(i-1) , ] + sqrt(2*dt)*mat%*%v
    if(i%% 1e3 ==0)
    {
      print(i)
    }
    
  } 
  return(list(x , ind, counter))
}
x0 = rep(1/2, 2)
try = aug_UBA(x0, dt = 1e-3 , N = 1e5) 
try[[3]]

try1 = try[[1]]
bd = numeric(1e4)
s = numeric(1e4)
##To check if every sample is inside the support or not 
for( i in 1: 1e4)
{
  r = try1[i,]
  bd[i] = sum(r>0 & r<1)
  s[i] = sum(r)
}
s

