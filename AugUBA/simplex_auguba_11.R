###AugUBA for Simplex
set.seed(123)

##Dimesion of the supprt
d = 2

#Parameters of Dirichlet Distribution
nl = rep(1000, d)#c(1e5,10,10,rep(0,8))
al= rep(1 , d) 
alpha <- al + nl
#eps = 1e-10

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
  
  u.sub <- u[1:(d-1)]
  x.sub <- x[1:(d-1)]
  
  dlogf = (nl+al -1)/x.sub
  rtn = ifelse( u.sub - x.sub == 0 , dlogf , sign(u.sub - x.sub)*Inf)
  return(rtn)
}

##Calculate probability function 
aug_probf = function(u , dt, v)
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
  u.sub <- u[1:(d-1)]
  x.sub <- x[1:(d-1)]
  s = matrix(u.sub - x.sub , nrow = length(u.sub) , ncol = 1)
  
  #Projection matrix
  den = t(s)%*%s
  den = den[1,1]
  
  foo = prod(s==0)
  #foo = prod(s < eps)
  count = 0
  if(foo>0)
  {
    I = diag(1,nrow =  d-1, ncol = d-1)
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
  x = matrix(0 ,nrow = N, ncol= d)
  
  #initializing the first row 
  x[1 , ] = x0
  ind = rep(FALSE, N)
  counter = 0
  for( i in 2: N)
  {
    v = rnorm(d-1)
    U = runif(d-1)
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
    x[i, -d] = x[(i-1), -d] + sqrt(2*dt)*mat%*%v
    x[i , d] = 1 - sum(x[i, -d]) #Using the dirichlet constraint
    #print(sum(x[i,]))
    if(i%% 1e3 ==0)
    {
      print(i)
    }
    
  } 
  return(list(x , ind, counter))
}

x0 = rep(1/d, d)
N <- 1e4
try = aug_UBA(x0, dt = .001 , N = N) 
c = try[[3]]
r = 1 - c/N
r
plot(try[[2]])

try1 = try[[1]]
out = try1[try[[2]] ,]

#plot(out , main = r)


# plot(density(out[,1]) , main = paste('X1, r = ',r))
# plot(density(out[,2]) , main = paste('X2, r = ',r) )
# print('r')
# r

foo.x <- seq(0,1, length = 1e3)
ind <- 2
plot(foo.x, dbeta(foo.x, shape1 = alpha[ind], shape2 = sum(alpha) - alpha[ind]), type = 'l')
lines(density(out[,ind]) , main = paste('X1, r = ',r), col = "red")

ind <- 100
plot(try[[1]][1:ind, ], col = try[[2]][1:ind] + 1, pch = 16, xlim = c(0,1), ylim = c(0,1))

plot(try[[1]][,1], col = try[[2]]+1, pch = 16)
