###############################
#### Comparison of MLD and AugUBA
################################set.seed(123)
set.seed(123)
##Dimesion of the supprt
d = 11

#Parameters of Dirichlet Distribution
nl = rep(10,11)
#nl = c(3,3)
al= rep(0.1, d) 
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
##############################################################################
##############################################################################
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
  nl.sub <- nl[1:(d-1)]
  al.sub<- al[1:(d-1)]
  
  dlogf = (nl.sub + al.sub -1)/x.sub - (nl[d]+ al[d] -1)/(1 -sum(x.sub))
  
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
    if(i%% 1e4 ==0)
    {
      print(i)
    }
    
  } 
  return(list(x , ind, counter))
}

##########################################################################
##########################################################################

##Experiments for comparison
x0 = rep(1/d, d)
dt = 1e-5
N <- 1e6
out_MLD = MLD(x0, dt= dt,N=N)
try = aug_UBA(x0, dt = dt , N = N) 
c = try[[3]]
r = 1 - c/N
try1 = try[[1]]
out_AugUBA = try1[try[[2]] ,]


###Marginal density plots comparison
ind = 4
foo.x = seq(0,1, length.out = 1e3)
plot(foo.x , dbeta(foo.x , shape1 = pars[ind],
                   shape2 = sum(pars) - pars[ind]),
     type ='l', ylab = ' ', xlab = 'x',ylim = c(0,16), 
     main = paste('X_', ind, 'r = ',round(r,4),' dt= ', dt,' N= ',N))
lines(density(out_MLD[,ind]), col = 'red')
lines(density(out_AugUBA[,ind]), col = "blue")
legend('topright', legend = c('Truth', 'MLD', 'AugUBA'), 
       col = c('black', 'red', 'blue'), lty = 1)
