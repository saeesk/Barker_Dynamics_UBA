set.seed(123)
#############################
###Mirror Langevin Dynamics
#############################


###Mirror Langevin Dynamics
MLD = function(x0,dt,N) 
{
  x = numeric(N)
  y = numeric(N)
  x[1] = x0
  for(i in 1:(N-1))
  {
    nu = rnorm(1)
    y[i]= 1 + log(x[i])
    y[i+1] = y[i] - lambda*dt*x[i]+dt + sqrt(2*dt)*nu
    x[i+1] = exp(y[i+1]-1)
  }
  return(x)
}
####################
## MLD h(x)= log x
####################
MLD1 = function(x0,dt,N) 
{
  x = numeric(N)
  y = numeric(N)
  x[1] = x0
  for(i in 1:(N-1))
  {
    nu = rnorm(1)
    y[i]= 1/x[i]
    y[i+1] = y[i] + dt*(lambda - 2/x[i])*(x[i])^2 + sqrt(2*dt)*nu
    x[i+1] = y[i+1]
  }
  return(x)
}

f1 = MLD1(x0 = 0.1, dt = 1e-3, N = 1e3)
x = seq(0,2,length.out = length(f1))
plot(density(f1), col = 'green',
     main ='dt = 1e-3, N = 1e3', xlim = c(0,3))
lines(x,dexp(x, rate = 10))

f2 = MLD1(x0 = 0.1, dt = 1e-3, N = 1e5)
x = seq(0,2,length.out = length(f2))
plot(density(f2), col = 'green',
     main ='dt = 1e-3, N = 1e5', xlim = c(0,3))
lines(x,dexp(x, rate = 10))


f3 = MLD1(x0 = 0.1, dt = 1e-4, N = 5e6)
x = seq(0,2,length.out = length(f3))
plot(density(f3), col = 'green',
     main ='dt = 1e-4, N = 5e6', 
     xlim = c(0,3), ylim = c(0,lambda))
lines(x,dexp(x, rate = 10))

##############################
###UBA with Augmented drift 
##############################
aug_drift = function(x)
{ 
  if(x>0)
  {
    return(-lambda)
  }
  else
    return(Inf)
}


###Calculate the probability function for jumps 
probf = function(y, dt , v )
{
  mu = aug_drift(y)
  temp = -sqrt(2*dt)*v*mu
  #print(v)
  return(1/(1 + exp(temp)))
}
# 
###Unadjusted Barker Algorithm 

Aug_UBA= function(x0,dt,N)
{
  x = numeric(N)
  x[1]= x0
  for( i in 2: N)
  {
    v = rnorm(1)
    U = runif(1)
    prob = probf(y = x[i-1], dt = dt , v = v)
    b = ifelse(U<=prob , 1, -1)
    x[i] = x[i-1] + sqrt(2*dt)*b*v
  }
  
  #remove negative values
  lost_values = length(x[x<0]) 
  ratio = lost_values/N
  x = x[x>0]
  rtn = c(ratio , x)
  return(rtn)
}

###################
### Experiments
###################
pdf('MLDvsUBA3.pdf', width = 10 , height = 10)
timesteps = c(1e-3,1e-4,1e-6)
steps = c(1e4,1e5,5e6) 
lambda = 10
par(mfrow = c(3,3))

for( dt in timesteps) {
   for ( N in steps)
   {
    ##Logs
    print(paste('Timestep : ' , dt , ', N : ',N))
   
    ##AugUBA 
    r1= Aug_UBA(x0 = 0.1 , dt ,N )
    ratio = r1[1]
    r1 = r1[-1]
    
    ##MLD - entropic map
    r2 = MLD(x0 = 0.1, dt , N)
    
    ##MLD - log map 
    r3 = MLD1(x0 = 0.1, dt, N)
    
    ##Plots
    vec = seq(0,max(r1),by = 0.01)
    plot(density(r1),
         col = 'red' , ylim = c(0 , lambda+2),
         main = paste('dt= ' , dt , ' N= ', N ,'r =' , ratio))
    lines(vec , dexp(vec,rate = lambda))
    lines(density(r2), col = 'blue')
    lines(density(r3) , col = 'green')
    legend('topright' , legend = c('AugUBA','MLD_entr','MLD_log','Truth') , col = c('red' ,'blue','green', 'black'), lty = 1)

    
  }
}
dev.off()



