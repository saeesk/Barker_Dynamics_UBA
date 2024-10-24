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
# pdf('MLDvsUBA2.pdf', width = 10 , height = 10)
timesteps = c(1e-2,1e-3, 1e-4, 1e-6)
steps = c(1e5) 
lambda = 10
par(mfrow = c(4,3))


for( dt in timesteps) {
   for ( N in steps)
   {
    ##Logs
    print(paste('Timestep : ' , dt , ', N : ',N))
   
    ##AugUBA 
    r1= Aug_UBA(x0 = 0.1 , dt ,N )
    ratio = r1[1]
    r1 = r1[-1]
    
    ##MLD
    r2 = MLD(x0 = 0.1, dt , N)
    
    ##Plots
    vec = seq(0,max(r1),by = 0.01)
    plot(density(r1),
         col = 'red' , ylim = c(0 , lambda+2),
         main = paste('dt= ' , dt , ' N= ', N ,'r =' , ratio))
    lines(vec , dexp(vec,rate = lambda))
    lines(density(r2), col = 'blue')
    legend('topright' , legend = c('AugUBA','MLD', 'Truth') , col = c('red' ,'blue', 'black'), lty = 1)

    acf(r2, col = "blue", lag.max = 100)
    lines(acf(r1, plot = FALSE, lag.max = 100)$acf, col = "red")

    plot.ts(r1, col = "red")
    lines(r2, col = "blue")
  }
}
#dev.off()


