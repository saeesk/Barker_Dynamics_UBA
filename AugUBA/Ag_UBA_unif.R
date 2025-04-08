##############################
###Example : U(0,1)
##############################
##############################
###UBA with Augmented drift 
##############################
aug_drift = function(x)
{ 
  if(x>1)
  {
    return(-Inf)
  }
  else if(x<0)
  {
    return(Inf)
  }
  else
  {
    return(0)
  }
    
}

p = function(y,dt,v)
{
  mu = aug_drift(y)
  temp = sqrt(pi*dt)*v*mu
  return(pnorm(temp/2)) 
}

Aug_UBA= function(x0,dt,N)
{
  x = numeric(N)
  x[1]= x0
  for( i in 2: N)
  {
    v = rnorm(1)
    U = runif(1)
    prob = p(y = x[i-1], dt = dt , v = v)
    b = ifelse(U<=prob , 1, -1)
    x[i] = x[i-1] + sqrt(2*dt)*b*v
  }
  
  #remove negative values
  lost_values = length(x[x<0 | x>1]) 
  ratio = lost_values/N
  x = x[x>0 & x<1]
  rtn = c(ratio , x)
  return(rtn)
}


r1= Aug_UBA(0.1,0.1,1e6)
ratio = r1[1]
r1 = r1[-1]

vec = seq(0, 1,by = 0.01)
plot(density(r1), col = 'red' ,main = ratio, ylim = c(0,1.2))
lines(vec , dunif(vec, min = 0 , max = 1))
legend('topright' , legend = c('AugUBA', 'Truth') , col = c('red' , 'black'), lty = 1)

pdf('Ag_Drift_Unif_Probit.pdf')
timesteps = c(1,1e-3, 1e-6)
steps = c(1e4, 5e5,1e7) 
par(mfrow = c(3 ,3), mar = c(4, 4, 2, 2), oma = c(1, 1, 1, 4))

sapply(timesteps , function(dt){
  sapply(steps , function(N){
    print(paste('Timestep : ' , dt , ', N : ',N))
    r1= Aug_UBA(x0 = 0.1 , dt ,N )
    ratio = r1[1]
    r1 = r1[-1]
    vec = seq(0, 1,by = 0.01)
    plot(density(r1), 
         col = 'red' , ylim = c(0 , 1.2),
         main = paste('r =' , round(ratio,2), 'dt=', dt, 'N= ', N))
    lines(vec , dunif(vec, min = 0 , max = 1))
    
    
    
  })
})
#par(xpd = NA)
#legend('topright' ,inset = c(-0.5, -2), legend = c('AugUBA', 'Truth') , col = c('red' , 'black'), lty = 1, bty = "n")
dev.off()



