
####################
###Utilis functions
####################
## Euler-Maruyama Numerical Scheme for SDE 
## Calculate the scaling constant for pdf 
calcz = function(eps) {
  x = rexp(1e5)
  Y = x - (x^(2+eps))/(2+eps)
  avg = mean(exp(Y))
  z = 2 * avg
  return(z)
}
## Calculate the pdf 
pif = function(t , eps) {
  z = calcz(eps)
  num = abs(t)^(2+eps)
  denom = 2 + eps
  power = exp(-num/denom)
  return(power/z)
}

## Calculate drift function 
drift = function(t, eps){
  absdrift = -abs(t)^(1 + eps)
  dir = sign(t)
  return(dir*absdrift)
}
##Calculate the probability function for jumps 

probf = function(y ,sigma,eps, dt , v )
{
  mu = drift(y ,eps)
  temp = -2*sqrt(dt)*v*mu/sigma
  return(1/(1 + exp(temp)))
}

#################
###To functions
##################
##Euler Maruyama to get solution paths 
Euler_Maruyama = function(x0, eps, sigma,N ,dt)
{
  x = numeric(N)
  x[1] = x0
  sqrt_dt_sigma = sigma * sqrt(dt)
  for( i in 2:N)
  {
    mu = drift(x[i-1],eps)
    x[i] =  x[i-1]+ dt*mu+  sqrt_dt_sigma*rnorm(1)
  }
  return(x)
}

##Do the Unadjusted Barker Procedure 
UBarker = function(x0,sigma,eps, dt,N)
{
  x = numeric(N)
  x[1]= x0
  for( i in 2: N)
  {
    v = rnorm(1)
    U = runif(1)
    prob = probf(y = x[i-1], sigma = sigma , eps = eps , dt = dt , v = v)
    b = ifelse(U<=prob , 1, -1)
    x[i] = x[i-1] + sqrt(dt)*b*sigma*v
  }
  return(x)
}
###############################
###FUnction to Save the plots
###############################
save_comparison_plots <- function(eps = 4, dt = 0.0001, reps = 1, n = 1e5) {
  
  ## Initialise and run the simulations
  store_mat <- sapply(1:reps, function(i) {
    #print(i)
    result1 <- UBarker(x0 = 0, eps = eps, sigma = sqrt(2), dt = dt, N = n)
    result2 <- Euler_Maruyama(x0 = 0, eps = eps, sigma = sqrt(2), dt = dt, N = n)
    return(list(result1, result2))
  })
  
  store_Barker <- matrix(store_mat[[1]], nrow = reps, ncol = n, byrow = TRUE)
  store_EM =  matrix(store_mat[[2]], nrow = reps, ncol = n, byrow = TRUE)
  
  ## Generate the sequence and equilibrium distribution
  v <- seq(-3, 3, length.out = 100)
  truth <- pif(v, eps = eps)
  
  Barker = store_Barker[1,]
  EM = store_EM[1,]
  
  ## Create the plot
  plot(v, truth, col = 'black', main = paste('eps =', eps, "dt=", dt, "reps=", reps),
       type = "l", ylim = c(0, range(truth)[2] + 0.35))
  if(any(is.na(Barker))== FALSE)
  {
    lines(density(store_Barker[1,]), col = 'red')
  }
  if(any(is.na(EM))== FALSE)
  {
    lines(density(store_EM[1,]), col = 'blue')
  }
  
  legend('topright', legend = c('Truth',"UBA", "ULA"), 
         col = c('black','red','blue'), lty = 1)
}

####################
###Eps, dt, reps
####################
eps = 0:6
timesteps = c(2 ,1 , 0.1,0.01,0.001,0.0001, 0.00001)
reps = c(1, 1e3)
pdf('comparison_rep1.pdf', width = 24, height = 24)  # Increase size as needed


par(mfrow = c(7,7))
for(e in eps)
{
  for(dt in timesteps)
  { 
    print(paste(e,"_",dt))
    save_comparison_plots(eps = e, dt= dt, reps = 1, n = 1e5)
  }
}

dev.off()


