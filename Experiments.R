
####################
###Utilis functions
####################
## Euler-Maruyama Numerical Scheme for SDE 
## Calculate the scaling constant for pdf 
calcz = function(eps) {
  X = rnorm(1e5)
  Y = (abs(X)^eps)/(1+eps/2)
  avg = mean(Y) 
  adj = sqrt(2*pi)
  z = adj * avg
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
  x = numeric(n)
  x[1] = x0
  sqrt_dt_sigma = sigma * sqrt(dt)
  for( i in 2:n)
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

###################
#### Experiments 
###################
##Initialise
eps = 4
dt = 0.0001
reps = 1
n = 1e5
##Get Results for both schemes
store_mat <- sapply(1:reps, function(i) {
    result1 <- UBarker(x0 = 0, eps = eps, sigma = sqrt(2), dt =dt, N = n)
    result2 <- Euler_Maruyama(x0 = 0, eps = 2, sigma = sqrt(2), dt = 0.01, N = n)
    return(list(result1,result2))
  })
store_Barker <- matrix(store_mat[[1]], nrow = 1000, ncol = n, byrow = TRUE)
store_EM =  matrix(store_mat[[2]], nrow = 1000, ncol = n, byrow = TRUE)

#Comparison Plots  
v = seq(-3,3,length.out = 100)
#Equilibrium distribution 
truth = pif(v , eps = eps)
plot(v ,truth, col = 'black' , main = paste('eps = ' ,eps," dt= ", dt, " reps= ",reps),
       type = "l" , ylim = c(0 ,range(truth)[2]+ 0.35))
lines(density(store_Barker[ ,]) , col = 'red') 
lines(density(store_EM[ , ]), col = 'blue')

  
###############################
###FUnction to Save the plots
###############################
save_comparison_plots <- function(eps = 4, dt = 0.0001, reps = 1, n = 1e5) {
  ## Initialise and run the simulations
  store_mat <- sapply(1:reps, function(i) {
    print(i)
    result1 <- UBarker(x0 = 0, eps = eps, sigma = sqrt(2), dt = dt, N = n)
    result2 <- Euler_Maruyama(x0 = 0, eps = 2, sigma = sqrt(2), dt = 0.01, N = n)
    return(list(result1, result2))
  })
  store_Barker <- matrix(store_mat[[1]], nrow = reps, ncol = n, byrow = TRUE)
  store_EM =  matrix(store_mat[[2]], nrow = reps, ncol = n, byrow = TRUE)
  ## Generate the sequence and equilibrium distribution
  v <- seq(-3, 3, length.out = 100)
  truth <- pif(v, eps = eps)
  
  ## Open a PNG graphics device
  filename = paste("eps",eps,"dt",dt,"reps",reps,".png" , sep = "")
  png(filename = filename)
  
  ## Create the plot
  plot(v, truth, col = 'black', main = paste('eps =', eps, "dt=", dt, "reps=", reps),
       type = "l", ylim = c(0, range(truth)[2] + 0.35))
  lines(density(store_Barker[,-1]), col = 'red')
  lines(density(store_EM[, -1 ]), col = 'blue')
  
  ## Close the graphics device
  dev.off()
  
  ## Print a message to indicate that the plot has been saved
  cat("Plot saved as", filename, "\n")
}

# ## Example usage of the function
# save_comparison_plots(eps = 5, dt = 0.1, reps = 1, n = 1e5)
# save_comparison_plots(eps = 5, dt = 0.01, reps = 1, n = 1e5)
# save_comparison_plots(eps = 5, dt = 0.001, reps = 1, n = 1e5)
# save_comparison_plots(eps = 5, dt = 0.0001, reps = 1, n = 1e5)
####################
###Eps, dt, reps
####################
eps = 3:5
timesteps = c(0.1,0.01,0.001,0.0001)
reps = c(1, 1e3)
for(e in eps)
{
  for(dt in timesteps)
  {
    save_comparison_plots(eps = e, dt= dt, reps = 500, n = 1e5)
  }
}


