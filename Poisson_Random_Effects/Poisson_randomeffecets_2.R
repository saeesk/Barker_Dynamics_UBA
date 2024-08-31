##############################
#Poisson Random Effects Model
##############################
set.seed(123)
##Data Generation 
mustar = 5
I =50
J = 5
sigma_mu = 10
eta = rnorm(n = I , mean = mustar,sd = 1)
y = matrix(nrow = I , ncol = J)
for(i in 1:I)
{
  y[i ,] = rpois(n = J , lambda = exp(eta[i]))
}

########################################b##
## Functions for numerical schemes
######################################
####EM scheme : ULA

##Calculate drift 
pois_drift = function(x)
{
  mu = x[1]
  eta = x[-1]
  #derivative wrt mu
  dU1 = mu/100- sum(eta-mu)  
  #derivative wrt eta_i
  dU = J*exp(eta) - rowSums(y) + eta - mu
  #derivative of U 
  dU = c(dU1, dU)
  return(-1*dU)
}

##Unadjusted Langevin Algorithm 
ULA = function(x0,dt ,N)
{
  x = matrix(0 ,nrow = N,ncol= (I + 1))
  x[1 , ] = x0
  for( i in 2:N)
  {
    mu = pois_drift(x[(i-1) , ])
    x[i, ] =  x[(i-1), ]+ dt*mu+ sqrt(2*dt)*rnorm(1)
    #print(i)
  }
  
  return(x)
}

####Skew-Symmetric Scheme : UBA

###Calculate the probability function for jumps 
probf = function(y, dt , v )
{
  mu = pois_drift(y)
  temp = -sqrt(2*dt)*v*mu
  #print(v)
  return(1/(1 + exp(temp)))
}
# 
###Unadjusted Barker Algorithm 

UBA= function(x0,dt,N)
{
  x = matrix(0 ,nrow = N,ncol= (I + 1))
  x[1 , ] = x0
  for( i in 2: N)
  {
    v = rnorm(I+1)
    U = runif(I+1)
    prob = probf(y = x[(i-1), ], dt = dt , v = v)
    b = ifelse(U<= prob , 1, -1)
    x[i , ] = x[(i-1) , ] + sqrt(2*dt)*b*v
  }
  return(x)
}

#################
##Experiments
#################

##Setting 1 : Warm start
#Initialize
x0 = c(rnorm(1, mean = 5, sd = 10),eta)
foo = 1:1e4
timesteps =seq(from = 0.005 , to = 0.05 , by = 0.005)
ULA_mu = matrix(0 , nrow = 100 ,ncol = length(timesteps) )
UBA_mu = matrix(0 , nrow = 100 ,ncol = length(timesteps) )

for( j in 1 : length(timesteps))
{
  for( i in 1: 100)
  {
    r1= ULA(x0,  dt = timesteps[j], N = 1e5)
    r2 = UBA(x0,  dt = timesteps[j], N = 6e4)
    ULA_mu[i ,j] = mean(r1[,1])
    UBA_mu[i , j] = mean(r2[-foo,1])
    print(paste('Timestep: ', timesteps[j],'repetition no: ', i))
  }
  
  
  print(paste('Timestep number', i , 'done!'))
  
}

##Find column wise mean square error here! Using apply! 
ULA_mu = ULA_mu - mustar 
UBA_mu = UBA_mu - mustar 
err_ULA = colSums(ULA_mu^2)/nrow(ULA_mu)
err_UBA = colSums(UBA_mu^2)/nrow(UBA_mu)
##Make plots and save in pdf 
pdf('warm_start_version3.pdf')
plot(timesteps , err_mat[,1], type = 'o' , ylab = 'mean sqaured error', main = 'Warm Start', ylim = c(0.01,0.425))
lines(timesteps , err_mat[ ,2], col = 'red' , type = 'o')
legend('topright' , legend = c('ULA' , 'UBA') , col = c('black' , 'red') ,lty = 1 )
dev.off()
write.csv(err_mat, file = 'errors_warm.csv')

##########################################

###Setting 2 : Fixed start 
x0 = c(5,eta)
ULA_mu = matrix(0 , nrow = 100 ,ncol = length(timesteps) )
UBA_mu = matrix(0 , nrow = 100 ,ncol = length(timesteps) )

for( j in 1 : length(timesteps))
{
  for( i in 1: 100)
  {
    r1= ULA(x0,  dt = timesteps[j], N = 1e5)
    r2 = UBA(x0,  dt = timesteps[j], N = 6e4)
    ULA_mu[i ,j] = mean(r1[,1])
    UBA_mu[i , j] = mean(r2[-foo,1])
    print(paste('Timestep: ', timesteps[i],'repetition no: ', j))
  }
  

  print(paste('Timestep number', i , 'done!'))
  
}

##Find column wise mean square error here! Using apply! 
ULA_mu = ULA_mu - mustar 
UBA_mu = UBA_mu - mustar 
err_ULA = colSums(ULA_mu^2)/nrow(ULA_mu)
err_UBA = colSums(UBA_mu^2)/nrow(UBA_mu)

#Make plot and save it 
pdf('fixed_start_version3.pdf')
plot(timesteps , err_ULA, type = 'o' , ylab = 'mean sqaured error', main = 'Fixed Start')
lines(timesteps , err_UBA, col = 'red' , type = 'o')
legend('topright' , legend = c('ULA' , 'UBA') , col = c('black' , 'red') ,lty = 1 )
dev.off()
write.csv(err_mat, file = 'errors_fixed.csv')
#




