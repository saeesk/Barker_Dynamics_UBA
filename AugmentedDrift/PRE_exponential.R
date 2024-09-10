#########################
##PRE positive parameter
#########################
set.seed(0809204)
#Data Generation 
lambdastar = 1e-4
I =50
J = 5
sigma_l = 1e-4
eta = rnorm(n = I , mean = lambdastar , sd = 1)
y = matrix(nrow = I , ncol = J)
for(i in 1:I)
{
  y[i ,] = rpois(n = J , lambda = exp(eta[i]))
}

##AgUBA 
aug_drift = function(x)
{
  #Augmented gradient for lambda
  lambda = x[1]
  dU1 = ifelse(lambda>0 , -sum(eta - lambda)+ sigma_l , -Inf)

  #Gradient for eta 
  eta = x[-1]
  dU2 = J*exp(eta) - rowSums(y)
  
  dU = c(dU1,dU2)
  return(-dU)
}

aug_probf = function(u, dt , v )
{
  mu = aug_drift(u)
  temp = -sqrt(2*dt)*v*mu
  return(1/(1 + exp(temp)))
}

aug_UBA= function(x0,dt,N)
{
  x = matrix(0 ,nrow = N,ncol= (I + 1))
  #initializing the first row 
  x[1 , ] = x0
  
  for( i in 2: N)
  {
    v = rnorm(I+1)
    U = runif(I+1)
    prob = aug_probf(u = x[(i-1), ], dt = dt , v = v)
    b = ifelse(U<= prob , 1, -1)
    x[i , ] = x[(i-1) , ] + sqrt(2*dt)*b*v
  } 
  pst_lambda= x[ ,1]
  pst_lambda =  pst_lambda[pst_lambda >0] 
  ratio = 1 - length(pst_lambda)/N
  return(c(ratio,pst_lambda))
}

# ##Tetsing if code is correct
# l_aug = aug_UBA(x0, 1e-5,5e5)
# r = l_aug[1]
# foo = 1: (1e4+1)
# l_aug = l_aug[-foo]
# plot(density(l_aug) , main = paste("Aug UBA , r = " ,r) )
# abline(v = lambdastar , col = 'blue', lty = 2)


##UBA 
drift = function(x)
{
  # gradient for lambda
  lambda = x[1]
  dU1 = sigma_l -sum(eta - lambda)
  
  #Gradient for eta 
  eta = x[-1]
  dU2 = J*exp(eta) - rowSums(y)
  
  dU = c(dU1,dU2)
  return(-dU)
}

probf = function(u , dt , v)
{
  mu = drift(u)
  temp = -sqrt(2*dt)*v*mu
  return(1/(1 + exp(temp)))
}

UBA = function(x0,dt,N)
{
  x = matrix(0 ,nrow = N,ncol= (I + 1))
  #initializing the first row 
  x[1 , ] = x0
  
  for( i in 2: N)
  {
    v = rnorm(I+1)
    U = runif(I+1)
    prob = probf(u = x[(i-1), ], dt = dt , v = v)
    b = ifelse(U<= prob , 1, -1)
    x[i , ] = x[(i-1) , ] + sqrt(2*dt)*b*v
  } 
  pst_lambda= x[ ,1]
  return(pst_lambda)
}
# 
# ##Tetsing if code is correct
# l_uba = UBA(x0, 0.005,6e4)
# ps = l_uba[l_uba<0]
# print(length(ps)/6e4)
# 
# foo = 1:1e4
# l_ula = l_aug[-foo]
# plot(density(l_ula) , main = " UBA")
# abline(v = lambdastar , col = 'red', lty = 2)
# 
# lambdastar
# median(l_aug)
# median(l_uba)

####################################
######################################
##Experiments####
##Setting 1 : Warm start
#Initialize
x0 = c(rexp(1, rate =1/ lambdastar),eta)
foo = 1:1e4
timesteps =seq(from = 0.005 , to = 0.05 , by = 0.005)
UBA_l = matrix(0 , nrow = 100 ,ncol = length(timesteps))
AgUBA_l = matrix(0 , nrow = 100,ncol = length(timesteps)) 
ratio_mat1 =  matrix(0 , nrow = 100,ncol = length(timesteps))
ratio_mat2 =  matrix(0 , nrow = 100,ncol = length(timesteps))

for( j in 1: length(timesteps))
{
  for( i in 1:100) 
  {
    r1 = UBA(x0,timesteps[j],6e4) 
    ratio_mat1[i,j]= length(r1[r1>0])/6e4
    
    r2 = aug_UBA(x0,timesteps[j],6e4) 
    ratio_mat[ i, j] = r2[1]  
    
    UBA_l[ i , j] = mean(r1[-foo])
    AgUBA_l[i,j] = mean(r2[-foo])
    
    print(paste('Timestep: ', timesteps[j],'repetition no: ', i))
  }
  print(paste('Timestep number', j , 'done!'))
} 

ratio_vec1 = colMeans(ratio_mat1) 
ratio_vec2 = colMeans(ratio_mat2)

UBA_l = UBA_l - lambdastar 
AgUBA_l = AgUBA_l - lambdastar  
err_UBA = colMeans(UBA_l^2)
err_AgUBA = colMeans(AgUBA_l^2) 
err_mat = cbind(err_UBA, err_AgUBA) 


write.csv(cbind(err_mat,ratio_vec1 , ratio_vec2) , 'errors_warmstart.csv')

###############
##Setting 2 : Fixed start 
x0 = c(lambdastar,eta) 
foo = 1:1e4
timesteps =seq(from = 0.005 , to = 0.05 , by = 0.005)
UBA_l = matrix(0 , nrow = 100 ,ncol = length(timesteps))
AgUBA_l = matrix(0 , nrow = 100,ncol = length(timesteps)) 
ratio_mat1 =  matrix(0 , nrow = 100,ncol = length(timesteps))
ratio_mat2 =  matrix(0 , nrow = 100,ncol = length(timesteps))
for( j in 1: length(timesteps))
{
  for( i in 1:100) 
  {
    r1 = UBA(x0,timesteps[j],6e4)
    ratio_mat1[i,j] = length(r1[r1<0])/6e4
    
    r2 = aug_UBA(x0,timesteps[j],6e4) 
    ratio_mat2[ i, j] = r2[1]  
    
    UBA_l[ i , j] = mean(r1[-foo])
    AgUBA_l[i,j] = mean(r2[-foo])
    
    print(paste('Timestep: ', timesteps[j],'repetition no: ', i))
  }
  print(paste('Timestep number', j , 'done!'))
} 

ratio_vec1 = colMeans(ratio_mat1) 
ratio_vec2 = colMeans(ratio_mat2)

UBA_l = UBA_l - lambdastar 
AgUBA_l = AgUBA_l - lambdastar  
err_UBA = colMeans(UBA_l^2)
err_AgUBA = colMeans(AgUBA_l^2) 
err_mat = cbind(err_UBA, err_AgUBA) 

write.csv(cbind(err_mat ,ratio_vec1,ratio_vec2), 'errors_fixedstart.csv')  
