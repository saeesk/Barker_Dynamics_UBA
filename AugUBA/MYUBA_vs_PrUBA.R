#####################################
###UBA with Moreau Yosida Envelopes
#####################################
set.seed(123)
library(mvtnorm)
library(ggplot2)

lambda = 1e-8


#Calculate projetction
projection  = function(u)
{
  if(all(u[1] < 0, na.rm = TRUE) && all(u[2] < 0, na.rm = TRUE))
  {
    return(c(0,0))
  } else if (u[1]>0 & u[2]<0)
  {
    return(c(u[1],0))
  } else if (u[1]<0 & u[2]>0)
  {
    return(c(0,u[2]))
  }
  else{
    return(u)
  }
}

logf = function(x)
{
  psi.grad = x 
  u = projection(x)
  my.grad =  (x - u)/lambda
  rtn = -psi.grad - my.grad
  return(rtn)
} 

probf = function(x , dt,v)
{
  mu = logf(x)
  temp = -sqrt(2*dt)*v*mu
  return(1/(1 + exp(temp)))
}


UBA = function(x0 ,dt, N)
{
  d = length(x0)
  x = matrix(nrow = N , ncol = d)
  x[1, ]= x0
  for( i in 2: N)
  {
    v = rnorm(d)
    U = runif(d)
    prob = probf( x[i-1, ],dt = dt , v = v)
    b = ifelse(U<=prob , 1, -1)
    x[i, ] = x[i-1, ] + sqrt(2*dt)*b*v
  }
  return(x) 
  
}

############################
#### Preconditioned AugUBA
############################
drift = function(x)
{
  u = projection(x)
  s = u-x 
  rtn = ifelse( s ==0 , logf(x) ,sign(s)*Inf )
  return(rtn)
}

###Calculate Preconditioning matrix 
C = function(x)
{
  u = projection(x)
  s = matrix(u-x , nrow = 2, ncol = 1)
  den = t(s)%*%s
  den = den[1,1]
  if(s[1]== 0 & s[2]==0)
  {
    I = diag(1,nrow = 2 , ncol = 2)
    return(I)
  }
  else
  {
    P = s%*% t(s)
    P = P/den
    return(P)
  }
  
}
PrAugUBA = function(x0, dt, N)
{
  out = matrix(nrow = N , ncol = 2)
  out[1 ,] = x0
  
  for(t in 2:N)
  {
    x = out[t-1, ]
    v = matrix(rnorm(n = 2), ncol = 1)
    U = runif(2)
    
    ##Augmented Drift and Preconditioning matrix
    mu = matrix(drift(x), nrow = 2)
    mat = C(x)
    
    ###Calculating probabilities
    den = exp(-sqrt(2*dt)*v*mu)
    p = 1/(1+den)
    
    ###Making the updates
    z = ifelse(U <= p , v , -v)
    out[t , ] = out[t-1 , ] + sqrt(2*dt)*mat%*%z
    #print('Current states: ')
    #print(out[t, ])
  }
  return(out)
}


#######################
### Comparison plots
#######################

x0 = c(0.1,0.1)
N = 1e6
dt = 1e-4

###Truth
n_samples <- 10000
data <- rmvnorm(n = n_samples, mean = c(0,0))
data_df <- as.data.frame(data)
colnames(data_df) <- c("x", "y")
truncated_data <- subset(data_df, x > 0 & y > 0)  # Keep only points in the positive quadrant

####MYUBA
out = UBA(x0 =x0 , dt = dt , N = N)
MYUBA = as.data.frame(out)
MYUBA= MYUBA[ MYUBA$V1 > 0 & MYUBA$V2>0 , ]
r = 1 - dim(MYUBA)[1]/N
#plot(density(MYUBA$V1))
#plot(density(MYUBA$V2))

###AugUBA
try1 = PrAugUBA(x0 = x0, dt = dt , N = N)
df = as.data.frame(try1)
df.clean = df[df$V1 > 0 & df$V2 > 0, ]
correct= nrow(df.clean)
ratio = 1 - (correct/N)


###Scatter plot comparison 
plot(truncated_data$x, truncated_data$y,
     col = 'green', xlab = 'x',ylab = 'y',
     main = paste('MYUBA: ', round(r ,2) ,
           'AgUBA: ',  round(ratio,2)))
points(df.clean$V1, df.clean$V2, col = 'blue')
points(MYUBA$V1 , MYUBA$V2 , col = 'pink')
legend('topright', legend = c('Truth',"AgUBA", "MYUBA"), 
       col =c( 'green' , 'blue' , 'pink'), 
       pch = 16)



###Marginal densities
plot(density(truncated_data$x) , ylim = c(0,0.8),
     main = paste('MYUBA: ', round(r ,2) ,
                   'AgUBA: ',  round(ratio,2)))
lines(density(df.clean$V1) , col = 'blue')
lines(density(MYUBA$V1), col = 'green')
legend('topright', legend = c('Truth',"AgUBA", "MYUBA"), 
       col =c('black', 'blue' , 'green'), 
       lty= 1)
