library(mvtnorm)
library(ggplot2)
###Augmented drift
drift = function(x)
{
  x1 = x[1]
  x2 = x[2]
  
  if( x1>0 && x2>0)
  {
    return(-x)
  }
  else if(x1 <0 && x2<0)
  {
    return(c(Inf, Inf))
  }
  else if(x1 <0 && x2>0 )
  {
    return(c(Inf,0))
  }
  else
  {
    return(c(0,Inf))
  }
}

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

# Generate random samples
# For reproducibility
n_samples <- 10000
data <- rmvnorm(n = n_samples, mean = c(0,0))

# Convert to a data frame and apply truncation
data_df <- as.data.frame(data)
colnames(data_df) <- c("x", "y")
truncated_data <- subset(data_df, x > 0 & y > 0)  # Keep only points in the positive quadrant
pdf(file="PreAugUBA.pdf")
par(mfrow = c(2,2))

try1 = PrAugUBA(x0 = c(-1,1), dt = 1e-3 , N = 1e4)
df = as.data.frame(try1)
df.clean = df[df$V1 > 0 & df$V2 > 0, ]
correct= nrow(df.clean)
ratio = 1 - (correct/1e4)
plot(truncated_data$x, truncated_data$y, col = 'green', xlab = 'x',ylab = 'y',
     main = paste( 'x0 = c(-1,1), ratio= ',ratio ))
points(df.clean$V1, df.clean$V2, col = 'blue')
legend('topright', legend = c('Truth',"PreAugUBA"), col =c('green', 'blue'), pch= 16)

try2 = PrAugUBA(x0 = c(-1,-1), dt = 1e-3 , N = 1e4)
df = as.data.frame(try2)
df.clean = df[df$V1 > 0 & df$V2 > 0, ]
correct= nrow(df.clean)
ratio = 1 - (correct/1e4)
plot(truncated_data$x, truncated_data$y, col = 'green', xlab = 'x',ylab = 'y',
     main = paste( 'x0 = c(-1,-1), ratio= ',ratio ))
points(df.clean$V1, df.clean$V2, col = 'blue')
legend('topright', legend = c('Truth',"PreAugUBA"), col =c('green', 'blue'), pch= 16)


try3 = PrAugUBA(x0 = c(1,-1), dt = 1e-3 , N = 1e4)
df = as.data.frame(try3)
df.clean = df[df$V1 > 0 & df$V2 > 0, ]
correct= nrow(df.clean)
ratio = 1 - (correct/1e4)
plot(truncated_data$x, truncated_data$y, col = 'green', xlab = 'x',ylab = 'y',
     main = paste( 'x0 = c(1,-1), ratio= ',ratio ))
points(df.clean$V1, df.clean$V2, col = 'blue')
legend('topright', legend = c('Truth',"PreAugUBA"), col =c('green', 'blue'), pch= 16)


try4 = PrAugUBA(x0 = c(1,1), dt = 1e-3 , N = 1e4)
df = as.data.frame(try1)
df.clean = df[df$V1 > 0 & df$V2 > 0, ]
correct= nrow(df.clean)
ratio = 1 - (correct/1e4)
plot(truncated_data$x, truncated_data$y, col = 'green', xlab = 'x',ylab = 'y',
     main = paste( 'x0 = c(1,1), ratio= ',ratio ))
points(df.clean$V1, df.clean$V2, col = 'blue')
legend('topright', legend = c('Truth',"PreAugUBA"), col =c('green', 'blue'), pch= 16)
dev.off()

