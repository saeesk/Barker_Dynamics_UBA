#####################################
### Bivariate Truncated Distribution
####################################
###################################
#### Function to find projection
###################################
set.seed(123) 
projection  = function(u)
{
  if(all(u[1] < 0, na.rm = TRUE) && all(u[2] < 0, na.rm = TRUE))
  {
    return(c(0,0))
  } else if (u[1]>0 & u[2]<0)
  {
    return(c(u[1],0))
  } else
  {
    return(c(0,u[2]))
  }
}

####################################
##### Shooting an arrow function
####################################
#u = (x,y): is current position
#c0 = (a,b): is projection from current position 
# v : random draw from standard gaussian distribution
arrow = function(u,c0,z)
{
  x = u[1]
  y = u[2]
  a = c0[1]
  b = c0[2]
  
  z1 = sqrt(abs(x)) + z*(x - a)
  z2 = sqrt(abs(y)) + z*(y - b)
  return(c(z1,z2))
}

# z = rnorm(1)
# a = 0 
# b = 0
# arrow(4,-5,z)
####################
### Augmented Drift
#################### 
aug_drift= function(x)
{
  x1 = x[1]
  x2 = x[2]
  
  if( x1>0 && x2>0)
  {
    return(-x)
  }
  else
  {
    return(c(-Inf, -Inf))
  }
}

#aug_drift(c(a,b))

########################
##Probability function
#######################
p = function(x,dt,v)
{
  mu = aug_drift(x) 
  temp = sqrt(pi*dt)*v*mu
  return(pnorm(temp/2)) 
}

#p(c(2,-3), 1e-2,z)

###############################################
### AugUBA2D : Truncated Bivariate Gaussian
################################################
AugUBA2d = function(x0,dt,N)
{
  #Initialize
  out = matrix(0 , nrow = N , ncol = 2)
  step = matrix(0 , nrow= N , ncol = 2)
  d = numeric(N)
  out[1, ] = x0
  step[1, ] = 0
  d[1] = 0
  for( i in 2:N)
  {
    x1 = out[(i-1),1]
    x2 = out[(i-1),2]
    U = runif(2)
    
    ##If point is inside the support, run UBA
    if(!is.na(x1) && x1 > 0 && !is.na(x2) && x2 > 0){
      print('Entering UBA')
      v = rnorm(n = 2 , mean = 0, sd = 1)
      prob = p(out[(i-1), ], dt = dt, v = v)
      b = ifelse(U <= prob, 1, -1)
    }else{ 
      ##Point is outside the support, run AugUBA
      print('Entering AugUBA')
      z = rnorm(n=1 , mean = 0 , sd =1)
      prob = p(out[(i-1), ], dt = dt , v= rep(z,2) )
      
      c0 = projection(out[(i-1), ]) #find projection
      v = arrow(out[(i-1), ], c0,z) #shoot the arrow
      # #Calculate the unit sized jump along the arrow
      dist = sqrt(sum(out[(i-1) , ] - c0)^2)
      jump =  (c0 - out[(i-1), ])/dist
      b = ifelse(U<= prob , jump, -jump)
    }
    
    out[i , ] = out[(i-1) , ] + sqrt(2*dt)*b*v
    step[i, ] = sqrt(2*dt)*b*v
    d[i]= 2*dt*sum( (b*v)^2)
    print(i)
  } 
  
  rtn = cbind(out, step , d)
  return(rtn)
}

### Example 

foo_AugUBA = AugUBA2d(x0 = c(1,-1) , dt = 1e-3, N = 1e5) 
df = data.frame(round(foo_AugUBA,4))

plot(df$V1 , df$V2 , main = 'Chain Movement')
df.clean = df[df$V1 > 0 & df$V2 > 0, ]
correct= nrow(df.clean)
ratio = 1 - (correct/1e5)
ratio
colnames(df.clean)
ts.plot(df.clean$V1)
plot.ts(df.clean$V2)

plot(df.clean$V1,df.clean$V2 , main = paste("Simulation",ratio))
##### Truth 
# Load required libraries
library(MASS)     # For mvrnorm to generate bivariate normal data
library(ggplot2)  # For visualization

# Set parameters for the standard bivariate Gaussian
mean <- c(0, 0)               # Mean vector (standard normal)
sigma <- matrix(c(1, 0,       # Covariance matrix (identity for standard case)
                  0, 1), 
                nrow = 2)


# Generate random samples
# For reproducibility
n_samples <- 100000
data <- mvrnorm(n = n_samples, mu = mean, Sigma = sigma)

# Convert to a data frame and apply truncation
data_df <- as.data.frame(data)
colnames(data_df) <- c("x", "y")
truncated_data <- subset(data_df, x > 0 & y > 0)  # Keep only points in the positive quadrant

# Plot the truncated data
ggplot(truncated_data, aes(x = x, y = y)) +
  geom_point(alpha = 0.3, color = "blue") +       # Scatterplot
  labs(title = "Truncated Standard Bivariate Gaussian (Positive Quadrant)",
       x = "X",
       y = "Y") +
  theme_minimal()


