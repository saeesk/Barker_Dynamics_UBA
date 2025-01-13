##############
###AugUBApD 
##############
set.seed(123)

###Find projection: Need to change, current verstion for 2D
projection = function(u)
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

###Shooting the arrow
arrow = function(x ,z, k)
{
  prj = projection(x)
  u = matrix(c(z , 1), ncol = 1 , byrow = TRUE)
  sq_x = sqrt(abs(x))
  A = cbind( x - prj, sq_x)
  v = A%*% u
  return(v)
}

z = rnorm(100)
vec = matrix(0 , nrow = 100 , ncol =2)
pt = c(-1,1)
for( i in 1:100)
{
  vec[i, ] = arrow(pt, z[i], k = 2)
}
plot(vec[ , 1], vec[,2])

##Augmented Drift: For truncated gaussian 
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

##Probability function 
p = function(x,dt,v)
{
  mu = aug_drift(x) 
  temp = sqrt(pi*dt)*v*mu
  return(pnorm(temp/2)) 
}

#p(c(2,-3), 1e-2,z)

AugUBA_k= function(k ,x0, dt, N)
{
  #Initialize
  out = matrix(0 , nrow = N , ncol = k)
  step = matrix(0 , nrow= N , ncol = k)
  d = numeric(N)
  out[1, ] = x0
  step[1, ] = 0
  d[1] = 0
  for( i in 2: N)
  {
  #Update 
  U = runif(k)
  ##If point is inside the support, run UBA
  if(!is.na(x1) && x1 > 0 && !is.na(x2) && x2 > 0){
    print('Entering UBA')
    v = rnorm(n = k , mean = 0, sd = 1)
    prob = p(out[(i-1), ], dt = dt, v = v)
    b = ifelse(U <= prob, 1, -1)
  }else{ 
    ##Point is outside the support, run AugUBA
    print('Entering AugUBA')
    z = rnorm(n=1 , mean = 0 , sd =1)
    prob = p(out[(i-1), ], dt = dt , v= rep(z,k) )
    
    c0 = projection(out[(i-1), ]) #find projection
    v = arrow(x = out[(i-1), ],z = z, k = 2) #shoot the arrow
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
foo_AugUBA = AugUBA_k(k = 2,x0 = c(-1,-1) , dt = 1e-3, N = 1e5) 
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

