#####################################
### Bivariate Truncated Distribution
####################################

###################################
#### Function to find projection
###################################
projection  = function(u)
{
  if(all(u[1] < 0, na.rm = TRUE) && all(u[2] < 0, na.rm = TRUE))
  {
    return(c(0,0))
  }
  else if (u[1]>0 & u[2]<0)
  {
    return(c(u[1],0))
  }
  else
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
  
  z1 = x + z*(x - a)
  z2 = y + z*(b - y)
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
  
  if( x1>0 & x2>0)
  {
    return(-x)
  }
  else
  {
    return(c(Inf, Inf))
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
  out[1, ] = x0
  
  for( i in 2:N)
  {
    x1=out[(i-1),1]
    x2= out[(i-1),2]
    
    U = runif(2)
    
    ##If point is inside the support, run UBA
    if (all(x1 > 0, na.rm = TRUE) && all(x2 > 0, na.rm = TRUE))
    {
      v = rnorm(2)
      prob = p(out[i, ], dt = dt , v= v ) 
      b = ifelse(U<= prob , 1, -1)
    }
    else ##Point is outside the support, run AugUBA
    { 
      z = rnorm(1)
      prob = p(out[(i-1), ], dt = dt , v= rep(z,2) )
       
      c0 = projection(out[(i-1),]) #find projection
      v = arrow(out[(i-1),], c0,z) #shoot the arrow
      
      # #Calculate the unit sized jump along the arrow
      m = (c0[2] - x2)/(x1 - c0[1])
      jump =  c(1/sqrt(1 + m^2), m/(1+ sqrt(m^2)))
      #ahead = out[(i-1) , ]+ jump
      #back =  out[(i-1) , ]- jump
      
      b = ifelse(U<= prob , jump, -jump)
    }
    out[i , ] = out[(i-1) , ] + sqrt(2*dt)*b*v
  }
  return(out)
}

### Example 
foo =AugUBA2d(c(-0.001,-0.001) , dt = 1e-3 , N = 1e6)
plot(foo[ ,2] , foo[ ,1])
