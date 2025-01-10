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
  
  z1 = x + z*(x - a)
  z2 = y + z*(y - b)
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
    x1=out[(i-1),1]
    x2= out[(i-1),2]
    print(c(x1,x2))
    
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
      m = (x2 - c0[2])/(x1 - c0[1])
      if(m == Inf || m == -Inf)
      {
        jump = c(0 , 1)
        #print('m is inf!')
      }else if(m == 0)
      {
        jump = c(1 ,0)
      }
      else{
        jump =  c(1/sqrt(1 + m^2), m/(1+ sqrt(m^2)))
      }
      
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

foo_AugUBA = AugUBA2d(x0 = c(-1,-1) , dt = 1e-3, N = 1e4) 
df = data.frame(round(foo_AugUBA,4))
plot(df$V1 , df$V2 , main = 'Chain Movement')
##Plotting arrow 
out =  matrix(0 , nrow = 1000 , ncol = 2)
a = 0 
b = 0
A = c(-1,-1) 
z = rnorm(1000)
for( i in 1: 1000)
{
  out[i, ] = arrow(A , c(a,b), z[i])
}
points(out[ , 1] , out[ ,2], col = 'green')
points(-1 ,-1 , col = 'blue' , pch = 16)
points(0,0, col = 'blue', pch = 16)
legend('bottomright' , legend = c('arrow' , 'chain') , 
       col = c('green' , 'black') , pch =1)



plot(df$V3 , df$V4 , main = 'Jumps' )
out =  matrix(0 , nrow = 1000 , ncol = 2)
a = 0 
b = 0
A = c(-1,-1) 
z = rnorm(1000)
for( i in 1: 1000)
{
  out[i, ] = arrow(A , c(a,b), z[i])
}
points(out[ , 1] , out[ ,2], col = 'green')
points(-1 ,-1 , col = 'blue' , pch = 16)
points(0,0, col = 'blue', pch = 16)
legend('bottomright' , legend = c('arrow' , 'jumps') , 
       col = c('green' , 'black') , pch =1)

dist = df$d 

plot( 1:150 , dist[1:150] , type = 'o' , main = 'Update Size')














#####################
# 
# plot(df$X3, df$X4 , type = 'l', col = 'red')
# plot(1:1e3,df$X3 , type = 'l')
# lines(1:1e3 , df$X1 , col = 'red')
# plot(1:1e3,df$X4 , type = 'l')
# lines(1:1e3 , df$X2 , col = 'blue')
# ts.plot(df$X4)
# abline(h = mean(df$X3), col = 'red')
# abline(h = 0 , col = 'grey')
# 
# df.clean = df[df$X1 > 0 & df$X2 > 0, ]
# correct= nrow(df.clean)
# ratio = 1 - (correct/1e3)
# ratio
# plot(df.clean$X1,df.clean$X2 , main = ratio)
# ts.plot(df.clean$X1)
# plot.ts(df.clean$X2)
#######################

##When in 3rd quadrant works fine 
##almost reaches the support
##but does not enter it.




