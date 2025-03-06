
############################################
##This code is AugUBA sampler for truncated 
##bivariate gaussian on circular suuport
## x^2 + y^2 <= 1
############################################
set.seed(123)
## Density
logf = function(x)
{
  if(sum(x^2)<1)
  {
    rtn = -log(2*pi)- log(1 - exp(-0.5))- sum(x^2)/2
  }
  else
  {
    rtn = -Inf
  }
  return(rtn)
}

#Calculate projection
projection = function(x)
{
  rtn = min( 1, 1/sum(x^2))*x
  return(rtn)
}


##Calculate Augmented Drift
aug_drift = function(x)
{
  u = projection(x) 
  drift = ifelse(u - x == 0, -x, sign(u-x)*Inf)
  return(drift)
}


## Calculate the probabiity fucntion
aug_probf = function(u, dt , v )
{
  mu = aug_drift(u)
  temp = -sqrt(2*dt)*v*mu
  return(1/(1 + exp(temp)))
}

###Calculate preconditioning matrix
C = function(x)
{
  u = projection(x)
  s = matrix(u-x , nrow =2 , ncol = 1)
  
  #Projection matrix
  den = t(s)%*%s
  den = den[1,1]
  
  foo = sum(s==0)
  if(foo>0)
  {
    I = diag(1,nrow = 2 , ncol = 2)
    rtn = I
  }else
  {
    P = s%*% t(s)
    P = P/den
    rtn = P 
  }
  return(rtn)
}

aug_UBA = function(x0,dt,N)
{
  
  x = matrix(0 ,nrow = N,ncol=2)
  #initializing the first row 
  x[1 , ] = x0
  
  for( i in 2: N)
  {
    v = rnorm(2)
    U = runif(2)
    prob = aug_probf(u = x[(i-1), ], dt = dt , v = v)
    v = ifelse(U<= prob , v, -v)
    mat = C(x[(i-1), ])
    x[i , ] = x[(i-1) , ] + sqrt(2*dt)*mat%*%v
    if(i%% 1e2 ==0)
    {
      print(i)
    }
    
  } 
  return(x)
}

### ACcept reject to compare
AR = function(N) {
  samples = matrix(NA, ncol = 2, nrow = N)  # Storage for accepted samples
  i = 1  # Counter for accepted samples
  count = 0
  while (i <= N) {
    y = rnorm(2)  # Sample from standard bivariate normal (N(0,1))
    U = runif(1)  # Uniform random number in [0,1]
    
    # Compute acceptance criterion
    if (sum(y^2) < 1 ) {
      samples[i, ] = y  # Accept the sample
      i = i + 1  # Move to next sample
      count = count + 1
      #print(count)
    }
  }
  
  return(list(count/N,samples))  # Return the accepted samples
}

samples = AR(1e4)
samples = samples[[2]]

#(1 - exp(-0.5))^(-1)

##Draw from aug UBA
try1 = aug_UBA(c(0,0.99),dt = 1e-3, N = 1e4)
ind = ((try1[ ,1])^2 + (try1[,2])^2) <=1
out = try1[ind,]
r = 1 - nrow(out)/nrow(try1)

###We run accept reject to compare with truth
# Define parameters
theta <- seq(0, 2 * pi, length.out = 500)  # Generate angles
x <- cos(theta)  # X-coordinates
y <- sin(theta)
frame = cbind(x,y)
pdf('tr_norm_circle.pdf')
plot(try1[,1],try1[,2], col = 'blue', main = paste('r = ',r))
points(samples[,1], samples[,2], col = 'green')
points(frame[,1] , frame[,2])
legend('topright', legend = c('boundary', 'AugUBA', 'AR '), 
                  col = c('black','blue','green'), pch = 1)

dev.off()

