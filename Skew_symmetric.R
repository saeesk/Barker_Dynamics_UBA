###################
#### SKew Symmetric 
###################
## Calculate the scaling constant for pdf 
calcz = function(eps) {
  x = rexp(1e5)
  Y = x - (x^(2+eps))/(2+eps)
  avg = mean(exp(Y))
  z = 2 * avg
  return(z)
}


## Calculate the pdf 
pif = function(t , eps) {
  z = calcz(eps)
  num = (abs(t))^(2+eps)
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

##To get equilibrium distribution 
n = 1e5
store_mat <- sapply(1:1000, function(i) {
  result <- UBarker(x0 = 0, eps = 2, sigma = sqrt(2), dt = 0.01, N = n)
  print(i)
  return(result)
})

store_mat <- matrix(store_mat, nrow = 1000, ncol = n, byrow = TRUE)


v = seq(-3,3,by = 0.01)
#Equilibrium distribution 
plot(density(store_mat[ ,-1]) , main = "Eps = 2" , ylim = c(0,0.9))
lines(v, pif(v,2), col = 'Blue')
#legend('topright' , legend = c('Euler-Maruyama' , 'Truth') , col = c('black' , 'red'), lty = 1 )


##Visualising the sample paths of the SDE 
t = seq(1 , 50)
plot( t , store_mat[1 , 1: 50], col = "black"  , lty = 1,
      main ='Sample Paths - UBarker' , type ='l' , 
      ylim =c(-6,6) ,ylab = 'X_t')

for( i in 2:10)
{
  lines(store_mat[i, 1:50] , col = i)
  Sys.sleep(1)
}



v = seq(-4,4,length.out = 100)
##Plots to understand why it is called equilibrium distribution 
plot(density(store_mat[,1]) , col = "black"  , lty = 2 , 
     xlim = c(-5,5) ,
     main ='Equilibrium Distribution : Lipschitz Drift' )
for( i in 2:100)
{
  lines(density(store_mat[ ,i]) , col = 'blue')
  Sys.sleep(1)
}

lines( v , dnorm(v) , col = 'red' , lty = 1 , lwd = 3) #true distribution 
legend('topright' , legend = c('Euler-Maruyama' , 'Truth') , col = c('black' , 'red'), lty = 1 )
