## Euler-Maruyama Numerical Scheme for Non Lipschitz drift 
##To calculate the scaling constant in pdf 
samp = rnorm(1e5)
c = mean(abs(samp)^4)
## To calculate pdf of distribution 
pdf = function(t)
{
  return(exp(-(abs(t)^4)/4)/c) 
}
fooc = pdf(0)


##Calculates Gradient
grad_f = function(t)
{
  return(-t^3)
}

##Euler Maruyama to get solution paths 
Euler_Maruyama = function(x0, eps, sigma,n)
{
  x = numeric(n)
  x[1] = x0
  for( i in 2:n)
  {
    x[i] =  x[i-1]+ eps/2*grad_f(x[i-1])+ sigma*sqrt(eps)*rnorm(1)
  }
  return(x)
}
foo = Euler_Maruyama(x0=3 , eps = 0.1 , sigma = 3, n = 1000)
n = 1000
store_mat = matrix(0 , nrow = 100, ncol = n)

##To get equilibrium distribution 

for( i in 1:100)
{
  store_mat[i,] = Euler_Maruyama(x0=3 , eps = 0.1 , sigma = 3, n = n)
}

##Visualising the sample paths of the SDE 
t = seq(1 , 50)
plot( t , store_mat[1 , 1: 50], col = "black"  , lty = 1,
      main ='Sample Paths of SDE' , type ='l' , 
      ylim =c(-5,5) ,ylab = 'X_t')

for( i in 2:50)
{
  lines(store_mat[i, 1:50] , col = i)
  Sys.sleep(1)
}



v = seq(-4,4,length.out = 100)
##Plots to understand why it is called equilibrium distribution 
plot(density(store_mat[,1]) , col = "black"  , lty = 2 , 
     xlim = c(-5,5) , ylim = c(0 , 1), 
     main ='Equilibrium Distribution :Non Lipschitz Drift' )
for( i in 2:100)
{
  lines(density(store_mat[ ,i]) , col = 'blue')
  Sys.sleep(1)
}

lines(v, pdf(v), lty = 1 , lwd = 2 , col = 'red')
legend('topright' , legend = c('Euler-Maruyama' , 'Truth') , col = c('black' , 'red'), lty = 1 )

#Equilibrium distribution 
plot(density(store_mat[ ,-1]) , main = "Non-Lipschitz Eq" , ylim = c(0,0.5))
lines( v , pdf(v) , col = 'red') #true distribution 
legend('topright' , legend = c('Euler-Maruyama' , 'Truth') , col = c('black' , 'red'), lty = 1 )

