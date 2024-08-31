###################
###Euler Maruyama
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





##Euler Maruyama to get solution paths 
Euler_Maruyama = function(x0, eps, sigma,n ,dt)
{
  x = numeric(n)
  x[1] = x0
  for( i in 2:n)
  {
    mu = drift(x[i-1],eps)
    x[i] =  x[i-1]+ dt*mu+ sigma*sqrt(dt)*rnorm(1)
  }
  return(x)
}


##To get equilibrium distribution 
n = 1e5
reps = 100
store_mat = matrix(0 , nrow = reps, ncol = n)
for( i in 1:reps)
{
  store_mat[i,] = Euler_Maruyama(x0=3 , eps = 2 ,dt = 0.01, sigma = sqrt(2), n = n)
  print(i)
}


#Equilibrium distribution 
plot(density(store_mat[ 1,]) , main = "eps = 2" , ylim = c(0,1))
lines( v , pif(v,eps = 2) , col = 'red') #true distribution 
#legend('topright' , legend = c('Euler-Maruyama' , 'Truth') , col = c('black' , 'red'), lty = 1 )




##Visualising the sample paths of the SDE 
t = seq(1 , 50)
plot( t , store_mat[1 , 1: 50], col = "black"  , lty = 1,
      main ='Sample Paths of SDE' , type ='l' , 
      ylim =c(-5,5) ,ylab = 'X_t')

for( i in 2:10)
{
  lines(store_mat[i, 1:50] , col = i)
  Sys.sleep(1)
}



v = seq(-4,4,length.out = 100)
##Plots to understand why it is called equilibrium distribution 
plot(density(store_mat[,1]) , col = "black"  , lty = 1 , 
     xlim = c(-5,5) , ylim = c(0 , 1), 
     main ='EM: Normal' )
for( i in 2:100)
{
  lines(density(store_mat[ ,i]) , col = 'blue')
  Sys.sleep(1)
}


lines( v , pif(v, eps = 2) , col = 'red' , lty = 1 , lwd = 3) #true distribution 
legend('topright' , legend = c('EM' , 'Truth') , col = c('black' , 'red'), lty = 1 )
