##################################
###How Epsilon changes PDFs
################################
## Euler-Maruyama Numerical Scheme for SDE 
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


##Comparing pdfs
v = seq(-3,3,by =0.01)
eps = 0:5

pdf('PDF_eps.pdf')
plot(v , pif(v,eps = 0), col = 1 , type = 'l' ,
     ylab ="f(x)" , xlab = 'x',
     main= "Effect of epsilon on PDF",
     title = expression(epsilon))

for(ep in 1:5)
 {
   lines(v,pif(v,ep), col = ep+1 )
 }
legend('topright' , legend = eps , col = 1:6 , lty = 1,
       title = expression(epsilon))

dev.off()
##Checking effect on tails
v = seq(1,4, by = 0.01)
eps = 0:5
pdf('Tails_eps.pdf')
plot(v , pif(v,eps = 0), col = 1 , type = 'l' ,
     ylab ="f(x)" , xlab = 'x',
     main = "Effect of epsilon on tails")
for(ep in 1:5)
{
  lines(v,pif(v,ep), col = ep+1 )
}
legend('topright' , legend = eps , col = 1:6 , lty = 1, 
       title = expression(epsilon))
dev.off()

##Comapring drifts
v = seq(-3,3,by =0.01)
eps = 0:5
pdf('Drift_eps.pdf')
plot(v , drift(v,eps = 0), col = 1 , type = 'l' ,
     ylab ="Drift" , xlab = 'X',
     main= "Effect of epsilon on Drift")
for(ep in 1:5)
{
  lines(v,drift(v,ep), col = ep+1 )
}
lines(v , rep(0, length(v)) , lty = 2 , col = 'grey')
lines( rep(0 , length(v)) , v, lty = 2 , col = 'grey')
legend('topright' , legend = eps , col = 1:6 , lty = 1, 
       title = expression(epsilon))
dev.off()


