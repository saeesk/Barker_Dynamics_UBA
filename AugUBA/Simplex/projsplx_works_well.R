####The following code proves that 
#### The function projsplx 
####is working fine 
samp = matrix(rnorm(1e4*2), ncol = 2)
pro_samp = matrix(nrow = 1e4 , ncol = 2)
for( i in 1: 1e4)
{
  pro_samp[i, ] = projsplx(samp[i,])
}
  
plot(samp, col = 'green')
points(pro_samp[,1], pro_samp[,2], col = 'blue', type = 'l')
abline(h = 1)
abline(h = 0)
abline(v = 0)
abline(v = 1)
