#######################
###Code to save plots
#######################
##Warm Start MSE comparisons
pdf('WarmStart_AgUBA.pdf')
warm = read.csv('errors_warmstart.csv')
UBA = warm[ ,2]
AgUBA = warm[ ,3]
timesteps = seq(0.005, 0.05, by = 0.005)
plot(timesteps , UBA , type = 'o', col = 'black',
     ylim =c(0.0008,0.0455), ylab = 'MSE', 
     main = 'Warm Start')
lines(timesteps, AgUBA, type = 'o',col = 'red')
legend('topleft' , legend = c('UBA' , 'AgUBA') , 
       col = c('black','red'), lty = 1)
dev.off()

##Fixed Strat MSE comparisons
pdf('FixedStart_AgUBA.pdf')
fixed= read.csv('errors_fixedstart.csv')
UBA = fixed[ ,2]
AgUBA = fixed[ ,3]
timesteps = seq(0.005, 0.05, by = 0.005)
plot(timesteps , UBA , type = 'o', col = 'black',
     ylim =c(0.0008,0.0455), ylab = 'MSE', 
     main = 'Fixed Start')
lines(timesteps, AgUBA, type = 'o',col = 'red')
legend('topleft' , legend = c('UBA' , 'AgUBA') , 
       col = c('black','red'), lty = 1)
dev.off()


###Plots for Ratios
rw_UBA = warm[ ,4]
rf_UBA=  fixed[ ,4]
rw_AgUBA = warm[ ,5]
rf_AgUBA = fixed[ ,5] 

##Warm star ratio comparison
pdf('ratio_warmstart.pdf')
plot(timesteps , rw_AgUBA , col ='red', main = 'Warm Start',
     ylim = c(0.2,0.6),type = 'o', ylab = 'ratio')
lines(timesteps , rw_UBA[,2] ,col = 'black',  type = 'o')
abline(h = 0.5 , col = 'blue', lty = 2)
legend('bottomright' , legend = c('UBA' , 'AgUBA'), lty = 1 , col = c('black', 'red'))
dev.off()

##Fixed start ratio comparison 
pdf('ratio_fixedstart.pdf')
plot(timesteps , rf_AgUBA , col ='red', main = 'Fixed Start',
     ylim = c(0.2,0.6),type = 'o', ylab = 'ratio')
lines(timesteps , rf_UBA[,2] ,col = 'black',  type = 'o')
abline(h = 0.5 , col = 'blue', lty = 2)
legend('bottomright' , legend = c('UBA' , 'AgUBA'), lty = 1 , col = c('black', 'red'))
dev.off()

