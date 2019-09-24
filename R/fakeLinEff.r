if(!exists('datF')|!exists('xF'))
  source('R/fakePrep.r')

functions <- currentCode()

##################################
######### linear effect ##########
##################################
datLin <- datF
datLin$te <- mean(effs$b0)+mean(effs$b1)*datF$U
datLin$Y[datLin$treatment==1] <- datLin$Y[datLin$treatment==1]+datLin$te[datLin$treatment==1]

sdatFlin <- makeStanDat(datLin,xF)

linEff <- stan('R/prinStratStan.stan',data=sdatFlin,iter=3000,chains=6)
##cat('\n\n\n\n',rep('-',40),'\n','TRUTH: LINEAR EFFECT b1=',round(mean(effs$b1),2),'\n',rep('-',40),'\n\n\n')
##print(linEff,c('a0','a1','b0','b1'),c(0.05,0.95))

modSum(linEff)

save(linEff,sdatFlin,functions,file='fitModels/linEff.RData');rm(linEff);gc();



