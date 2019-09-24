if(!exists('datF')|!exists('xF'))
  source('R/fakePrep.r')

functions <- currentCode()

########### quadratic effects
datQuad <- within(datF,{
    te <- -(U-mean(U))^2
    te <- te/sd(te)*0.1
    te <- te-mean(te)+0.13
    Y[treatment==1] <- Y[treatment==1]+te[treatment==1]
})

sdatFquad <- makeStanDat(datQuad,xF)



quadEff <- stan('R/prinStratStan.stan',data=sdatFquad,iter=3000,chains=6)

##cat('\n\n\n\n',rep('-',40),'\n','TRUTH: QUADRATIC EFFECT',rep('-',40),'\n\n\n')
##print(constEff,c('a0','a1','b0','b1'),c(0.05,0.95))

modSum(quadEff)

save(quadEff,sdatFquad,functions,file='fitModels/quadEff.RData'); rm(quadEff); gc();








