if(!exists('datF')|!exists('xF'))
  source('R/fakePrep.r')

sdatF <- makeStanDat(datF,xF)

functions <- currentCode()
######## constant effect
########
sdatFConst <- within(sdatF,{
 te <- rnorm(sum(Z),0.18,0.1)
 Y[Z==1] <- Y[Z==1]+ te
 }
)

constEff <- stan('R/prinStratStan.stan',data=sdatFConst,iter=3000,chains=6)

modSum(constEff)

##cat('\n\n\n\n',rep('-',40),'\n','TRUTH: CONSTANT EFFECT ATE=0.18\n',rep('-',40),'\n\n\n')
##print(constEff,c('a0','a1','b0','b1'),c(0.05,0.95))
save(constEff,sdatF,functions,file='fitModels/constEff.RData')
