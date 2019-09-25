if(!exists('datF')|!exists('xF'))
  source('R/fakePrep.r')

functions <- currentCode()

#########################################################
### Compile Stan data as before
#####################################################
sdatF <- makeStanDat(datF,xF)


noEff <- stan('R/prinStratStan.stan',data=sdatF,iter=3000,chains=6)

modSum(noEff)

##cat('\n\n\n\n',rep('-',40),'\n','TRUTH: NO EFFECT\n',rep('-',40),'\n\n\n')
##print(noEff,c('a0','a1','b0','b1'),c(0.05,0.95))
save(noEff,sdatF,functions,file='fitModels/noEffect.RData'); rm(noEff); gc()




