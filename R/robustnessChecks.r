require(rstan)
require(tidyverse)

source('R/helperFunctions.r')

if(is.null(getOption('mc.cores'))) options(mc.cores = 6)
if(getOption('mc.cores')<6)  options(mc.cores = 6)

functions <- currentCode()

############################
## robustness checks
############################

load('data/hintsData.RData')

### drop treatment cases w/o usage data (same dataset as obs study and mediation analysis)
sdatFull <- makeStanDat(dat,x,missingUsage=TRUE)

psmodFull <- stan('R/prinStratStan.stan',data=sdatFull,chains=6,iter=6000)

modSum(psmodFull)

psmodFullSumm <- summary(psmodFull)$summary

save(psmodFull,sdatFull,functions,file='fitModels/psmodFull.RData'); rm(psmodFull); gc()
