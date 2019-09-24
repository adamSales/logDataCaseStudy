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
sdatDrop <- makeStanDat(dat,x,missingUsage=FALSE)

psmodDrop <- stan('R/prinStratStan.stan',data=sdatDrop,chains=6,iter=3000)

modSum(psmodDrop)

psmodDropSumm <- summary(psmodDrop)$summary

save(psmodDrop,sdatDrop,functions,file='fitModels/psmodDrop.RData'); rm(psmodDrop); gc()
