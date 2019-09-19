require(rstan)
require(tidyverse)

source('R/helperFunctions.r')

if(is.null(getOption('mc.cores'))) options(mc.cores = 4)
if(getOption('mc.cores')<4)  options(mc.cores = 4)

functions <- currentCode()

############################
## robustness checks
############################

load('data/hintsData.RData')

### drop treatment cases w/o usage data (same dataset as obs study and mediation analysis)
sdatDrop <- makeStanDat(dat,x,missingUsage=FALSE)

psmodDrop <- stan('R/prinStratStan.stan',data=sdatDrop)

modSum(psmodDrop)

psmodDropSumm <- summary(psmodDrop)$summary

save(psmodDrop,sdatDrop,functions,file='output/psmodDrop.RData'); rm(psmodDrop); gc()
