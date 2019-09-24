library(rstan)
library(tidyverse)

options(mc.cores = 6)

#source('R/prelimStan.r')
source('R/helperFunctions.r')

load('data/hintsData.RData')
################################################
## Fake Models
##############################################

### use data, E[studEff] from main model
load('fitModels/hintPS.RData')
dat$U <- colMeans(rstan::extract(psmod1,'studEff')[[1]])
### U has a couple outliers:
## newRng <- quantile(dat$U,c(0.25,0.75))+c(-1,1)*2*IQR(dat$U)
## dat$U[dat$U>newRng[2]] <- newRng[2]
## dat$U[dat$U<newRng[1]] <- newRng[1]

effs <- rstan::extract(psmod1,c('b0','b1'))

rm(psmod1); gc()

##### make fake data:

### first delete control schools
datF <- subset(dat,treatment==1)

### now double the dataset
datF$schoolid2 <- as.character(datF$schoolid2)
datF$teachid2 <- as.character(datF$teachid2)
dat2 <- datF
dat2$schoolid2 <- paste0(dat2$schoolid2,'Fake')
dat2$teachid2 <- paste0(dat2$teachid2,'Fake')
dat2$field_id <- dat2$field_id*100+99
dat2$treatment <- 0

## ### take bootstrap samples within each school
## for(scl in unique(datF$schoolid2)){
##     ind <- which(datF$schoolid2==scl)
##     datF[ind,] <- datF[sample(ind,length(ind),replace=TRUE),]
##     dat2[ind,] <- dat2[sample(ind,length(ind),replace=TRUE),]
## }

datF <- rbind(datF,dat2)

datF$schoolid2 <- factor(datF$schoolid2)
datF$teachid2 <- factor(datF$teachid2)

### now delete usage data for "control" group
xF <- x[x$field_id%in%datF$field_id[datF$treatment==1],]


