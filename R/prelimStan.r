### only look at year 2
### and only at students who weren't in  year 1
library(splines)
library(rstan)
library(tidyverse)
memory.limit(50000)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

discard <- TRUE
afterError <- FALSE

#save(list=ls(),file=paste0('prev',Sys.time(),'.RData'))
#rm(list=ls())

load('~/Box Sync/CT/data/RANDstudyData/HSdata.RData')
load('~/Box Sync/CT/data/problemLevelUsageData/probLevelData.RData')

for(nn in c('curriculum','unit','section','Prob1')) x[[nn]] <- tolower(x[[nn]])

datOrig <- dat
xOrig <- x


## ### estimated eta using only treatment group
## library(jagstools)
## load('~/Google Drive/CTmodels/fullUsage.RData')
## usageMod <- mod
## U <- jagsresults(usageMod,'studEff')[,1]

#dataPrep <- function(dat,x,discard=TRUE,afterError=FALSE){
 dat <- dat[!dat$field_id%in%dat$field_id[dat$year==1],]
 #print(table(dat$year))

 dat <- droplevels(dat)

 x <- x[x$field_id%in%dat$field_id[dat$treatment==1],]

 ### just look at algebra I sections--- likely different advance patterns in other curricula
### make sure to keep algebra i units that are also part of other curricula
 algUnit <- unique(x$unit[x$curriculum=='algebra i'])
 x <- subset(x,unit%in%algUnit)

### take out bridge to algebra, algebra ii, geometry
 x <- x[-grep('ii|geo|bridge',x$curriculum),]

 x$hint <- x$nhints1>0 #| x$nerrs1>0

 x <- subset(x,nhints1>0 |x$nerrs1>0)

 x <- droplevels(x)


 ### discard some pairs
 ## discard treatment schools with no usage data
 ## (do a robustness check with everything left in afterwards)
 percUse <- function(scl)
    length(intersect(unique(dat$field_id[dat$schoolid2==scl]),unique(x$field_id)))/
        length(unique(dat$field_id[dat$schoolid2==scl]))

 obsUse <- vapply(unique(dat$schoolid2[dat$treatment==1]),percUse,1)

 obsUse <- unique(dat$schoolid2[dat$treatment==1])[obsUse>0.1]
 obsUse <- c(as.character(obsUse),as.character(unique(dat$schoolid2[dat$treatment==0])))

 if(discard)
  dat <- dat[dat$schoolid2%in%obsUse,]

 ## discard pairs with only a treatment or a control school
 pairTrtTab <- with(dat,table(pair,treatment))
 trtVar <- apply(pairTrtTab,1,prod)
 dat <- dat[trtVar[dat$pair]>0,]

 x <- x[x$field_id%in%dat$field_id,]

 aaa <- aggregate(x$hint,by=list(section=x$section),FUN=mean)
 aaa$n <- as.vector(table(x$section))

 if(discard) x <- subset(x,section%in%aaa$section[aaa$n>100 & aaa$x<1])

# x <- x%>%group_by(field_id,section)%>%summarize(grad=min(ass,na.rm=TRUE))

 x <- droplevels(x)
 dat <- droplevels(dat)

 if(afterError)
    x <- x%>%
    mutate(ts1=as.POSIXct(ts1,format='%m/%d/%y %H:%M'))%>%
    arrange(ts1)%>%
    group_by(field_id,section)%>%
    filter(n()>1)%>%
    mutate(prevErr=lag(nerrs1))%>%
    filter(!is.na(prevErr),prevErr>0)

 save(dat,x,file='hintsData.RData')





makeStanDat <- function(dat,x,xInteract=FALSE,missingUsage=TRUE){
  stanDat <- list()

  if(!missingUsage)
    dat <- filter(dat,field_id%in%x$field_id)

  x <- droplevels(x)
 dat <- droplevels(dat)

 x$section <- as.factor(x$section)

 stanDat$nsecWorked <- nrow(x)
 stanDat$nstud <- nrow(dat)
 stanDat$nteacher <- nlevels(dat$teachid2)
 stanDat$nschool <- nlevels(dat$schoolid2)
 stanDat$nsec <- nlevels(x$section)
 stanDat$npair <- nlevels(dat$pair)

 stanDat$teacher <- as.numeric(as.factor(dat$teachid2))
 stanDat$pair <- as.numeric(dat$pair)
 stanDat$school <- as.numeric(dat$school)
 stanDat$studentM <- seq(stanDat$nstud)[match(x$field_id,dat$field_id)]
 stanDat$section <- as.numeric(x$section)

 stanDat$grad <- as.numeric(x$hint)

 X <- model.matrix(~poly(xirt,2)+race+sex+spec,data=dat)[,-1]
 if(xInteract) X <- model.matrix(~poly(xirt,2)*(race+sex+spec)+(race+sex+spec)^2+state,data=dat)[,-1]
 stanDat$X <- scale(X)
 stanDat$ncov <- ncol(X)


 stanDat$Z <- as.numeric(dat$treatment)

 stanDat$Y <- dat$Y

 stanDat
}


sdat <- makeStanDat(dat,x)




