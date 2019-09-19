
makeStanDat <- function(dat,x,xInteract=FALSE,missingUsage=TRUE){
  stanDat <- list()

  if(!missingUsage)
    dat <- filter(dat,field_id%in%x$field_id|dat$treatment==0)

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


pdMod <- function(mod,row=1,column=1,func){
    draws <- rstan::extract(mod)
    samp <- seq(1,length(draws$b1),length=1000)
    Usamp <- draws$studEff[samp,]
    iqr <- apply(Usamp,1,IQR)
    studEff95 <- quantile(Usamp,c(0.025,0.975))
    Usamp[Usamp<studEff95[1] | Usamp>studEff95[2]] <- NA
    trtEff <- sweep(sweep(Usamp,1,draws$b1[samp],'*'),1,draws$b0[samp],'+')



    if(missing(func)){
        func <- function(x) mean(draws$b0)+mean(draws$b1)*x
        knownTruth <- FALSE
    } else knownTruth <- TRUE
    truth <- curve(func,from=studEff95[1],to=studEff95[2],n=length(samp)/3)
    avg <- curve(mean(draws$b0)+x*mean(draws$b1),
                 from=studEff95[1],to=studEff95[2],n=length(samp)/3)
    postDraw <- curve(mean(draws$b0)+x*mean(draws$b1),
                      from=studEff95[1],to=studEff95[2],n=length(samp)-length(truth$x)-length(avg$x))
    x <- c(postDraw$x,truth$x,avg$x)
    y <- c(postDraw$y,truth$y,avg$y)
    if(knownTruth) truthOrAvg <- c(rep('Posterior\nDraws',length(postDraw$x)),rep('True\nEffect',length(truth$x)),rep('Posterior\nAverage',length(avg$x))) else
     truthOrAvg <- c(rep('Posterior\nDraws',length(postDraw$x)),rep('Posterior\nAverage',length(avg$x)+length(truth$x)))

#    if(knownTruth) title <- paste('True Effe

    pd <- data.frame(b0=draws$b0[samp],b1=draws$b1[samp],id=1:length(samp),row=row,column=column,xmin=studEff95[1],xmax=studEff95[2],ymin=min(trtEff,na.rm=T),ymax=max(trtEff,na.rm=T),x=x,y=y,
                     truthOrAvg=truthOrAvg,
                     iqr=iqr)
    pd
}

modSum <- function(mod){
  summ <- summary(mod)$summary
  dims <- mod@par_dims
  for(pp in names(dims)){
    print(pp)

    if(!length(dims[[pp]])){
      print(summ[pp,])
    } else if(prod(dims[[pp]])<10){
      print(summ[startsWith(rownames(summ),pp),])
    } else{
      print(summary(summ[startsWith(rownames(summ),pp),'Rhat']))
    }
  }
}
