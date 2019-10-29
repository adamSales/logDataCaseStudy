require(tidyverse)
source('R/helperFunctions.r')

library(bayesplot)

if(!exists('draws')){
  load('fitModels/hintPS.RData')
  draws <- rstan::extract(psmod1)
}
if(!exists('sdat')){
  load('data/hintsData.RData')
  sdat <- makeStanDat(dat,x)
}

### make Yrep
set.seed(613)
samp <- sample(1:nrow(draws$studEff),size=1000)
draws1k <- lapply(draws, function(x) if(NCOL(x)==1) x[samp] else x[samp,])

ppc_dens_overlay(sdat$Y,draws1k$Yrep)+ggtitle('Pooled Treatment and Control')
ggsave('output/ppcYoverall.jpg')

ppc_dens_overlay(sdat$Y[sdat$Z==1],draws1k$Yrep[,sdat$Z==1])+ggtitle('Treatment Group')
ggsave('output/ppcYtrt.jpg')

ppc_dens_overlay(sdat$Y[sdat$Z==0],draws1k$Yrep[,sdat$Z==0])+ggtitle('Control Group')
ggsave('output/ppcYctl.jpg')

### residual plots
fitd <- colMeans(draws$Yrep)
qplot(fitd,sdat$Y-fitd,xlab='Fitted Values',ylab='Residuals',main='Overall')
ggsave('output/residPlotOverall.jpg')

qplot(fitd[sdat$Z==1],sdat$Y[sdat$Z==1]-fitd[sdat$Z==1],xlab='Fitted Values',ylab='Residuals',main='Treatment Group')
ggsave('output/residPlotTrt.jpg')

qplot(fitd[sdat$Z==0],sdat$Y[sdat$Z==0]-fitd[sdat$Z==0],xlab='Fitted Values',ylab='Residuals',main='Control Group')
ggsave('output/residPlotCtl.jpg')


pdf('output/qqPlots.pdf', height=3,width=6)
par(mfrow=c(1,2))
qqnorm(sdat$Y[sdat$Z==1]-fitd[sdat$Z==1],main='Treatment Group')
qqline(sdat$Y[sdat$Z==1]-fitd[sdat$Z==1])
qqnorm(sdat$Y[sdat$Z==0]-fitd[sdat$Z==0],main='Control Group')
qqline(sdat$Y[sdat$Z==0]-fitd[sdat$Z==0])
dev.off()


### usage model
set.seed(613)
samp <- sample(1:1000,9)
lp <- with(c(sdat,draws1k),studEff[samp,studentM]+secEff[samp,section])
prob <- exp(lp)/(1+exp(lp))
p <- ppc_error_binned(sdat$grad,prob)
ggplot2::ggsave('output/binnedplot.pdf',p)

## ### ppc m-bar
## mbarRep <- apply(prob,1,function(p) aggregate(rbinom(length(sdat$grad),1,p),by=list(sdat$studentM),FUN=mean)$x)
## mbar <- aggregate(sdat$grad,by=list(sdat$studentM),FUN=mean)$x
## ppc_dens_overlay(mbar,t(mbarRep))
## ggsave('output/mbarPPC.jpg')

