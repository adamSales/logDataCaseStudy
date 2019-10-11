library(gridExtra)
library(optmatch)
library(RItools)
library(sandwich)
library(clubSandwich)
library(lmtest)
library(mgcv)
library(tidyverse)
library(rstan)
library(lme4)
traceplot <- rstan::traceplot
select <- dplyr::select

omar <- par()$mar


source('R/helperFunctions.r')
source('~/Box Sync/rcode/matchingFunctions.r')
source('~/Box Sync/rcode/hhh.r')
source('~/Box Sync/rcode/printXbal.r')

source('R/prelimStan.r')

## usage data availability

Nt <- sum(dat$treatment)

obsUse=sum(dat$field_id[dat$treatment==1]%in%x$field_id)

sdat <- makeStanDat(dat,x,missingUsage=FALSE)

source('R/covariateTable.r')

source('R/evaluateSmooths.r')

overallHintInfo <-
  xOrig%>%
  filter(field_id%in%dat$field_id)%>%
  mutate(hint=nhints1>0,err=nerrs1>0)%>%
  xtabs(~hint+err,data=.)



sdat2 <- with(
  sdat,
  list(
    N=nsecWorked,
    nstud=sum(Z),
    nsec=nsec,
    studentM=as.numeric(as.factor(studentM)),
    section=section,
    hint=grad
  )
)

###########################################################################################
### defining "high hint users"
###########################################################################################

### raw means
hintAvg <-
  with(sdat2,tibble(student=studentM,section=section,hint=hint,field_id=x$field_id))%>%
  group_by(section)%>%
  mutate(nprobSec=n(),nstudSec=n_distinct(student),propHintSec=mean(hint))%>%
  group_by(student)%>%
  mutate(nprobStud=n(),nsecStud=n_distinct(section),propHintStud=mean(hint),secDiffStud=mean(propHintSec))

hintAvgStud <- hintAvg%>%group_by(student)%>%slice(1L)

hintAvgStud%>%
  gather("xaxis","x",nprobStud,secDiffStud)%>%
  mutate(nsec=ifelse(xaxis=='nprobStud', median(nsecStud), nsecStud))%>%
  ggplot(aes(x,propHintStud,size=nsec))+
  geom_point()+
  facet_wrap(~xaxis,scales="free_x")


### rasch mixture model?
##NOT RUN
## source('R/raschMixture.r')
load('output/raschMixOutput.RData')
mixProp <- round(mean(hintAvgStud$classProb),1)

hintAvgStud$M <- hintAvgStud$propHintStud>=quantile(hintAvgStud$propHintStud,1-mixProp)
hintAvgStud$Mmod <- hintAvgStud$classProb>quantile(hintAvgStud$classProb,1-mean(hintAvgStud$classProb))
agree=with(hintAvgStud,mean(M==Mmod))


dat1 <- inner_join(dat,select(hintAvgStud,field_id,propHintStud,M))#%>%
#  mutate(M=propHintStud>.6)


##########################################################################################
### propensity score matching
##########################################################################################

##NOT RUN
## source('R/matching.r')

 psmod3.2 <- ## stan_
   glmer(M~(1|state)+(1|schoolid2)+(1|classid2)+sex+grade+race+ns(xirt,5)+spec+esl+gradeMIS+raceMIS+sexMIS+frlMIS,family=binomial,data=dat1)
save(psmod3.2,file='output/psmod.RData')
load('output/psmod.RData')

dist3 <- match_on(M~predict(psmod3.2,type='link'),within=exactMatch(M~schoolid2,data=dat1),data=dat1)
m3 <- fullmatch(dist3,data=dat1)

balMod <- xBalance(M~grade+race+xirt+spec+esl+state,strata=list(unmatched=NULL,matched=~m3),report=c('std.diff','z.score','chisquare.test'),data=droplevels(dat1))

pdf('output/balancePlot.pdf')
plot(balMod,colors=subwayPalette)
abline(v=c(-.05,-.25,.05,.25),lty='dotted')
dev.off()

dat1$match <- m3

save(dat1,file='data/dat1.RData')

##########################################################################################
### estimate treatment effects (obs. study)
##########################################################################################


out1 <- lm(Y~M+match,data=dat1)
out1vcv <- vcovHC(out1)
out1t <- coeftest(out1,out1vcv)[2,]
out1ci <- coefci(out1,'MTRUE',vcov.=out1vcv)

summary(out1.1 <- lm(Y~M+splines::ns(xirt,5)+race+grade+esl+frl+sex+match,data=dat1))$coef[1:5,]
out1.1vcv <- vcovHC(out1.1)
out1.1t <- coeftest(out1.1,out1.1vcv)[2,]
out1.1ci <- coefci(out1.1,'MTRUE',vcov.=out1.1vcv)

### ATE and TOT weights
dat1 <-
  dat1%>%
  group_by(match)%>%
  mutate(ntot=n(),ntrt=sum(M),nctl=ntot-ntrt,eff=mean(Y[M])-mean(Y[!M]),winv=(1/nctl+1/ntrt))


out2 <- update(out1.1,weights=winv)
out2vcv <- vcovHC(out2)
out2t <- coeftest(out2,out2vcv)[2,]
out2ci <- coefci(out2,'MTRUE',vcov.=out2vcv)


out3 <- update(out1.1,weights=ntrt*winv)
out3vcv <- vcovHC(out3)
out3t <- coeftest(out3,out3vcv)[2,]
out3ci <- coefci(out3,'MTRUE',vcov.=out3vcv)


### sensitivity analysis
X <- as.data.frame(model.matrix(psmod3.2)[,-1])
X <- X[,-grep('xirt',names(X))]
X$xirt <- dat1$xirt
rrr <- rhos(dat1$Y,X=X)

X$M <- dat1$M
ttt <- Tz(X=X,treatment='M')

interval(ttt['xirt'],rrr['xirt'],b=out1.1t[1], se=out1.1t[2],df=out1.1$df)
interval(ttt['xirt']/2,rrr['xirt'],b=out1.1t[1], se=out1.1t[2],df=out1.1$df)


##########################################################################################
### mediation analysis
##########################################################################################



### direct effects
yhat <- sapply(levels(dat1$match),function(m) with(subset(dat1,match==m & !M), mean(Y)))
dat1$yhat <- ifelse(dat1$M,yhat[as.character(dat1$match)],dat1$Y)
datC <- subset(dat,treatment==0)
datC$yhat <- datC$Y
dat2 <- rbind(dat1[,names(datC)],datC)

dir1 <- lmer(yhat~treatment+pair+(1|classid2)+(1|schoolid2),data=dat2)
dir2 <- update(dir1,.~.+poly(xirt,2)+race+grade+esl+frl+sex)

dir1.1 <- lm(yhat~treatment+pair,data=dat2)
dir1.1vcv <- vcovCR(dir1.1,dat2$schoolid2,type='CR2')
dir1.1t <- coeftest(dir1.1,dir1.1vcv)[2,]
dir1.1ci <- coefci(dir1.1,'treatment',vcov.=dir1.1vcv)

dir2.1 <- update(dir1.1,.~.+poly(xirt,2)+race+grade+esl+frl+sex)
dir2.1vcv <- vcovCR(dir2.1,dat2$schoolid2,type='CR2')
dir2.1t <- coeftest(dir2.1,dir2.1vcv)[2,]
dir2.1ci <- coefci(dir2.1,'treatment',vcov.=dir2.1vcv)








########## principal stratification
if(is.null(getOption('mc.cores'))) options(mc.cores = 6)
if(getOption('mc.cores')<6)  options(mc.cores = 6)

functions <- currentCode()

### NOT RUN
psmod1 <- stan('R/prinStratStan.stan',data=sdat,chains=6,iter=6000)
save(psmod1,sdat,functions,file='fitModels/hintPS.RData')
load('fitModels/hintPS.RData')



draws <- rstan::extract(psmod1)


############################
######## main effect plot
#########################


pdMain <- pdMod(psmod1)
## pdMain <- within(pdMain,
## {
##     b0 <- b0/pooledSD
##     b1 <- b1/pooledSD*iqr
##     xmin <- xmin/mean(iqr)
##     xmax <- xmax/mean(iqr)
##     ymin <- ymin/pooledSD
##     ymax <- ymax/pooledSD
## }
## )

#tikzDevice::tikz('mainEffects.tex', standAlone=T,
#     width=6,height=5)
p1 <- ggplot(pdMain)+
    geom_abline(aes(intercept=b0,slope=b1,group=id),color='red')+
    coord_cartesian(xlim=c(min(pdMain$xmin),max(pdMain$xmax)),
                    ylim=c(min(pdMain$ymin),max(pdMain$ymax)),expand=FALSE)+
    geom_line(aes(x=x,y=y,group=truthOrAvg,linetype=truthOrAvg,color=truthOrAvg,alpha=truthOrAvg),size=1.5)+
    xlab('$\\eta(1)$')+ylab('CTA1 Treatment Effect')+
    labs(group=NULL,color=NULL,linetype=NULL)+
    scale_color_manual(values=c('black','red','black'))+scale_linetype_manual(values=c('solid','solid','dotted'))+
    scale_alpha_manual(values=c(1,0,1),guide=FALSE)+theme(legend.position='top')+
    theme(text=element_text(size=15),legend.key.width=unit(.5,'in'))
#dev.off()
#tools::texi2dvi('mainEffects.tex', pdf = T, clean = T)


#########################################
### eta_T vs Y
#########################################
pooledSD <- with(sdat, sqrt(((sum(Z)-1)*var(Y[Z==1])+(sum(1-Z)-1)*var(Y[Z==0]))/(nstud-2)))
draw <- which.min(abs(draws$b1-mean(draws$b1)))
plotDat <- with(sdat,data.frame(Y=scale(Y,center=mean(Y[Z==0]),scale=pooledSD),
                                   eta=scale(draws$studEff[draw,],scale=IQR(draws$studEff[draw,])),
                                   Z=Z))

plotDat$treat <- ifelse(plotDat$Z==1,'Treatment','Control')
plotDat$slope <- (draws$a1[draw]+ifelse(plotDat$treat=='Control',0,draws$b1[draw]))*IQR(draws$studEff[draw,])/pooledSD
plotDat$int <- ifelse(plotDat$treat=='Control',0,draws$b0[draw]/pooledSD)

#plotDat <- within(plotDat, int <- int-( mean(int+slope*eta)-mean(plotDat$Y)))
#plotDat <- plotDat[order(plotDat$treat),]
plotDat$treat2 <- plotDat$treat

p2 <- ggplot(plotDat,aes(eta,Y,fill=treat,group=treat,color=treat))+geom_point(size=1)+
      coord_cartesian(xlim=quantile(plotDat$eta,c(0.005,0.995)),ylim=quantile(plotDat$Y,c(0.005,0.995)))+
      geom_abline(aes(intercept=int,slope=slope,color=treat),size=4,alpha=1)+#+scale_alpha_discrete(range=c(0.4,.8))+
      geom_abline(aes(intercept=int,slope=slope),color='black',size=2,alpha=1)+#+scale_alpha_discrete(range=c(0.4,.8))+
    scale_colour_manual(values=c('red','blue'))+
    labs(group=NULL,fill=NULL,alpha=NULL)+xlab('$\\eta(1)$')+
    ylab('Posttest Score')+theme(legend.position='top',text=element_text(size=15))+
        guides(color = guide_legend(title=NULL,override.aes=list(alpha=1,size=3),keywidth=3),linetype=guide_legend(title=NULL,keywidth=1,override.aes=list(size=1)))

tikzDevice::tikz(file = "output/psModel.tex",
  standAlone = T,
  width  = 6.5, height  = 4)
grid.arrange(p1,p2,nrow=1)
dev.off()
tools::texi2dvi('output/psModel.tex', pdf = T, clean = T)




