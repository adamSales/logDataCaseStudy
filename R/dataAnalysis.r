library(splines)
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
library(estimatr)
traceplot <- rstan::traceplot
select <- dplyr::select
options(stringsAsFactors=FALSE)

omar <- par()$mar

### this is a tool to skip the most computationally-expensive steps but run everything else:
runFull <- FALSE
cat('\n\n----------------------------------------------------\n')
cat('Warning: runFull=',runFull,'\n')
cat('----------------------------------------------------\n\n')

source('R/helperFunctions.r')
source('~/Box Sync/rcode/matchingFunctions.r')
source('~/Box Sync/rcode/hhh.r')
source('~/Box Sync/rcode/printXbal.r')

load('~/Box Sync/CT/data/RANDstudyData/HSdata.RData') ## for "error"
if(runFull) source('R/prelimStan.r') else load('data/hintsData.RData')

## for printing nice CIs
printci <- function(ci) paste0('[',paste(sprintf('%.2f',ci),collapse=','),']')

### I don't know why this is necessary...
map_dfr <- function(.x, .f, ...){
  out <- map(.x,.f,...)
  as.data.frame(do.call('rbind',out))
}

citopm <- function(ci) paste(sprintf('%.2f',c(median(ci),(ci[2]-ci[1])/2)),collapse='$\\pm$')

## usage data availability

Nt <- sum(dat$treatment)

obsUse=sum(dat$field_id[dat$treatment==1]%in%x$field_id)

sdat <- makeStanDat(dat,x,missingUsage=FALSE)

Ntot <- sdat$nstud
NtFinal <- sum(sdat$Z)
nteacher <- sdat$nteacher
nschool <- sdat$nschool

source('R/covariateTable.r')

#source('R/evaluateSmooths.r')

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

cat('
###########################################################################################
### defining "high hint users"
###########################################################################################
',as.character(Sys.time()),'\n')

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

#runFull <- FALSE
### rasch mixture model?
##NOT RUN
if(runFull) source('R/raschMixture.r') else load('output/raschMixOutput.RData')

#runFull <- TRUE

mixProp <- round(mean(hintAvgStud$classProb),1)

hintAvgStud$M <- hintAvgStud$propHintStud>=quantile(hintAvgStud$propHintStud,1-mixProp)
hintAvgStud$Mmod <- hintAvgStud$classProb>quantile(hintAvgStud$classProb,1-mean(hintAvgStud$classProb))
agree=with(hintAvgStud,mean(M==Mmod))


dat1 <- inner_join(dat,select(hintAvgStud,field_id,propHintStud,M))#%>%
#  mutate(M=propHintStud>.6)

cat('
##########################################################################################
### propensity score matching
##########################################################################################
',as.character(Sys.time()),'\n')

if(runFull){
  psmod3.2 <- ## stan_
    glmer(M~(1|state)+(1|schoolid2)+(1|classid2)+sex+grade+race+ns(xirt,5)+spec+esl+gradeMIS+raceMIS+sexMIS+frlMIS,family=binomial,data=dat1)
  save(psmod3.2,file='output/psmod.RData')
} else load('output/psmod.RData')

dist3 <- match_on(M~predict(psmod3.2,type='link'),within=exactMatch(M~schoolid2,data=dat1),data=dat1)
m3 <- fullmatch(dist3,data=dat1)

dat1$pretest <- dat1$xirt

balMod <- xBalance(M~grade+race+pretest+spec+esl+state,strata=list(`Before Matching`=NULL,`After Matching`=~m3),report=c('std.diff','z.score','chisquare.test'),data=droplevels(dat1))

pdf('output/balancePlot.pdf')
plot(balMod,colors=subwayPalette)
abline(v=c(-.05,-.25,.05,.25),lty='dotted')
dev.off()

dat1$match <- m3

save(dat1,file='data/dat1.RData')

cat('
##########################################################################################
### estimate treatment effects (obs. study)
##########################################################################################
',as.character(Sys.time()),'\n')

## anova <- lm(Y~M+match,data=dat1)
## anovavcv <- vcovHC(anova,'HC2')
## anovat <- coeftest(anova,anovavcv)[2,]
## anovaci <- coefci(anova,'MTRUE',vcov.=anovavcv)

## summary(ancova <- lm(Y~M+splines::ns(xirt,5)+race+grade+esl+frl+sex+match,data=dat1))$coef[1:5,]
## ancovavcv <- vcovHC(ancova,'HC2')
## ancovat <- coeftest(ancova,ancovavcv)[2,]
## ancovaci <- coefci(ancova,'MTRUE',vcov.=ancovavcv)


## ### ben's method
## olsC <- update(ancova,.~.-M-match,subset=!M)
## yhat <- predict(olsC,dat1)
## dat1$e <- dat1$Y-yhat

effWithWeights <-
  dat1%>%
  group_by(match)%>%
  summarize(
    eff=mean(Y[M])-mean(Y[!M]),
    #adj=mean(e[M])-mean(e[!M]),
    n=n(),
    ntrt=sum(M),
    nctl=n-ntrt,
    prec=1/(1/ntrt+1/nctl)
  )

effects <- list()
for(est in c('ate','tot','prec'))
  for(ca in c('eff'))#,'adj'))
    effects[[paste0(est,'_',ca)]] <-
      lm_robust(
        as.formula(paste(ca,'~1')),
        effWithWeights,
        weights=if(est=='ate') n else if(est=='tot') ntrt else prec
      )[c('coefficients', 'std.error','p.value', 'conf.low', 'conf.high')]


### sensitivity analysis
X <- as.data.frame(model.matrix(psmod3.2)[,-1])
X <- X[,-grep('xirt',names(X))]
X$xirt <- dat1$xirt
rrr <- rhos(dat1$Y,X=X)

X2 <- X
X2$M <- dat1$M
ttt <- Tz(X=X2,treatment='M')

multipliers <- map2_dbl(rrr,ttt,~MEmult(abs(.y),.x,Inf,1.96))

raceAnova <- anova(lm(M~.,data=X2),lm(M~.-raceBlackMulti-raceHispAIAN,data=X2))

m1 <- lm(dat1$Y~as.matrix(X))
m2 <- lm(dat1$Y~as.matrix(select(X,-raceBlackMulti,-raceHispAIAN)))
r1 <- summary(m1)$r.squared
r2 <- summary(m2)$r.squared
raceR <- ((1-r2)-(1-r1))/(1-r2)

sens <- map(effects,
  ~rbind(
    pre=interval(ttt['xirt'],rrr['xirt'],b=.$coefficients,se=.$std.error,df=Inf),
      #interval(ttt['xirt']/2,rrr['xirt']/2,b=.$coefficients,se=.$std.error,df=Inf),
      race=interval(sqrt(raceAnova$F)[2],raceR,b=.$coefficients,se=.$std.error,df=Inf))
)


cat('
##########################################################################################
### mediation analysis
##########################################################################################
',as.character(Sys.time()),'\n')


### direct effects
yhat <- sapply(levels(dat1$match),function(m) with(subset(dat1,match==m & !M), mean(Y)))
dat1$yhat <- ifelse(dat1$M,yhat[as.character(dat1$match)],dat1$Y)
datC <- subset(dat,treatment==0)
datC$yhat <- datC$Y
dat2 <- rbind(dat1[,names(datC)],datC)

#dir1 <- lmer(yhat~treatment+pair+(1|classid2)+(1|schoolid2),data=dat2)
#dir2 <- update(dir1,.~.+poly(xirt,2)+race+grade+esl+frl+sex)

dir1.1 <- lm(yhat~treatment+pair,data=dat2)
dir1.1vcv <- vcovCR(dir1.1,dat2$schoolid2,type='CR2')
dir1.1t <- coef_test(dir1.1,dir1.1vcv,coefs='treatment')
dir1.1ci <- conf_int(dir1.1,dir1.1vcv,coefs='treatment')
dir1.1ci <- c(dir1.1ci$CI_L[1],dir1.1ci$CI_U[1])

dir2.1 <- update(dir1.1,.~.+poly(xirt,2)+race+grade+esl+frl+sex)
dir2.1vcv <- vcovCR(dir2.1,dat2$schoolid2,type='CR2')
dir2.1t <- coef_test(dir2.1,dir2.1vcv,coefs='treatment')
dir2.1ci <- conf_int(dir2.1,dir2.1vcv,coefs='treatment')
dir2.1ci <- c(dir2.1ci$CI_L[1],dir2.1ci$CI_U[1])





runFull <- FALSE

cat('
############################
########## principal stratification
############################
',as.character(Sys.time()),'\n')

if(is.null(getOption('mc.cores'))) options(mc.cores = 6)
if(getOption('mc.cores')<6)  options(mc.cores = 6)

functions <- currentCode()

cat(ifelse(runFull,'running','loading'),'ps model',as.character(Sys.time()),'\n')
### NOT RUN
if(runFull){
  psmod1 <- stan('R/prinStratStan.stan',data=sdat,chains=6,iter=6000)
  save(psmod1,sdat,functions,file='fitModels/hintPS.RData')
} else load('fitModels/hintPS.RData')

cat('extracting ps model',as.character(Sys.time()),'\n')
draws <- rstan::extract(psmod1)

cat('
############################
######## main effect plot
#########################
',as.character(Sys.time()),'\n')

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

cat('
#########################################
### eta_T vs Y
#########################################
',as.character(Sys.time()),'\n')

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

cat('
#########################################
### saving
#########################################
',as.character(Sys.time()),'\n')

big <- sapply(ls(),function(x) object.size(get(x))>1e7)
save(list=names(big)[!big],file='output/analysis.RData')

cat('done', as.character(Sys.time()),'\n')
