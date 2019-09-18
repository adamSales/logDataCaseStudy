library(bayesplot)
### make Yrep
set.seed(613)
samp <- sample(1:nrow(draws$studEff),size=1000)
draws1k <- lapply(draws, function(x) if(NCOL(x)==1) x[samp] else x[samp,])

ppc_dens_overlay(sdat$Y,draws1k$Yrep)+ggtitle('Pooled Treatment and Control')
ggsave('output/ppcYoverall.jpg')

ppc_dens_overlay(sdat$Y[sdat$Z==1],t(Yrep)[,sdat$Z==1])+ggtitle('Treatment Group')
ggsave('output/ppcYtrt.jpg')

ppc_dens_overlay(sdat$Y[sdat$Z==0],t(Yrep)[,sdat$Z==0])+ggtitle('Control Group')
ggsave('output/ppcYctl.jpg')

### residual plots
fitd <- rowMeans(Yrep)
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
lp <- with(c(sdat,draws1k),studEff[,studentM]+secEff[,section])
prob <- exp(lp)/(1+exp(lp))
p <- ppc_error_binned(sdat$grad,prob[samp,])
ggplot2::ggsave('output/binnedplot.pdf',p)

### ppc m-bar
mbarRep <- apply(prob,1,function(p) aggregate(rbinom(length(sdat$grad),1,p),by=list(sdat$studentM),FUN=mean)$x)
mbar <- aggregate(sdat$grad,by=list(sdat$studentM),FUN=mean)$x
ppc_dens_overlay(mbar,t(mbarRep))
ggsave('output/mbarPPC.jpg')



############################
## robustness checks
############################

### drop treatment cases w/o usage data (same dataset as obs study and mediation analysis)
sdatDrop <- makeStanDat(dat,x,missingUsage=FALSE)
psmodDrop <- stan('R/prinStratStan.stan',data=sdatDrop)
save(psmodDrop,file='output/psmodDrop.RData')

mainb1 <- summary(psmod1,par='b1')[[1]]['b1',]
robustness <- data.frame(Value=mainb1['mean'],Coefficient='Main Model',HighInner=mainb1['75%'],LowInner=mainb1['25%'],
                                      HighOuter=mainb1['97.5%'],LowOuter=mainb1['2.5%'],Model='model')
ests <- c(mainb1['mean'],mainb1['sd'])
modFiles <- c('xInteractions','modNoTeacher','pooledU','stanMod2pl','stanMod3pl','bcModel','hardSections','fullMod','rawMod','rawModBC2')

for(mmm in modFiles){
    modName <- load(paste0('output/',mmm,'.RData'))
    modName <- modName[1]
    modSumm <- summary(get(modName),par='b1')[[1]]['b1',]
    ests <- rbind(ests,c(modSumm['mean'],modSumm['sd']))
    robustness <- rbind(robustness,
                        data.frame(Value=modSumm['mean'],Coefficient=modName,HighInner=modSumm['75%'],
                                   LowInner=modSumm['25%'],HighOuter=modSumm['97.5%'],LowOuter=modSumm['2.5%'],
                                   Model='model'))
    rm(list=modName); gc()
}

#robustness <- subset(robustness,Coefficient!='rawMod')

modelNames <- c('Main Model','Covariate Interations','No Teacher Effects','Pooled Rasch Model','2PL Mastery', '3PL Mastery','Power Transform Y','Exclude Sections w/o Hints','Include All Data','Raw Scores',
                'Raw Scores (BC Trans)')
robustness$Coefficient <- factor(modelNames,levels=rev(modelNames))

coefplot.data.frame(subset(robustness,modelNames!='Raw Scores'),xlab=expression(hat(b)[1]),ylab=NULL,lwdOuter=0.5,lwdInner=1.5)+theme(text=element_text(size=15))

rob2 <- subset(robustness,modelNames!='Raw Scores')
rob2$estX <- min(rob2$LowOuter)-0.05
rob2$seX <- min(rob2$LowOuter)-0.02
rob2$est <- paste0('  ',round(ests[-which(robustness$Coefficient=='Raw Scores'),'mean'],3),' (',round(ests[-which(robustness$Coefficient=='Raw Scores'),'sd'],3),')')
#rob2$se <- paste0(

ggplot(rob2,aes(x=Value,y=Coefficient))+geom_point(size=2)+
    geom_errorbarh(aes(xmin=LowInner,xmax=HighInner),size=1.5,height=0)+
    geom_errorbarh(aes(xmin=LowOuter,xmax=HighOuter),size=0.5,height=0)+
    geom_text(aes(x=estX,label=est),hjust='left')+
                                        #geom_text(aes(x=seX,label=se))+
    theme_minimal()+
    ylab('')+xlab(expression(b[1]))+scale_x_continuous(breaks=c(-0.1,-0.05,0,0.05))+geom_vline(xintercept=0,linetype='dotted')+theme(text=element_text(size=15))

ggsave('output/robustness.pdf',height=8,width=6)



##################################
### fake data models
################################
