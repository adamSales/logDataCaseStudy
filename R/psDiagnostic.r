
source('R/diagnosticPlots.r')

source('R/robustnessChecks.r')

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
