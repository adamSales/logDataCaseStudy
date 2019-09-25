library(tikzDevice)
### run all four fake models

source('R/fakePrep.r')
source('R/fakeNoEff.r')
source('R/fakeConstEff.r')
source('R/fakeLinEff.r')
source('R/fakeQuadEff.r')

extract <- rstan::extract

######################
### results from fake models
#####################
estEff <- list()
print(load('fitModels/noEffect.RData'))
pd <- pdMod(noEff,1,1,function(x) x*0)
estEff$noEff <- summary(noEff,par=c('b0','b1'))[[1]]
rm(noEff);gc()

print(load('fitModels/constEff.RData'))
pd <- rbind(pd,pdMod(constEff,1,2,function(x) x*0+0.18))
estEff$constEff <- summary(constEff,par=c('b0','b1'))[[1]]
rm(constEff); gc()

print(load('fitModels/linEff.RData'))
pd <- rbind(pd,pdMod(linEff,2,1,function(x) mean(effs$b0)+x*mean(effs$b1)))
estEff$linEff <- summary(linEff,par=c('b0','b1'))[[1]]
rm(linEff);gc()

print(load('fitModels/quadEff.RData'))
U <- datF$U
mux <- mean(U)
te <- -(U-mux)^2
sigte <- sd(te)
mute <- mean(te/sigte*0.1)


pd <- rbind(pd,pdMod(quadEff,2,2,function(x) -(x-mux)^2*0.1/sigte-mute+0.13))
estEff$quadEff <- summary(quadEff,par=c('b0','b1'))[[1]]
rm(quadEff);gc()

pd$title <- NA
pd <- within(pd,{
    title[row==1 & column==1] <- paste0('$\\tau=0$\n$\\hat{\\tau}=',sprintf("%.2f",estEff$noEff['b0',1]),
                                        ifelse(estEff$noEff['b1',1]>0,'+',''),
                                        sprintf("%.2f",estEff$noEff['b1',1]),'\\eta_T$')
    title[row==1 & column==2] <- paste0('$\\tau=0.18+\\epsilon$\n$\\hat{\\tau}=',
                                        sprintf("%.2f",estEff$constEff['b0',1]),
                                        ifelse(estEff$constEff['b1',1]>0,'+',''),
                                        sprintf("%.2f",estEff$constEff['b1',1]),'\\eta_T$')
    title[row==2 & column==1] <- paste0('$\\tau=',
                                        round(mean(effs$b0),2),'+',
                                        round(mean(effs$b1),2),
                                        '\\eta_T$\n$\\hat{\\tau}=',
                                        sprintf("%.2f",estEff$linEff['b0',1]),
                                        ifelse(estEff$linEff['b1',1]>0,'+',''),
                                        (sprintf("%.2f",estEff$linEff['b1',1])),'\\eta_T$')
    title[row==2 & column==2] <- paste0('$\\tau=-0.04\\eta_T^2+0.01\\eta_T+0.02$\n$\\hat{\\tau}=',
                                        sprintf("%.2f",estEff$quadEff['b0',1]),
                                        ifelse(estEff$quadEff['b1',1]>0,'+',' '),
                                        sprintf("%.2f",estEff$quadEff['b1',1]),'\\eta_T$')})

pd <- within(pd, {
    title <- factor(title,levels=c(title[row==1 & column==1][1],
                                   title[row==1 & column==2][1],
                                   title[row==2 & column==1][1],
                                   title[row==2 & column==2][1]))})

tikz('output/fakePlots.tex',standAlone=TRUE,width=6,height=6)
print(ggplot(pd)+
    geom_abline(aes(intercept=b0,slope=b1,group=id),color='red')+
    coord_cartesian(xlim=c(min(pd$xmin),max(pd$xmax)),ylim=c(min(pd$ymin),max(pd$ymax)),expand=FALSE)+
    geom_line(aes(x=x,y=y,group=truthOrAvg,linetype=truthOrAvg,color=truthOrAvg,alpha=truthOrAvg),size=1.5)+
    facet_wrap(~title,ncol=2)+xlab('$\\eta_T$')+ylab('$\\hat{\\tau}(\\eta_T)$')+
    labs(group=NULL,color=NULL,linetype=NULL)+
    #theme(strip.background = element_blank(),strip.text.x = element_blank(),strip.text.y=element_blank())+
    scale_color_manual(values=c('black','red','black'))+scale_linetype_manual(values=c('solid','solid','dotted'))+scale_alpha_manual(values=c(1,0,1),guide=FALSE)+theme(legend.position='top')+theme(text=element_text(size=15))+theme(legend.key.width=unit(.5,'in')))
    dev.off()
setwd('output'); tools::texi2dvi('fakePlots.tex', pdf = T, clean = T); setwd('..')

