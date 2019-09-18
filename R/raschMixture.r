
raschMix <- stan('R/raschMix.stan',data=sdat2,chains=5,iter=2000)
save(raschMix,file='output/mixtureModel2.RData')
ss <- summary(raschMix)
studEffSumm <- ss$summary[grep('studEff',rownames(ss$summary)),]
mean(studEffSumm[,'Rhat']<1.01)
draws <- rstan::extract(raschMix,permute=TRUE)

rhat <- function(mat){
  n <- nrow(mat)
  W <- apply(mat,2,var)
  W <- mean(W)
  B <- var(colMeans(mat))
  sqrt(((1-1/n)*W+B)/W)
}


post <- function(theta,p1,mu0,mu1,sigma){
  p0 <- 1-p1
  post0 <- dnorm(theta,mu0,sigma[1])*p0
  post1 <- dnorm(theta,mu1,sigma[2])*p1
  post1/(post0+post1)
}

classProbs <- with(draws,
  vapply(1:nrow(studEff),
    function(i)
      post(
        theta=studEff[i,],
        p1=p1[i],
        mu0=mu0[i],
        mu1=mu1[i],
        sigma=sigma[i,]
      ),
    numeric(ncol(studEff))
  )
)

classProbsOverall <- rowMeans(classProbs)


hintAvgStud$classProb <- classProbsOverall[1:nrow(hintAvgStud)]
hintAvgStud$studEff <- colMeans(draws$studEff[,1:nrow(hintAvgStud)])

## mm <- glmer(hint~(1|field_id)+(1|section),family=binomial,data=x)
## save(mm,file='hardHintRaschNew.RData')

## re <- ranef(mm)
## bbb <- x%>%group_by(field_id)%>%summarize(nprob=n(),hint=mean(nhints1>0),err=mean(nerrs1>0), hintHard=mean(nhints1[nerrs1>0]>0),nhard=sum(nerrs1>0|nhints1>0))

## bbb$theta <- re[[1]][as.character(bbb$field_id),1]
## bbb$theta <- bbb$theta+fixef(mm)[1]
## med <- median(bbb$theta,na.rm=TRUE)
## bbb$M <- bbb$theta>med
## dat <- merge(dat,bbb,all.x=TRUE)

save(hintAvgStud,file='output/raschMixOutput.RData')
dat1 <- inner_join(dat,select(hintAvgStud,field_id,classProb,studEff,propHintStud))

#dat1$theta <- colMeans(draws$studEff)
#dat1$prob <- classProbsOverall
dat1 <- mutate(dat1,
  M0=classProb> quantile(classProb,1-mean(draws$p1)),
  M1=propHintStud>quantile(propHintStud,1-mean(draws$p1)),
  M2=propHintStud>quantile(propHintStud,2/3),
  M=propHintStud>.6)

