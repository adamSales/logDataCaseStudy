#### mediation?

### estimate E[Y(1,M(0))]
dat10 <- subset(dat1,!M)
tab <- xtabs(~match+M,data=dat1)
K <- 1+tab[,2]/tab[,1]
dat10$ww <- K[dat10$match]
out5 <- lm(Y~1,data=dat10,weights=ww)


### SE depends on question. for
vvv <- vcovCR(out5,dat10$schoolid2,'CR2')
se1 <- sqrt(vvv[1,1])


out5 <- lmer(Y~(1|schoolid2)+(1|classid2),data=dat10)
muhat <- weighted.mean(dat10$Y,dat10$ww)
cls <- dat10%>%group_by(classid2)%>%summarize(ww=sum(ww)^2)
scl <- dat10%>%group_by(schoolid2)%>%summarize(ww=sum(ww)^2)
vvv <- as.data.frame(VarCorr(out5))
varMu <- (vvv$vcov[1]*sum(cls$ww)+vvv$vcov[2]*sum(scl$ww)+vvv$vcov[3]*sum(dat10$ww^2))/nrow(dat1)^2
se2 <- sqrt(varMu)


### direct effect
ddd1 <- dat10
ddd1$Y <- with(ddd1,Y*ww/sum(ww)*nrow(ddd1))
ddd1$xirt <- with(ddd1,xirt*ww/sum(ww)*nrow(ddd1))
datDir1 <- rbind(subset(dat,treatment==0),ddd1[,names(dat)])
dir1 <- lm(Y~treatment+pair,data=datDir1)
vvv <- vcovCR(dir1,datDir1$schoolid2,'CR2')
se1 <- sqrt(vvv[1,1])

dir2 <- lmer(Y~treatment+pair+(1|classid2)+(1|schoolid2),data=datDir1)
dir3 <- lmer(Y~treatment+poly(xirt,2)+race+grade+esl+frl+sex+pair+(1|classid2)+(1|schoolid2),data=datDir1)

### mediated effect
### E[Y(1,m(1))-Y(1,m(0))]
ddd2 <- dat1
adj <- lm(Y~M+schoolid2+poly(xirt,2)+race+grade+esl+frl+sex,data=dat1)
ddd2$Ytilde <- ddd2$Y-model.matrix(adj)[,-c(1:2)]%*%coef(adj)[-c(1:2)]
ddd2$ww <- with(ddd2,1-(1-M)*K[dat1$match])
ddd2$Yw <- with(ddd2,Y*ww)
ddd2$Ytw <- with(ddd2,Ytilde*ww)
med1 <- lm(Yw~1,data=ddd2)
med2 <- lm(Ytw~1,data=ddd2)

### alternatively, remove what's outside of common support
datC <- subset(dat,treatment==0)
dc <- model.matrix(~grade+race+poly(xirt,2)+spec+esl,datC)
ps0 <- dc%*%fixef(psmod3.2)

ps0 <- predict(psmod3.3,datC,type='link')
datC$theta <- NA
datC$M2 <- NA
datC$match <- NA
datC$ps <- ps0
dat1$ps <- psmod3.3$linear
dat2 <- rbind(datC[,names(dat1)],dat1)

distC <- match_on(treatment~ps+xirt+strata(pair),data=dat2)
distC <- as.matrix(distC)
distC <- distC[!is.na(m4.4),]
distC[distC>max(dmat[is.finite(dmat)])] <- Inf
keep <- which(apply(distC,2,function(x) any(is.finite(x))))


dummy <- fullmatch(distC,data=dat2)

datC <- subset(datC,!is.na(dummy[dat2$treatment==0]))
dat2.2 <- rbind(datC[,names(dat1)],dat1[!is.na(m4.4),])

dat2.2$qq <- ifelse(dat2.2$treatment==1,
                    ifelse(dat2.2$M,1,2),0)
