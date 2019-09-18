
#### record of some of the process of coming up with a match

psmod0 <- lmer(qlogis(propHint)~(1|state)+(1|schoolid2)+(1|classid2)+sex+grade+race+ns(xirt,5)+spec+esl+gradeMIS+raceMIS+sexMIS+frlMIS,data=subset(dat1,propHint<1&propHint>0))

psmod3.2 <- ## stan_
  glmer(M~(1|state)+(1|schoolid2)+(1|classid2)+sex+grade+race+ns(xirt,5)+spec+esl+gradeMIS+raceMIS+sexMIS+frlMIS,family=binomial,data=dat1)

psmod1 <-   glmer(M~(1|state)+(1|schoolid2)+(1|classid2)+grade+race+xirt+spec+esl,family=binomial,data=dat1)
#save(psmod3.2,file='stanPSmod2.RData')
#load('stanPSmod.RData')


dist1 <- match_on(M~predict(psmod3.2,type='link'),data=dat1)
m1 <- fullmatch(dist1,data=dat1)

dist2 <- match_on(M~predict(psmod3.2,type='link')+xirt,data=dat1)
m2 <- fullmatch(dist2,data=dat1)
balMod <- xBalance(M~grade+race+xirt+spec+esl,strata=list(unmatched=NULL,m2=~m2),report=c('std.diff','z.score','chisquare.test'),data=droplevels(dat1))

dist3 <- match_on(M~predict(psmod3.2,type='link'),within=exactMatch(M~schoolid2,data=dat1),data=dat1)
m3 <- fullmatch(dist3,data=dat1)

balMod <- xBalance(M~grade+race+xirt+spec+esl+schoolid2+state,strata=list(unmatched=NULL,m2=~m2,m3=~m3),report=c('std.diff','z.score','chisquare.test'),data=droplevels(dat1))

plot(balMod,colors=subwayPalette)
abline(v=c(-.05,-.25,.05,.25),lty='dotted')


dist4.3 <- match_on(M~predict(psmod3.2,type='link'),within=exactMatch(M~schoolid2,data=dat1),data=dat1)
m4.3 <- fullmatch(dist4.3,data=dat1)
balMod <- xBalance(M~grade+race+xirt+spec+esl,strata=list(unmatched=NULL,m4.3=~m4.3),report=c('std.diff','z.score','chisquare.test'),data=droplevels(dat1))

plot(balMod,colors=subwayPalette)
abline(v=c(-.05,-.25,.05,.25),lty='dotted')

psmod3.3 <- arm::bayesglm(M~grade+race+poly(xirt,2)+spec+esl,family=binomial,data=dat1)
dist4 <- match_on(M~psmod3.3$linear,data=dat1)
m4 <- fullmatch(dist4,data=dat1)
balMod <- xBalance(M~grade+race+xirt+spec+esl,strata=list(unmatched=NULL,m4=~m4),report=c('std.diff','z.score','chisquare.test'),data=droplevels(dat1))

dist5 <- match_on(M~psmod3.3$linear+xirt,data=dat1)
m5 <- fullmatch(dist5,data=dat1)
balMod <- xBalance(M~grade+race+xirt+spec+esl,strata=list(unmatched=NULL,m4=~m4,m5=~m5),report=c('std.diff','z.score','chisquare.test'),data=droplevels(dat1))


dist6 <- match_on(M~predict(psmod3.2,type='link')+xirt,within=exactMatch(M~schoolid2,data=dat1),data=dat1)
m6 <- fullmatch(dist6,data=dat1)

balMod <- xBalance(M~grade+race+xirt+spec+esl,strata=list(unmatched=NULL,m5=~m5,m6=~m6),report=c('std.diff','z.score','chisquare.test'),data=droplevels(dat1))

psmod4 <-  arm::bayesglm(M~grade+race+splines::ns(xirt,5)+spec+esl,family=binomial,data=dat1)
dist7 <- match_on(M~psmod4$linear+xirt,within=exactMatch(M~schoolid2,data=dat1),data=dat1)
m7 <- fullmatch(dist7,data=dat1)

balMod <- xBalance(M~grade+race+xirt+spec+esl,strata=list(unmatched=NULL,m3=~m3,m6=~m6),report=c('std.diff','z.score','chisquare.test'),data=droplevels(dat1))





dist4.4 <- match_on(M~psmod3.3$linear+xirt+strata(pair),data=dat1)#,caliper=0.1)
m4.4 <- fullmatch(dist4.4,data=dat1)
balMod <- xBalance(M~grade+race+xirt+spec+esl,strata=list(unmatched=NULL,m1=~m1,m4.4=~m4.4,m4.3=~m4.3),report=c('std.diff','z.score','chisquare.test'),data=droplevels(dat1))
