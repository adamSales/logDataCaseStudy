## covariate means/balance
miss <- NULL
for(i in c('race','sex','spec','xirt')) miss <- rbind(miss,
 c(sum(is.na(covs[[i]])),mean(is.na(covs[[i]])),error[i,'error']))
miss <- as.data.frame(miss)
miss$`Error Type` <- c('PFC','PFC','PFC','SRMSE')
rownames(miss) <- c('Race/Ethnicity','Sex','Special Education','Pretest')
names(miss)[1:3] <- c('# Missing','% Missing','Imputation Error')
miss[,2] <- as.integer(round(miss[,2]*100))
miss[,1] <- as.integer(miss[,1])
miss['Pretest','Imputation Error'] <- sqrt(miss['Pretest','Imputation Error'])/sd(covs$xirt,na.rm=TRUE)


covs2 <- covs[rownames(dat),] ### same as dat but with NAs

balTest <- balanceTest(treatment~race+sex+spec+xirt+cluster(schoolid2)+strata(pair),data=dat,report=c('adj.means','std.diffs','chisquare.test'),include.NA.flags=T)

## use "unstrat" for group means cuz no one knows/cares about "adjusted means"
trtMeans <- round(balTest$results[,'Treatment','Unstrat']*100)
ctlMeans <- round(balTest$results[,'Control','Unstrat']*100)
stdDiff <- sprintf('%.2f',balTest$results[,'std.diff','pair'])
names(stdDiff) <- dimnames(balTest$results)[[1]]

impErr <- sprintf('%.2f',miss[,3])
names(impErr) <- rownames(miss)

results <- list()
results$balP <- sprintf('%.2f',balTest$overall['pair','p.value'])


cat('
\\begin{tabular}{cccrccc}
  \\hline
  \\hline
 &\\% Miss.& Imp. Err.&Levels& BaU & CTA1 & Std. Diff.\\\\
  \\hline
  \\hline
\\multirow{3}{*}{Ethnicity}&\\multirow{3}{*}{',miss['Race/Ethnicity',2],'\\%}
  &\\multirow{3}{*}{',impErr['Race/Ethnicity'],'}& White/Asian & ',
ctlMeans['raceWhiteAsian'],'\\% &',trtMeans['raceWhiteAsian'],'\\% & ',stdDiff['raceWhiteAsian'],'\\\\
&&&Black/Multi &',
ctlMeans['raceBlackMulti'],'\\% &',trtMeans['raceBlackMulti'],'\\% & ',stdDiff['raceBlackMulti'],'\\\\
 & &&Hispanic/Nat.Am. & ',
ctlMeans['raceHispAIAN'],'\\% &',trtMeans['raceHispAIAN'],'\\% & ',stdDiff['raceHispAIAN'],' \\\\
\\hline
\\multirow{2}{*}{Sex}&\\multirow{2}{*}{',miss['Sex',2],'\\%}&\\multirow{2}{*}{',impErr['Sex'],'}& Female & ',ctlMeans['sexF'],'\\% &',trtMeans['sexF'],'\\% & ',stdDiff['sexF'],'\\\\
  &&&Male &',
ctlMeans['sexM'],'\\% &',trtMeans['sexM'],'\\% & ',stdDiff['sexM'],'\\\\
\\hline
\\multirow{3}{*}{Sp. Ed.}&\\multirow{3}{*}{',miss['Special Education',2],'\\%}&\\multirow{3}{*}{',impErr['Special Education'],'}&Typical &',
ctlMeans['spectypical'],'\\% &',trtMeans['spectypical'],'\\% & ',stdDiff['spectypical'],'\\\\
&&&Spec. Ed & ',
ctlMeans['specspeced'],'\\% &',trtMeans['specspeced'],'\\% & ',stdDiff['specspeced'],' \\\\
&&&Gifted &',
ctlMeans['specgifted'],'\\% &',trtMeans['specgifted'],'\\% & ',stdDiff['specgifted'],'\\\\
\\hline
Pretest&',miss['Pretest',2],'\\%&',impErr['Pretest'],'&   &',
sprintf('%.2f',balTest$results['xirt','Control','Unstrat']),'& ',
sprintf('%.2f',balTest$results['xirt','Treatment','Unstrat']),'& ',stdDiff['xirt'],'\\\\
   \\hline
&&&\\multicolumn{4}{c}{Overall Covariate Balance: p=',results$balP,'}\\\\
\\hline
\\hline
\\end{tabular}\n',sep='',file='output/covariateTable.tex')
