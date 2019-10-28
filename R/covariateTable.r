## covariate means/balance

covs2 <- covs[rownames(dat),] ### same as dat but with NAs

miss <- NULL
varbs <- c('grade','esl','frl','race','sex','spec','xirt')
for(i in varbs)
  miss <- rbind(miss,
    c(sum(is.na(covs2[[i]])),mean(is.na(covs2[[i]])),error[i,'error']))
miss <- as.data.frame(miss)
#miss$`Error Type` <- c('PFC','PFC','PFC','SRMSE')
#rownames(miss) <- c('Race/Ethnicity','Sex','Special Education','Pretest')
rownames(miss) <- varbs
names(miss)[1:3] <- c('# Missing','% Missing','Imputation Error')
miss[,2] <- as.integer(round(miss[,2]*100))
miss[,1] <- as.integer(miss[,1])
miss['xirt','Imputation Error'] <- sqrt(miss['xirt','Imputation Error'])/sd(covs2$xirt,na.rm=TRUE)




balTest <- balanceTest(treatment~grade+esl+frl+race+sex+spec+xirt+cluster(schoolid2)+strata(pair),data=dat,report=c('adj.means','std.diffs','chisquare.test'),include.NA.flags=T)

## use "unstrat" for group means cuz no one knows/cares about "adjusted means"
trtMeans <- round(balTest$results[,'Treatment','Unstrat']*100)
ctlMeans <- round(balTest$results[,'Control','Unstrat']*100)
stdDiff <- sprintf('%.2f',balTest$results[,'std.diff','pair'])
names(stdDiff) <- dimnames(balTest$results)[[1]]

impErr <- sprintf('%.2f',miss[,3])
names(impErr) <- rownames(miss)

results <- list()
results$balP <- sprintf('%.2f',balTest$overall['pair','p.value'])

varNames <- c(
  grade='Grade',
  esl='ELL',
  frl='FRL',
  race='Ethnicity',
  sex='Sex',
  spec='Sp. Ed.',
  xirt='Pretest'
)

levNames <- c(
  `9`='9th',
  higher='$>$9th',
  `0`='No',
  `1`='Yes',
  WhiteAsian='White/Asian',
  BlackMulti='Black/Multi',
  HispAIAN='Hispanic/Nat.Am.',
  F='Female',
  M='Male',
  gifted='Gifted',
  speced='Sp. Ed.',
  typical='Typical'
)

nlev <- sapply(dat[,varbs],nlevels)

cat0 <- function(...) cat(...,sep='')

sink('output/covariateTable.tex')
cat0('
\\begin{tabular}{cccrccc}
  \\hline
  \\hline
 &\\% Miss.& Imp. Err.&Levels& BaU & CTA1 & Std. Diff.\\\\
  \\hline
  \\hline
')
for(vv in varbs){
  if(nlev[vv]){
    cat0('\\multirow{',nlev[vv],'}{*}{',varNames[vv],'}&\\multirow{',nlev[vv],'}{*}{',miss[vv,2],'\\%} &\\multirow{',nlev[vv],'}{*}{',impErr[vv],'}&')
    levs <- levels(dat[[vv]])
    for(ll in 1:nlev[vv]){
      if(ll> 1) cat0('&&&')
      cat0(levNames[levs[ll]],'&',
      ctlMeans[paste0(vv,levs[ll])],'\\% &',
      trtMeans[paste0(vv,levs[ll])],'\\% &',
      stdDiff[paste0(vv,levs[ll])],'\\\\\n')
    }
  } else cat0(
    varNames[vv],'&',miss[vv,2],'\\%&',impErr[vv],'& &',
    sprintf('%.2f',balTest$results[vv,'Control','Unstrat']),'& ',
    sprintf('%.2f',balTest$results[vv,'Treatment','Unstrat']),'& ',stdDiff[vv],'\\\\\n')
  cat0('\\hline\n')
}
cat0('&&&\\multicolumn{4}{c}{Overall Covariate Balance: p=',results$balP,'}\\\\
\\hline
\\hline
\\end{tabular}\n')
sink()
