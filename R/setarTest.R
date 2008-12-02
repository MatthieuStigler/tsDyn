
###tests
#SSR of linear
 # check for stationarty
#call selectSETAR: obtain 2 best thresh


#buildXth1Common / buildXth1NoCommon / buildXth2Common / buildXth2NoCommon
 #   a<-lm()
 # is.InUnitCircle(a) to check
 # deviance(a) for SSR

#compute/return stat

###bootstrap
#bootstrap data
#selectSETAR()

#SSR of those models: SSR_1thresh / SSR_2threshCommon / SSR_2threshNoCommon

#replicate()


setarTest <- function (x, m, d = 1, steps = d, series, thDelay = 0, nboot=10, plot=FALSE, trim=0.1, test=c("1vs", "2vs3"), check=FALSE)
{

include<-"const" #other types not implemented in setar.sim

###setarTest 1: SSR and check linear model
  linear<-linear(x, m, d = 1, steps = d, series)
  SSR<-deviance(linear)
  
 
###setarTest 2: search best thresholds: call selectSETAR
  search<-selectSETAR(x, m=m,d=1, steps=d, thDelay=thDelay, trace=FALSE, include =include, common=FALSE, model=c("TAR", "MTAR"),nthresh=2,trim=trim,criterion = c("SSR"),thSteps = 7,ngrid="ALL",  plot=FALSE,max.iter=3) 
  
  firstBests<-search$firstBests
  bests<-search$bests
  

  
  ### Obtain infos for the two thresh models
  set1<-setar(x, m, d=d, steps=steps, thDelay=0, th=firstBests["th"], trace=FALSE, nested=FALSE,include = c("const", "trend","none", "both"), common=FALSE, model=c("TAR", "MTAR"), nthresh=1,trim=trim)
  
   set2<-setar(x, m, d=d, steps=steps, thDelay=0, th=bests[c("th1", "th2")], trace=FALSE, nested=FALSE,include = c("const", "trend","none", "both"), common=FALSE, model=c("TAR", "MTAR"), nthresh=2,trim=trim)
  #RESTRICTED currently:
  #common: FALSE
  #nthresh=2
  #max.iter=1
 n<-length(na.omit(residuals(set2)))
  ###setarTest 10: SSR of TAR(+) and (2),compute F stat and print
  SSR1thresh<-search$firstBests["SSR"]
  SSR2thresh<-search$bests["SSR"]
  SSRs<-matrix(c(SSR, SSR1thresh, SSR2thresh), ncol=3, dimnames=list("SSR", c("AR", "TAR(1)", "TAR(2)")))
  
  ###F test for original data
  Ftest12<-as.numeric(n*(SSR-SSR1thresh)/SSR1thresh)
  Ftest13<-as.numeric(n*(SSR-SSR2thresh)/SSR2thresh)
  Ftest23<-as.numeric(n*(SSR1thresh-SSR2thresh)/SSR2thresh)
  Ftests<-matrix(c(Ftest12, Ftest13, Ftest23),ncol=3, dimnames=list("Ftest", c("1vs2", "1vs3", "2vs3")))
  
########
###Boot
########
bootHoAr<-function(linearObject, type){
  ### Bootstrap under null: ar model
  bootLin<-setar.sim(setarObject=linearObject, type=type)$serie
  
###setarTest 1: SSR and check linear model
  linearBoot<-linear(bootLin, m=m)
  SSR<-deviance(linearBoot)
   
###setarTest 2: search best thresholds: call selectSETAR
  searchBoot<-selectSETAR(bootLin, m=m,d=1, steps=d, thDelay=thDelay, trace=FALSE, include =include, common=FALSE, model=c("TAR", "MTAR"),nthresh=2,trim=trim,criterion = c("SSR"),thSteps = 7,ngrid="ALL",  plot=FALSE,max.iter=3) 
  
  firstBests<-searchBoot$firstBests
  bests<-searchBoot$bests
  
###setarTest 10: SSR of all boot models
  SSR1thresh<-searchBoot$firstBests["SSR"]
  SSR2thresh<-searchBoot$bests["SSR"]
    
###F test for boot data
  Ftest12<-as.numeric(n*(SSR-SSR1thresh)/SSR1thresh)
  Ftest13<-as.numeric(n*(SSR-SSR2thresh)/SSR2thresh)
  return(c(Ftest12, Ftest13))
}

bootHoSetar1<-function(setarObject, type ){
  ### Bootstrap under null: Setar(1t) model
  bootSet<-setar.sim(setarObject=set1, type=type, nthresh=1)$serie
  if(check)
    print(all(bootSet-set1$str$x<0.00005))


###setarTest 2: search best thresholds: call selectSETAR
  searchBoot<-selectSETAR(bootSet, m=m,d=1, steps=d, thDelay=thDelay, trace=FALSE, include =include, common=FALSE, model=c("TAR", "MTAR"),nthresh=2,trim=trim,criterion = c("SSR"),thSteps = 7,ngrid="ALL",  plot=FALSE,max.iter=3) 
  
  firstBests<-searchBoot$firstBests
  bests<-searchBoot$bests
  
###setarTest 10: SSR of all boot models
  SSR1thresh<-searchBoot$firstBests["SSR"]
  SSR2thresh<-searchBoot$bests["SSR"]
    
###F test for boot data
  Ftest23<-as.numeric(n*(SSR1thresh-SSR2thresh)/SSR2thresh)
  

return(Ftest23)
}


### Run the function (boot, search best SSR for lin, set 1 and set2 ) and extract results


type<-ifelse(check, "check", "boot")

probs<-c(0.9, 0.95, 0.975,0.99)
if(test=="1vs"){
  Ftestboot<-replicate(nboot, bootHoAr(linear, type))
  Ftestboot12<-Ftestboot[1,]
  Ftestboot13<-Ftestboot[2,]
  PvalBoot12<-mean(ifelse(Ftestboot12>Ftest12,1,0))
  CriticalValBoot12<-quantile(Ftestboot12, probs=probs)
  PvalBoot13<-mean(ifelse(Ftestboot13>Ftest13,1,0))
  CriticalValBoot13<-quantile(Ftestboot13, probs=probs)
  CriticalValBoot<-matrix(c(CriticalValBoot12,CriticalValBoot13), nrow=2, dimnames=list(c("1vs2", "1vs3"), probs))
  PvalBoot<-c(PvalBoot12,PvalBoot13)
}
else{
  Ftestboot<-replicate(nboot, bootHoSetar1(set1, type))
  Ftestboot23<-Ftestboot
  print(Ftestboot23)
  PvalBoot23<-mean(ifelse(Ftestboot23>Ftest23,1,0))
  CriticalValBoot23<-quantile(Ftestboot23, probs=probs)
  CriticalValBoot<-matrix(CriticalValBoot23, nrow=1, dimnames=list("2vs3", probs))
  PvalBoot<-PvalBoot23
}


###res
  res<-list(Ftests=Ftests, SSRs=SSRs, firstBests=search$firstBests["th"],secBests=search$bests[c("th1", "th2")], CriticalValBoot=CriticalValBoot,PvalBoot=PvalBoot, test=test, m=m, Ftestboot=Ftestboot, nboot=nboot)
  class(res)<-"Hansen99Test"
  return(res)
}
  
print.Hansen99Test<-function(x,...){
  if(x$test=="1vs"){
    cat("Test of linearity against setar(2) and setar(3)\n\n")
    print(matrix(c(x$Ftests[-3], x$PvalBoot), ncol=2, dimnames=list(c("1vs2", "1vs3"), c("Test", "Pval"))))
  }
   else{
    cat("Test of setar(2) against  setar(3)\n\n")
    print(matrix(c(x$Ftests[3], x$PvalBoot), ncol=2, dimnames=list(c("2vs3"), c("Test", "Pval"))))
   }
}

summary.Hansen99Test<-function(x, ...){
  print.Hansen99Test(x)
  cat("\nCritical values:\n")
  print(x$CriticalValBoot)
  cat("\nSSR of original series:\n")
  print(matrix(x$SSRs, ncol=1, dimnames=list(c("AR", "SETAR(2)", "SETAR(3)"), "SSR")))
  cat("\nThreshold of original series:\n")
  print(matrix(c(x$firstBests, NA, x$secBests),byrow=TRUE, ncol=2, dimnames=list(c("SETAR(2)", "SETAR(3)"), c("th1", "th2"))))
  cat("\nNumber of bootstrap replications: ", x$nboot, "\n")
}

plot.Hansen99Test<-function(x,...){
  m<-x$m
  if(x$test=="1vs"){
    layout(c(1,2))
    Ftestboot12<-x$Ftestboot[1,]
    Ftestboot13<-x$Ftestboot[2,]
    Ftest12<-x$Ftests[1]
    Ftest13<-x$Ftests[2]

    #Plot 1vs2
    plot(density(Ftestboot12, from=0), xlab="Ftest12", xlim=c(0,max(Ftest12+1,max(Ftestboot12))),ylim=c(0,max(density(Ftestboot12)$y,dchisq(0:Ftest12, df=1+m))), main="")
    title("Test linear AR vs 1 threshold SETAR")
    abline(v=Ftest12, lty=2, col=2)
    curve(dchisq(x, df=1+m, ncp=0), from=0, n=Ftest12+5, add=TRUE, col=3)
    legend("topright", legend=c("Asymptotic Chi 2", "Bootstrap", "Test value"), col=c(3,1,2), lty=c(1,1,2))
  
    #Plot 1vs3
    plot(density(Ftestboot13, from=0), xlab="Ftest13", xlim=c(0,max(Ftest13+1,max(Ftestboot13))),ylim=c(0,max(density(Ftestboot13)$y,dchisq(0:Ftest12, df=2*(1+m)))), main="")
    title("Test linear AR vs 2 thresholds SETAR")
    abline(v=Ftest13, lty=2, col=2)
    cat("ok\n")
    curve(dchisq(x, df=2*(1+m), ncp=0), from=0, n=Ftest13+5, add=TRUE, col=3)
    legend("topright", legend=c("Asymptotic Chi 2", "Bootstrap", "Test value"), col=c(3,1,2), lty=c(1,1,2))
  }
  #plot 2vs3
  else {
    Ftestboot23<-x$Ftestboot
    Ftest23<-x$Ftests[3]
    
    plot(density(Ftestboot23, from=0), xlab="Ftest23", xlim=c(0,max(Ftest23+1,Ftestboot23)), ylim=c(0,max(density(Ftestboot23)$y)), main="")
    title("Test 1 threshold SETAR vs 2 thresholds SETAR")
    abline(v=Ftest23, lty=2, col=2)
    curve(dchisq(x, df=1+m, ncp=0), from=0, n=Ftest23+5, add=TRUE, col=3)
    legend("topright", legend=c("Asymptotic Chi 2", "Bootstrap", "Test value"), col=c(3,1,2), lty=c(1,1,2))
  }
}


if(FALSE){  
library(tsDyn)
sun<-(sqrt(sunspot.year+1)-1)*2

###Sunsport test
#Test 1vs2 and 1vs3
environment(setarTest)<-environment(star)
Han1<-setarTest(sun, m=11, thDelay=0:1, nboot=2, plot=TRUE, trim=0.1, test="1vs")

print(Han1)
summary(Han1)
plot(Han1)

#Test 2vs3
environment(setarTest)<-environment(star)
Han2<-setarTest(sun, m=11, thDelay=0:1, nboot=10, plot=TRUE, trim=0.1, test="2vs3")

print(Han2)
summary(Han2)
plot(Han2)

###two probs
#- deviance of selectSETAR and setar is nto the same
# selectSETAR does not select the good one...
environment(selectSETAR)<-environment(star)
selectSETAR(lynx, m=11, thDelay=1, nthresh=2, th=list(exact=c(5.3,8)), criterion="SSR")
selectSETAR(sun, m=11, thDelay=1, nthresh=2, th=list(exact=c(5.3,8)), criterion="SSR")

###US IP
IP<-read.table(file="/media/sda5/Mes documents/Uni/Statistique/Time Series/Handbooks/datasets/Hansen/Hansen1999Linearity/GAUSS/IP.DAT")
IPts<-ts(IP, start=c(1960,1), freq=12)
IPts

IP<-read.table(file="/media/sda5/Mes documents/Uni/Statistique/Time Series/Handbooks/datasets/Hansen/Hansen1999Linearity/Matlab/ipdat.txt")
IPts<-ts(IP, start=c(1960,1), freq=12)
IPts
dat<-log(IPts)
dat2<-(dat[13:length(dat)]-dat[1:(length(dat)-12)])*100
dat2bis<-diff(dat, 12)*100
cbind(dat2, dat2bis)
# dat=(dat(13:length(dat(:,1)))-dat(1:length(dat(:,1))-12))*100
dat3<-dat2[157:length(dat)]
# dat=dat(157:length(dat(:,1)));
dat4<-ts(dat3, start=c(1960,1), freq=12)
dat4
plot(dat4)
plot(dat)
}
