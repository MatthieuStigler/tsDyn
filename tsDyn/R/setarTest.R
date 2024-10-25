#' Test of linearity against threshold (SETAR)
#' 
#' Test of linearity against threshold of Hansen (1999) with bootstrap
#' distribution
#' 
#' Estimation of the first threshold parameter is made with CLS, a conditional
#' search with one iteration is made for the second threshold. The Ftest
#' comparing the residual sum of squares (SSR) of each model is computed.
#' 
#' \deqn{ F_{ij}=T( (S_{i}-S_{j})/S_{j} )}
#' 
#' where \eqn{S_{i}} is the SSR of the model with i regimes (and so i-1
#' thresholds).
#' 
#' Three test are available. The both first can be seen as linearity test,
#' whereas the third can be seen as a specification test: once the 1vs2 or/and
#' 1vs3 rejected the linearity and henceforth accepted the presence of a
#' threshold, is a model with one or two thresholds preferable?
#' 
#' Test \bold{1vs}2: Linear AR versus 1 threshold TAR
#' 
#' Test \bold{1vs}3: Linear AR versus 2 thresholds TAR
#' 
#' Test \bold{2vs3}: 1 threshold TAR versus 2 thresholds TAR
#' 
#' The two first tests are computed together and available with test="1vs". The
#' third test is available with test="2vs3".
#' 
#' The homoskedastic bootstrap distribution is based on resampling the
#' residuals from H0 model (ar for test 1vs, and setar(1) for test 2vs3),
#' estimating the threshold parameter and then computing the Ftest, so it
#' involves many computations and is pretty slow.
#' 
#' @aliases setarTest setartest
#' @param x time series
#' @param m,thDelay lag and 'time delay' for the threshold variable
#' @param test whether to test AR against SETAR, or SETAR(1 reg) against SETAR(2 reg)
#' @template param_include
#' @param nboot number of bootstrap replications
#' @param trim trimming parameter indicating the minimal percentage of
#' observations in each regime
#' @template param_boot.scheme
#' @param hpc Possibility to run the bootstrap on parallel core. See details in
#' @param seed Seed used in the bootstrap resampling
#' \code{\link{TVECM.HStest}}
#' @return A object of class "Hansen99Test" containing:
#' 
#' \item{SSRs}{The residual Sum of squares of model AR, 1 threshold TAR and 2
#' thresholds TAR}
#' 
#' \item{Ftests}{The Ftest statistic for the test}
#' 
#' \item{PvalBoot}{The bootstrap p-values for the test selected}
#' \item{CriticalValBoot}{The critical values for the test selected}
#' \item{Ftestboot}{All the F-test computed} \item{firstBests, secBests}{The
#' thresholds for the original series, obtained from search for 1 thresh
#' (\code{firstBests}) and conditional search for 2 thresh (secBests)}
#' \item{nboot,m}{The number of bootstrap replications (\code{nboot}), the
#' lags used (\code{m})}
#' @author Matthieu Stigler
#' @seealso \code{\link{TVAR.LRtest}} for the multivariate version.
#' \code{\link{SETAR}} for estimation of the model.
#' @references Hansen (1999) Testing for linearity, Journal of Economic
#' Surveys, Volume 13, Number 5, December 1999 , pp. 551-576(26) available at:
#' \url{http://www.ssc.wisc.edu/~bhansen/papers/cv.htm}
#' @keywords ts
#' @examples
#' 
#' #Data used by Hansen
#' sun <- (sqrt(sunspot.year + 1) - 1) * 2
#' 
#' #Test 1vs2 and 1vs3
#' #setarTest(sun, m=11, thDelay=0:1, nboot=5,trim=0.1, test="1vs")
#' @name setarTest
NULL

test_AR_SETAR <-  function(x, m, thDelay, include, trim, return_th = FALSE) {
  
  ### setarTest 1: SSR and check linear model
  linear <- linear(x=x, m=m, include =include)
  
  ###setarTest 2: search best thresholds: call selectSETAR
  search <- selectSETAR(x, m=m, thDelay=thDelay, trace=FALSE, 
                        include =include, nthresh=2,
                        trim=trim, criterion = "SSR",
                        plot=FALSE, max.iter=3) 
  
  SSR <- deviance(linear)
  SSR1thresh <- search$firstBests["SSR"]
  SSR2thresh <- search$bests["SSR"]
  SSRs <- matrix(c(SSR, SSR1thresh, SSR2thresh), ncol=3, dimnames=list("SSR", c("AR", "TAR(1)", "TAR(2)")))
  if(return_th) SSRs <- cbind(SSRs, 
                              th_1 = search$firstBests["th"], 
                              th_2_1 = search$bests["th1"],
                              th_2_2 = search$bests["th2"])
  SSRs
}

get_tests <- function(object, n) {
  SSR <-  object[,"AR"]
  SSR_1thresh <- object[,"TAR(1)"]
  SSR_2thresh <- object[,"TAR(2)"]
  
  Ftest12 <- as.numeric(n * (SSR - SSR_1thresh) / SSR_1thresh)
  # Bound <- (n / (2 * log(log(n)))) * Ftest12
  Ftest13 <- as.numeric(n * (SSR - SSR_2thresh) / SSR_2thresh)
  Ftest23 <- as.numeric(n * (SSR_1thresh - SSR_2thresh) / SSR_2thresh)
  # Ftests <- matrix(c(Ftest12, Ftest13, Ftest23),ncol=3, dimnames=list("Ftest", c("1vs2", "1vs3", "2vs3")))
  Ftests <- c("1vs2" = Ftest12, "1vs3" = Ftest13, "2vs3" = Ftest23)
  Ftests
}

#' @rdname setarTest
#' @export
setarTest <- function(x, m, thDelay = 0, trim=0.1, 
                      include = c("const", "trend","none", "both"),
                      nboot = 10, 
                      test=c("1vs", "2vs3"), hpc=c("none", "foreach"), 
                      boot.scheme = c("resample", "resample_block", "wild1", "wild2", "check"),
                      seed = NULL){
  
  test <- match.arg(test)
  hpc <- match.arg(hpc)
  boot.scheme <- match.arg(boot.scheme)
  include <- match.arg(include)
  
  N <-  length(x)
  n <- N - m
  
  ### run test
  SSR_orig <- test_AR_SETAR(x=x, m=m, thDelay=thDelay, include= include, trim=trim, 
                            return_th = TRUE)
    
  Ftest_orig <- get_tests(object=SSR_orig, n = n)
  
  ## estimate H0 model (for later bootstrap)
  null_model <- switch(test, 
                       "1vs" = linear(x=x, m=m, include= include, warn_root=FALSE),
                       "2vs3"=  setar_low(x=x, m=m, include= include, 
                                          nthresh=1, thDelay = thDelay, th =SSR_orig[,"th_1"], warn_min_obs = FALSE))
  
  ## old check
  # a <- setar(x=x, m=m, include= include, nthresh=1, thDelay = thDelay, th =SSR_orig[,"th_1"])
  # aa <- setar(x=x, m=m, include= include, nthresh=2, thDelay = thDelay, th =SSR_orig[, c("th_2_1", "th_2_2")])
  # c(deviance(a), deviance(aa))
  
  
  ################
  ### boot
  ################
  if(nboot>0){   
    
    ## define function
    boot_test <- function(null_model, boot.scheme) {
      null_model_boot <- setar.boot(null_model, boot.scheme=boot.scheme)
      SSR_boot <- test_AR_SETAR(x=null_model_boot, m=m, thDelay=thDelay, include= include, trim=trim)
      get_tests(object=SSR_boot, n = n)
    }
    
    # boot_test(null_model, boot.scheme=boot.scheme)
    
    ### run bootstrap
    if(hpc=="none"){
      if(!is.null(seed)) set.seed(seed)
      Ftestboot <- replicate(nboot, boot_test(null_model, boot.scheme=boot.scheme), simplify = FALSE)
      Ftestboot <- do.call("rbind", Ftestboot)
    } else {
      if(!is.null(seed)) warning("Argument 'seed' will not work with arg. 'hpc'")
      Ftestboot <- foreach(i=1:nboot, .export="boot_test", .combine="rbind") %dopar% boot_test(null_model, boot.scheme)
    }
    
    ### analyse bootstrap stats
    probs <- c(0.9, 0.95, 0.975, 0.99)      
    if(test=="1vs"){
      Ftest_boot_12 <- Ftestboot[, 1]
      Ftest_boot_13 <- Ftestboot[, 2]
      
      Pval_boot_12 <- mean(Ftest_boot_12 > Ftest_orig["1vs2"])
      Pval_boot_13 <- mean(Ftest_boot_13 > Ftest_orig["1vs3"])
      
      Cval_boot_12 <- quantile(Ftest_boot_12, probs = probs)
      Cval_boot_13 <- quantile(Ftest_boot_13, probs = probs)
      
      CriticalValBoot <- matrix(c(Cval_boot_12,Cval_boot_13), nrow=2, byrow=TRUE, dimnames=list(c("1vs2", "1vs3"), probs))
      PvalBoot <- c(Pval_boot_12, Pval_boot_13)
    } else {
      Ftest_boot_23 <- Ftestboot[,3]
      Pval_boot_23 <- mean(Ftest_boot_23 > Ftest_orig["2vs3"])
      Cval_boot_23 <- quantile(Ftest_boot_23, probs = probs)
      CriticalValBoot <- matrix(Cval_boot_23, nrow=1, dimnames=list("2vs3", probs))
      PvalBoot <- Pval_boot_23
    }
  } else{ #nboot=0
    CriticalValBoot <- NULL
    PvalBoot <- NULL
    Ftestboot <- NULL
  }
  
  args <- list(x=x, m=m, thDelay = thDelay, nboot=nboot, trim=trim, boot.scheme=boot.scheme,
               test = test)
  
  ###res
  res <- list(Ftests = Ftest_orig,
              SSRs = SSR_orig[1:3],
              firstBests = SSR_orig[,"th_1"],
              secBests = SSR_orig[,c("th_2_1", "th_2_2")],
              CriticalValBoot = CriticalValBoot,
              PvalBoot = PvalBoot,
              Ftestboot = Ftestboot,
              nboot = nboot, args = args, 
              # Bound = Bound,
              updated = NULL)
  class(res) <- "Hansen99Test"
  return(res)
}
  
#' @export
print.Hansen99Test<-function(x,...){
  if(x$args$test=="1vs"){
    cat("Test of linearity against setar(2) and setar(3)\n\n")
    print(matrix(c(x$Ftests[-3], x$PvalBoot), ncol=2, dimnames=list(c("1vs2", "1vs3"), c("Test", "Pval"))))
  }
   else{
    cat("Test of setar(2) against  setar(3)\n\n")
    print(matrix(c(x$Ftests[3], x$PvalBoot), ncol=2, dimnames=list(c("2vs3"), c("Test", "Pval"))))
   }
}

#' @export
summary.Hansen99Test<-function(object, ...){
  print.Hansen99Test(object)
  cat("\nCritical values:\n")
  print(object$CriticalValBoot)
  cat("\nSSR of original series:\n")
  print(matrix(object$SSRs, ncol=1, dimnames=list(c("AR", "SETAR(2)", "SETAR(3)"), "SSR")))
  cat("\nThreshold of original series:\n")
  print(matrix(c(object$firstBests, NA, object$secBests),byrow=TRUE, ncol=2, dimnames=list(c("SETAR(2)", "SETAR(3)"), c("th1", "th2"))))
  cat("\nNumber of bootstrap replications: ", object$nboot, "\n")
  # if(object$args$test=="1vs"){
  #   cat("Asymptotic bound: ", object$Bound, "\n")
  # }
    
}

#' @export
plot.Hansen99Test<-function(x,show.extended=TRUE, ...){
  m <- x$args$m
  test <- x$args$test
  leg <- c("Asymptotic Chi 2", "Bootstrap", "Test value")
  col <- c(3, 1, 2)
  updated <- ifelse(!is.null(x$updated) & show.extended, TRUE, FALSE)
  if(updated){
    leg <- c(leg[1:3], "Old Bootstrap")
    col <- c(col, 4)
   }
  
  
  if(test=="1vs"){
    layout(c(1,2))
    Ftest_boot_12 <- x$Ftestboot[1, ]
    Ftest_boot_13 <- x$Ftestboot[2, ]
    Ftest12 <- x$Ftests[1]
    Ftest13 <- x$Ftests[2]

    #Plot 1vs2
    plot(density(Ftest_boot_12, from=0), xlab="Ftest12", xlim=c(0,max(Ftest12+1,max(Ftest_boot_12))),ylim=c(0,max(density(Ftest_boot_12)$y,dchisq(0:Ftest12, df=1+m))), main="")
    title("Test linear AR vs 1 threshold SETAR")
    abline(v=Ftest12, lty=2, col=2)
    curve(dchisq(x, df=1+m, ncp=0), from=0, n=Ftest12+5, add=TRUE, col=3)
    if(updated)
      lines(density(Ftest_boot_12[1:x$updated], from=0), col=4)
    legend("topright", legend=leg, col=col, lty=c(1,1,2), bg="white")
  
    #Plot 1vs3
    plot(density(Ftest_boot_13, from=0), xlab="Ftest13", xlim=c(0,max(Ftest13+1,max(Ftest_boot_13))),ylim=c(0,max(density(Ftest_boot_13)$y,dchisq(0:Ftest12, df=2*(1+m)))), main="")
    title("Test linear AR vs 2 thresholds SETAR")
    abline(v=Ftest13, lty=2, col=2)
    curve(dchisq(x, df=2*(1+m), ncp=0), from=0, n=Ftest13+5, add=TRUE, col=3)
    if(updated)
      lines(density(Ftest_boot_13[1:x$updated], from=0), col=4)
    legend("topright", legend=leg, col=col, lty=c(1,1,2), bg="white")
  }
  #plot 2vs3
  else {
    Ftest_boot_23 <- x$Ftestboot
    Ftest23 <- x$Ftests[3]
    
    plot(density(Ftest_boot_23, from=0), xlab="Ftest23", xlim=c(0,max(Ftest23+1,Ftest_boot_23)), ylim=c(0,max(density(Ftest_boot_23)$y)), main="")
    title("Test 1 threshold SETAR vs 2 thresholds SETAR")
    abline(v=Ftest23, lty=2, col=2)
    curve(dchisq(x, df=1+m, ncp=0), from=0, n=Ftest23+5, add=TRUE, col=3)
    if(updated)
      lines(density(Ftest_boot_23[1:x$updated], from=0), col=4)
    legend("topright", legend=leg, col=col, lty=c(1,1,2), bg="white")
  }
}


if(FALSE){  
  library(tsDyn)
  sun <- (sqrt(sunspot.year + 1) - 1) * 2
  
  test_1vs <- setarTest(x=sun, m=2, thDelay = 0, trim = 0.15, nboot = 10)
  test_2vs3 <- setarTest(x=sun, m=2, thDelay = 0, trim = 0.15, nboot = 10, test = "2vs3")
  
  print(test_1vs)
  summary(test_1vs)
  plot(test_1vs)
  
  print(test_2vs3)
  summary(test_2vs3)
  plot(test_2vs3)
  

  
  ###two probs
  #- deviance of selectSETAR and setar is nto the same
  # selectSETAR does not select the good one...
  environment(selectSETAR)<-environment(star)
  selectSETAR(lynx, m=11, thDelay=1, nthresh=2, th=list(exact=c(5.3,8)), criterion="SSR")
  selectSETAR(sun, m=11, thDelay=1, nthresh=2, th=list(exact=c(5.3,8)), criterion="SSR")
  
  ###US IP
  IP<-read.table(file="/media/sda5/Mes documents/Uni/Statistique/Time Series/Handbooks/datasets/Hansen/Hansen1999Linearity/Matlab/ipdat.txt")
  
  #transform as in Hansen 1999
  dat<-log(IP[,1])
  dat2<-diff(dat, 12)*100 # dat=(dat(13:length(dat(:,1)))-dat(1:length(dat(:,1))-12))*100
  dat3<-dat2[157:length(dat2)] #dat=dat(157:length(dat(:,1)));
  dat4<-ts(dat3, start=c(1960,1), freq=12)
  
  end(dat4)
  length(dat4)
  #plot
  plot(dat4)
  
  #save
  IIPUs<-dat4
  save(IIPUs, file="IIPUs.rda")
  
  ###Test extendBoot
  # test with 10 bootstrap replications:
  a<-setarTest(sun[1:100], m=1, nboot=10)
  plot(a)
  
  #use old results and compue 20 new replications
  b<-extendBoot(a, n=20)
  #see the different distributions:
  plot(b)
  
  plot(b, show.extended=FALSE)
}
