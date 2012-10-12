#############################################
############ Rank test
#############################################

rank.test <- function(vecm, type=c("trace", "eigen"), r_null=0:(vecm$k-1), cval=0.05){

  type <- match.arg(type)
  lambda <- vecm$model.specific$lambda
  t <- vecm$t
  k <- vecm$k

## trace test:
  trace <- -t*rev(cumsum(rev(log(1-lambda))))

## eigen 
  eigen <- -t*log(1-lambda)

## p-values:
  trace_pval <- gamma_doornik_all(trace, nmp=length(trace):1, test="H_lc", type="trace")
  adjT <- vecm$t -floor(npar(vecm)/vecm$k)
  trace_pval_T <- gamma_doornik_all(trace, nmp=length(trace):1, test="H_lc", type="trace", smallSamp=TRUE, T=adjT)
  eigen_pval <- gamma_doornik_all(eigen, nmp=length(eigen):1, test="H_lc", type="eigen")

## select r
  pvals <- switch(type, trace=trace_pval, eigen=eigen_pval)
  w.pvals <- which(pvals>cval)[1]
  rank <- w.pvals-1
#   if(length(r_null)==1)

## assemble
  res <- list()
  res_df <- data.frame(r=0:(k-1), trace=trace, trace_pval=trace_pval, trace_pval_T=trace_pval_T, eigen=eigen, eigen_pval=eigen_pval)
  res$res_df <- res_df
  res$r <- rank
  res$cval <- cval
  class(res) <- "rank.test"
  return(res)
}

print.rank.test <- function(x, ...) {

  cat("Rank selected:", x$r, "(first above critical value of", 100*x$cval, "%)\n")

}

summary.rank.test <- function(object, ...) {

  x<- object$res_df
  trace_pvalf <- format.pval(x$trace_pval, eps=1e-03,digits = max(1, getOption("digits") - 3))
  trace_pval_Tf <- format.pval(x$trace_pval_T, eps=1e-03,digits = max(1, getOption("digits") - 3))
  eigen_pvalf <- format.pval(x$eigen_pval, eps=1e-03,digits = max(1, getOption("digits") - 3))
  x_show <- data.frame(r=x$r, trace=x$trace, trace_pval=trace_pvalf,trace_pval_T =trace_pval_Tf,  eigen=x$eigen, eigen_pval=eigen_pvalf)

  print(x_show)

}

#############################################
############ P val approximation
#############################################

gamma_doornik <- function(x, nmp,q,  test=c("H_c", "H_lc", "H_l"), type=c("trace", "eigen"), T, smallSamp=FALSE){

  test <- match.arg(test)
  type <- match.arg(type)
  miss <- c(missing(x), missing(q))
  if(all(miss)|all(!miss)) stop("Please provide only one of args 'x' or 'q'\n")

### TRACE CASE

if(type=="trace"){
## Mean
  me_np <- switch(test, H_c=2.01, H_lc=1.05, H_l=4.05)
  one   <- switch(test, H_c=0, H_lc=-1.55, H_l=0.5)
  add   <- if(nmp ==1) switch(test, H_c=0.06, H_lc=-0.5, H_l=-0.23) else if(nmp==2) switch(test, H_c=0.05, H_lc=-0.23, H_l=-0.07) else 0
  mean <- 2*(nmp)^2 +me_np*(nmp) +one+add

## Var
  var_np <- switch(test, H_c=3.6, H_lc=1.8, H_l=5.7)
  var_one   <- switch(test, H_c=0.75, H_lc=0, H_l=3.2)
  var_add   <- if(nmp ==1) switch(test, H_c=-0.4, H_lc=-2.8, H_l=-1.3) else if(nmp==2) switch(test, H_c=-0.3, H_lc=-1.1, H_l=-0.5) else 0
  var <- 3*(nmp)^2 +var_np*(nmp) +var_one+ var_add

# print(c(mean, var))
### Small sample case:
  if(smallSamp){
# print(c(mean, var))
    paras_mean <- doornik_tab9_mean[,test]
    paras_var  <- doornik_tab9_var[,test]
    vals <- c(sqrt(nmp)/T,  nmp/T,  nmp^2/T^2, ifelse(nmp==1, 1/T,0), ifelse(nmp==1,1,0), ifelse(nmp==2,1,0), ifelse(nmp==3,1,0))
    mean_corr <- vals %*% paras_mean
    var_corr <- vals %*% paras_var
    mean <- exp(log(mean) +mean_corr)
    var <-  exp(log(var) +var_corr)
# print(c(mean, var))
  }

} else {

### EIGEN CASE

## Mean
  me_np <- switch(test, H_c=5.9498, H_lc=5.8271 , H_l=5.8658)
  one   <- switch(test, H_c=0.43402, H_lc=-1.6487 , H_l=2.5595 )
  add1   <- if(nmp== 1) switch(test, H_c=0.048360, H_lc= -1.6118 , H_l=-0.34443) else 0
  add2   <- if(nmp== 2) switch(test, H_c=0.018198, H_lc=-0.25949, H_l=-0.077991) else 0
  square <- switch(test, H_c=-2.3669, H_lc=-1.5666 , H_l=-1.7552 )
  mean <- me_np*(nmp) +one+add1+add2 +sqrt(nmp)*square 

## Var
  var_np <- switch(test, H_c=2.2231, H_lc=2.0785, H_l=1.9955 )
  var_one   <- switch(test, H_c=-7.9064, H_lc=-9.7846 , H_l=-5.5428 )
  var_add1   <- if(nmp== 1) switch(test, H_c=0.58592, H_lc=-3.3680, H_l=1.2425) else 0
  var_add2   <- if(nmp== 2) switch(test, H_c=-0.034324, H_lc=-0.24528, H_l=0.41949 ) else 0
  var_square <- switch(test, H_c=12.058, H_lc=13.074, H_l=12.841 )
  var <- var_np*(nmp) +var_one+ var_add1+var_add2 +sqrt(nmp)*var_square 
}

## dfs
  df1 <- mean^2/var
  df2 <- mean/var

## Compute p val or quantiles
  if(!missing(x)) res <- 1-pgamma(x, df1, df2)
  if(!missing(q)) res <- qgamma(q, df1, df2)
  names_res <- if(!missing(x)) "pval" else paste(100*q, "%", sep="")
  names(res) <- names_res

## return res
  return(res)
}


gamma_doornik_all <- function(x, nmp,q,  test=c("H_c", "H_lc", "H_l"), type=c("trace", "eigen"), smallSamp=FALSE, T){

## small checks
  test <- match.arg(test)
  type <- match.arg(type)
  miss <- c(missing(x), missing(q))
  if(all(miss)|all(!miss)) stop("Provide only one of args 'x' or 'q'\n")

## many x
  if(!missing(x)){
    res <- vector("numeric", length(x))
    for(i in 1:length(x)) res[i] <- gamma_doornik(x=x[i], nmp=nmp[i], test=test, type=type, smallSamp= smallSamp, T=T)
    names(res) <-paste("pval n-p", nmp, sep="=")
  } else {
    res <- matrix(NA, ncol=length(q), nrow=length(nmp))
    for(i in 1:length(nmp)) res[i,] <- gamma_doornik(q=q, nmp=nmp[i], test=test, type=type, smallSamp= smallSamp, T=T)
    rownames(res) <- paste("n-p", nmp, sep="=")
    colnames(res) <- paste(100*q, "%", sep="")
  }
return(res)
}


### Critical values tables (taken from Doornik, in gretl plugin/johansen.c 
Hnames <- c("H_z", "H_c", "H_lc", "H_l", "H_ql")

doornik_tab9_mean <- matrix(c(-0.101, 0.499, 0.896, -0.562, 0.00229, 0.00662, 0, 
0, 0.465, 0.984, -0.273, -244, 0, 0, 0.134, 0.422, 1.02, 2.17, 
-0.00182, 0, -0.00321, 0.0252, 0.448, 1.09, -0.353, 0, 0, 0, 
-0.819, 0.615, 0.896, 2.43, 0.00149, 0, 0), ncol=5, dimnames=list(1:7, Hnames))

doornik_tab9_var <- matrix(c(-0.204, 0.98, 3.11, -2.14, 0.0499, -0.0103, -0.00902, 
0.224, 0.863, 3.38, -0.807, 0, 0, -0.0091, 0.422, 0.734, 3.76, 
4.32, -0.00606, 0, -0.00718, 0, 0.836, 3.99, -1.33, -0.00298, 
-0.00139, -0.00268, -1.29, 1.01, 3.92, 4.67, 0.00484, -0.00127, 
-0.0199), ncol=5, dimnames=list(1:7, Hnames))





##########################################################################################
############ TEST: reproduce table in sec 10.1 of Doornik (1998) (with corr Doornik 1999) 
##########################################################################################

if(FALSE){
## Individual test of p-vals:
gamma_doornik(49.14, nmp=4)
gamma_doornik(19.06, nmp=3)
gamma_doornik(8.89, nmp=2)
gamma_doornik(2.35, nmp=1)

## Individual test of quantiles:
gamma_doornik(q=c(0.9, 0.95, 0.99), nmp=4)

## Vector test of p-vals/quantiles:
gamma_doornik_all(x=c(49.14, 19.06, 8.89,2.35), nmp=4:1)
gamma_doornik_all(q=c(0.9, 0.95, 0.99), nmp=4:1)

}


##########################################################################################
############ TEST: rank.test 
##########################################################################################

if(FALSE){
library(tsDyn)
library(vars)
library(cajo)

environment(rank.test) <- environment(star)

data(Canada)
data(denmark)

ve_can <- VECM(Canada, lag=1, estim="ML")
ve_can_l2 <- VECM(Canada, lag=2, estim="ML")

r_l1<- rank.test(ve_can)
r_l1
summary(r_l1)


r_l2<- rank.test(ve_can_l2)
r_l2b<- rank.test(ve_can_l2, cval=0.15)
r_l2_ei<- rank.test(ve_can_l2, type="eigen", cval=0.25)
r_l2
r_l2b
r_l2_ei
summary(r_l2)



}
