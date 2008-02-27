extend <- function(...)
	UseMethod("extend")

extend.list <- function(this, subclass, ..., listV) {
	class(this) <- c(subclass, class(this))
	if(missing(listV))
		return(structure(c(this, list(...)), class=class(this)))
	else
		return(structure(c(this, listV), class=class(this)))
}

MyEnv <- function(...) {
	this <- new.env()
	vars <- list(...)
	nms <- names(vars)
	if(length(nms)) for(i in 1:length(nms))
		assign(nms[i], vars[[i]], env=this)
	structure(this, class=c("MyEnv",class(this)))
}

extend.MyEnv <- function(this, subclass, ...) {
	class(this) <- c(subclass, class(this))
	newvars <- list(...)
	if(length(names(newvars))) for(nm in names(newvars))
		assign(nm, newvars[[nm]], env=this)
	return(this)
}

availableModels <- function()
	fitters

latex <- function(obj, ...)
	UseMethod("latex")

formatSignedNum <- function(x, ...) {
  signChar <- c("-","+")[(x>=0)+1]
  nm <- names(x)
  res <- paste(signChar,format(abs(x), ...) )
  names(res) <- nm
  return(res)
}

build <- function(...)
	UseMethod("build")

add <- function(...)
	UseMethod("add")

sigmoid <- function(x) 1/(1 + exp(-x))

dsigmoid <- function(x) sigmoid(x) * (1 - sigmoid(x))

d2sigmoid <- function(x) dsigmoid(x) * (1 - 2 * sigmoid(x))

repmat <- function(a, b, c) kronecker(matrix(1,b,c), a)

###Function to create the threshold in TVAR
TVAR_thresh<-function(mTh,thDelay,thVar=NULL,y, p){
T <- nrow(y) 
k <- ncol(y) 
if (!is.null(thVar)) {		
        if (length(thVar) > T) {
		z <- thVar[seq_len(T)]
		warning("The external threshold variable is not of same length as the original variable")
        }
        else
		z <- thVar
	z<-embed(z,p+1)[,seq_len(max(thDelay))+1]		#if thDelay=2, ncol(z)=2
} ###Combination (or single value indicating position) of contemporaneous variables
else {
	if (!length(mTh)%in%c(1,k))
		stop("length of 'mTh' should be equal to the number of variables, or just one")
	if(!all(mTh%in%seq_len(k)))
		stop("Unable to select the variable ",mTh[which(mTh%in%seq_len(k)==FALSE)], " for the threshold. Please see again mTh ")
	if(length(mTh)==1) {
		combin <- matrix(0,ncol=1, nrow=k)
		combin[mTh,]<-1
	}
	else 
		combin<-matrix(mTh,ncol=1, nrow=k)
	zcombin <- y %*% combin
	z <- embed(zcombin,p+1)[,seq_len(max(thDelay))+1]		#if thDelay=2, ncol(z)=2
}
zcombin <- y %*% combin
z <- embed(zcombin,p+1)[,seq_len(max(thDelay))+1]	
trans<-as.matrix(z)
list(trans=trans, combin=combin)
}

#Function to obtain the number of digits of a value
getndp <- function(x, tol=2*.Machine$double.eps)
{
internal<-function(x){
  	ndp <- 0
  	while(!isTRUE(all.equal(x, round(x, ndp), tol=tol))) ndp <- ndp+1 
  	if(ndp > -log10(tol)) warning("Tolerance reached, ndp possibly underestimated.")
  	ndp 
}
if(!is.null(dim(x)))
	x<-c(x)
if(length(x)==1)
	return(internal(x))
else if(length(x) %in% c(2:20))
	return(max(apply(matrix(x),1,internal)))
else
	return(max(apply(matrix(sample(x,size=20)),1,internal)))
}
###Paramater matrix of tar with 1 threshold,  given the thresh and the delay
TAR1t_B<-function(Delay,gamma,yy, xxl,xxh,z,m) {#
        isL <- ifelse(z[, Delay + 1]<= gamma,1,0)	### isL: dummy 
	ndown<-mean(isL)
        xxthresh <- cbind(xxl * isL,xxh * (1 - isL))	### Lower matrix
	B<-round(matrix(solve(crossprod(xxthresh))%*%crossprod(xxthresh,yy), nrow=1),5)
	Bcolnames <- c("Trend", c(paste("t -", seq_len(m))))
	colnames(B)<-rep(Bcolnames,2)
	Bdown <- B[,seq_len(ncol(B)/2)]
	Bup <- B[,-seq_len(ncol(B)/2)]
	nobs <- c(ndown=ndown, nup=1-ndown)	
	list(Bdown=Bdown, Bup=Bup, nobs=nobs)
        }
###Paramater matrix of tar with 2 thresholds,  given the 2 thresh and the delay
TAR2t_B <- function(gam1,gam2,Delay, yy, xx,z,m){
	##Threshold dummies
	dummydown <- ifelse(z[, Delay + 1]<=gam1, 1, 0)
	regimedown <- dummydown*xx
	ndown <- mean(dummydown)
	dummyup <- ifelse(z[, Delay + 1]>gam2, 1, 0)
	regimeup <- dummyup*xx
	nup <- mean(dummyup)
	##SSR from TAR(3)
	XX <- cbind(regimedown, (1-dummydown-dummyup)*xx, regimeup)		# dim k(p+1) x t	
	B <- round(matrix(solve(crossprod(XX))%*%crossprod(XX,yy),nrow=1),5)	#SSR
	Bcolnames <- c("Trend", c(paste("t -", seq_len(m))))
 	colnames(B)<-rep(Bcolnames,3)
	npar<-ncol(B)/3
	Bdown <- B[,c(1:npar)]
	Bmiddle <- B[,c(1:npar)+npar]
	Bup <- B[,c(1:npar)+2*npar]
	nobs <- c(ndown=ndown, nmiddle=1-ndown-nup,nup=nup)	
	list(Bdown=Bdown, Bmiddle=Bmiddle, Bup=Bup, nobs=nobs)
}