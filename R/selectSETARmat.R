selectSETARmat<- function (x, m, d = 1, steps = d, thSteps = 7, mL = 1:m, mH = 1:m, 
    thDelay = seq_len(m)-1, criterion = c("pooled-AIC", "AIC", "SSR_OLS"), trim=0.15, ngrid="ALL", around, plot=TRUE,demean = c( "const", "trend","none", "both"), common=FALSE, model=c("TAR", "MTAR")) {
demean<-match.arg(demean)
if (max(thDelay) >= m)
	stop(paste("thDelay too high: should be < m (=", m, ")"))
str <- nlar.struct(x, m, d, steps)
xx <- getXX(str)
yy <- getYY(str)
model<-match.arg(model)

### Grid for threshold values
if(model=="TAR"){
	trans<-embed(x,m+1)}
else{
	if(max(thDelay)<m)
		trans<-embed(diff(x),m)
	else if(max(thDelay)==m){
		trans<-embed(diff(x),m+1)
		trans<-rbind(0,as.matrix(trans))}
}

allgammas <- sort(unique(trans[,1]))
ng <- length(allgammas)
if(ngrid=="ALL"){
	gammas <- allgammas[round(trim*ng):round((1-trim)*ng)]
} else {
	if(ngrid>ng)
		warning("There a not so many threshold values. Only ", ng, "values will be computed")
 	gammas <- allgammas[seq(from=trim, to=1-trim, length.out=ngrid)*ng]
}

if(!missing(around)){
	if(missing(ngrid)) {
		ngrid <- 20;
		cat("Searching for 20 values around\n")
	}
	wh.around <- which( round(allgammas,5) == round(around,5) )
	if(length(wh.around)==0)
		stop("Sorry, the value you gave for the around argument did not match")
	ar <- c((wh.around-round(ngrid/2)) : (wh.around+round(ngrid/2)))	#Values around the point
	gammas <- allgammas[ar[ar>=round(trim*ng)&ar<=round((1-trim)*ng)]]	#grid
}

th <- gammas

###Function pooled AIC
pooledAIC <- function(parms) {	
	thDelayVal <- parms[1] + 1
	mLVal <- parms[2]
	mHVal <- parms[3]
	m <- max(thDelayVal, mLVal, mHVal)
	lags <- c( (seq_len(m)-1) * (-d), steps)
	xxyy <- embedd(x, lags = lags)
	z <- xxyy[, thDelayVal]
	isLow <- (z <= parms[4])
	if ((sum(isLow) < mLVal) | (sum(!isLow) < mHVal)) 
	    return(NA)
	xx <- xxyy[isLow, seq_len(mLVal)]
	y <- xxyy[isLow, m + 1]
	AIC1 <- AIC(lm(y ~ xx))
	xx <- xxyy[!isLow, seq_len(mHVal)]
	y <- xxyy[!isLow, m + 1]
	AIC2 <- AIC(lm(y ~ xx))
	return(AIC1 + AIC2)
}

###Function for usual AIC
parsToModel <- function(parms) {
	thDelayVal <- parms[1]
	mLVal <- parms[2]
	mHVal <- parms[3]
	thVal <- parms[4]
	m <- max(thDelayVal + 1, mLVal, mHVal)
	return(setar(x, m = m, d = d, steps = steps, mL = mLVal, 
		mH = mHVal, th = thVal, thDelay=thDelayVal, demean=demean, common=common))
}

### SSR function
if(demean=="const"){
	const<-rep(1,nrow(xx))
	nconst<-"const"}
else if(demean=="trend"){
	const<-seq_len(nrow(xx))
	nconst<-"trend"}
else if(demean=="both"){
	const<-cbind(rep(1,nrow(xx)),seq_len(nrow(xx)))
	nconst<-c("const", "trend")}
else{
	const<-NULL
	nconst<-NULL}

SSRestim <- function(parameters) {
	Delay <- parameters[1]
	gammai <- parameters[2] 
        isL <- ifelse(trans[, Delay + 1]< gammai,1,0)	### isL: dummy variable
	if(common==FALSE){
		xxL <- cbind(const,xx[,seq_len(mL)])*isL
		xxH <- cbind(const,xx[,seq_len(mH)])*(1-isL)
		xxLH<-cbind(xxL,xxH)}
	else
		xxLH<-cbind(const,xx[,seq_len(mL)]*isL,xx[,seq_len(mH)]*(1-isL))
	crossprod(yy - xxLH %*% chol2inv( chol( crossprod(xxLH) ) ) %*% crossprod(xxLH, yy) )
}

###Grid for computation
    x <- str$x
    IDS <- as.matrix(expand.grid(thDelay,  mL, mH, th))			###Matrix of all combinations
	colnames(IDS)<-c("thDelay", "mL", "mH", "th")
    IDS2 <- as.matrix(expand.grid(thDelay, th))
	colnames(IDS2)<-c("thDelay", "th")
###Computation
criterion <- match.arg(criterion)
IDS <- switch(criterion, AIC = IDS, "pooled-AIC" = IDS, "SSR_OLS" = IDS2)	###Selection of the grid
if (criterion == "pooled-AIC") {
    computedCriterion <- apply(IDS, 1, pooledAIC)
} else if(criterion=="AIC"){
    critFun <- switch(criterion, AIC = AIC)
    computedCriterion <- apply(IDS, 1, function(x) critFun(parsToModel(x)))
} else if(criterion=="SSR_OLS"){
    computedCriterion <- apply(IDS, 1, SSRestim )
}

###Results
res <- cbind(IDS, computedCriterion)

###Graphical outuput
if(plot==TRUE){
	allcol <- seq_len(max(thDelay+1)*max(mL)*max(mH))
	col <- switch(criterion, AIC=allcol, "pooled-AIC"=allcol,"SSR_OLS"=(thDelay+1) )
	big <- apply(expand.grid(thDelay,mL, mH),1,function(a) paste("Th:", a[1],"mL:", a[2], "mH:", a[3]))
	legend <- switch(criterion, "AIC"=big, "pooled-AIC"=big, "SSR_OLS"=paste("Threshold Delay", thDelay))

	plot(res[,"th"], res[,"computedCriterion"], col=col, xlab="Treshold Value",ylab=criterion, main="Results of the grid search")
	legend("topleft", pch=1, legend=legend, col=col, bg=0)
}
### Results
colnames(res) <- c(colnames(IDS), criterion)
idSel <- sort(computedCriterion, index = TRUE)$ix
idSel <- idSel[seq_len(min(10, length(idSel)))]
res <- data.frame(res[idSel, ], row.names = NULL)

return(res)
}

if(FALSE) { #usage example
library(tsDyn)
environment(selectSETARmat)<-environment(selectNNET)


#Transformation like in Hansen 1999
sun<-(sqrt(sunspot.year+1)-1)*2		

###Full grid search with OLS
selectSETARmat(sun, m=3, criterion="SSR_OLS", d=1, thDelay=0:2,model="MTAR")

###restricted search with AIC or AIC pooled around the max selected by OLS
selectSETARmat(sun, m=2, criterion="AIC", d=1, thDelay=0:1, around=7.444575)
}

