linear2<-function(data, lag, include = c( "const", "trend","none", "both"), model=c("VAR", "VECM"), I=c("level", "diff"),beta=NULL, coINT=FALSE){
y <- as.matrix(data)
Torigin <- nrow(y) 	#Size of original sample
T <- nrow(y) 		#Size of start sample
p <- lag
t <- T-p 		#Size of end sample
k <- ncol(y) 		#Number of variables
t<-T-p			#Size of end sample
ndig<-getndp(y)
if(is.null(colnames(data)))
	colnames(data)<-paste("Var", c(1:k), sep="")

include<-match.arg(include)
model<-match.arg(model)
I<-match.arg(I)
Y <- y[(p+1):T,] #
X <- embed(y, p+1)[, -seq_len(k)]	#Lags matrix

DeltaY<-diff(y)[(p+1):(T-1),]
Xminus1<-embed(y,p+2)[,(k+1):(k+k)]
DeltaX<-embed(diff(y),p+1)[,-(1:k)]

if(model=="VAR"){
	Z<-X
	Y<-Y}
if(model=="VECM"|I=="diff"){
	Z<-DeltaX
	Y<-DeltaY
	t<-t-1}

##Long-run relationship OLS estimation
if(model=="VECM"){
	if(is.null(beta) ){
		if(coINT)
			coint<-lm(y[,1]~ y[,-1])
		else
			coint<-lm(y[,1]~ y[,-1]-1)
	}
	else
		coint<-c(1, -beta)

# beta0<-(beta0%*%coint$coef[-1])[-c(1:p,T),]}

	betaLT<-coint$coef
	betaLT_std <- sqrt(diag(summary(coint)$sigma*summary(coint)$cov))


# ECT<-y%*%c(1,-betaLT)
# ECT<-round(ECT,ndig)
#ECTminus1<-ECT[-c(1:p,T)]

	ECTminus1<-round(Xminus1%*%c(1,-betaLT),ndig)
	Z<-cbind(ECTminus1,Z)
}

###Regressors matrix


if(include=="const")
	Z<-cbind(1, Z)
else if(include=="trend")
	Z<-cbind(seq_len(t), Z)
else if(include=="both")
	Z<-cbind(rep(1,t),seq_len(t), Z)
##################
###Linear model
#################
 B<-t(Y)%*%Z%*%solve(t(Z)%*%Z)		#B: OLS parameters, dim 2 x npar


npar<-ncol(B)*nrow(B)

rownames(B)<-paste("Equation",colnames(data))
LagNames<-c(paste(rep(colnames(data),p), -rep(seq_len(p), each=k)))


if(model=="VECM")
	ECT<-"ECT"
else
	ECT<-NULL


if(include=="const")
	Bnames<-c("Intercept",ECT, LagNames)
else if(include=="trend")
	Bnames<-c("Trend",ECT, LagNames)
else if(include=="both")
	Bnames<-c("Intercept","Trend",ECT, LagNames)
else 
	Bnames<-c(LagNames)

colnames(B)<-Bnames
fitted<-Z%*%t(B)
res<-Y-fitted

model.specific<-list()
model.specific$nthresh<-0
if(model=="VECM"){
	model.specific$betaLT<-betaLT
	model.specific$betaLT_std<-betaLT_std}


z<-list(residuals=res,  coefficients=B,  k=k, t=t,T=T, npar=npar, nparB=ncol(B), type="linear", fitted.values=fitted, model.x=Z, include=include,lag=lag, model=model, model.specific=model.specific)
class(z)<-c("VAR","nlVar")
return(z)
}


if(FALSE) { #usage example
###Hansen Seo data
library(tsDyn)
environment(linear2)<-environment(star)
data(zeroyld)
dat<-zeroyld

aVAR<-linear2(dat[1:100,], lag=2, include="both", model="VECM")
#lag2, 2 thresh, trim00.05: 561.46
class(aVAR)
aVAR
print(aVAR)
logLik(aVAR)
AIC(aVAR)
BIC(aVAR)
deviance(aVAR)
coef(aVAR)
summary(aVAR)
toLatex(aVAR)
toLatex(summary(aVAR))
}




print.VAR<-function(x,...){
	print(coef(x))
}

summary.VAR<-function(object, digits=4,...){
	x<-object
	r<-4
	t<-x$t
	k<-x$k

	Sigma<-matrix(1/t*crossprod(x$residuals),ncol=k)
	VarCovB<-solve(crossprod(x$model.x))%x%Sigma
	StDevB<-matrix(diag(VarCovB)^0.5, nrow=k)

	Tvalue<-x$coefficients/StDevB

	Pval<-pt(abs(Tvalue), df=(nrow(x$model.x)-ncol(x$model.x)), lower.tail=FALSE)+pt(-abs(Tvalue), df=(nrow(x$model.x)-ncol(x$model.x)), lower.tail=TRUE)
	#Pval<-round(Pval,4)
	symp <- symnum(Pval, corr=FALSE,cutpoints = c(0,  .001,.01,.05, .1, 1), symbols = c("***","** ","*  ",".  ","    "))
	stars<-matrix(symp, nrow=nrow(Pval))
	ab<-matrix(paste(myformat(x$coefficients,digits),"(", myformat(StDevB,digits),")",stars,sep=""), nrow=nrow(Pval))
	dimnames(ab)<-dimnames(x$coefficients)		

	x$bigcoefficients<-ab
	x$Sigma<-Sigma
	x$StDevB<-StDevB
	x$Pvalues<-Pval
	x$aic<-AIC.nlVar(x)
	x$bic<-BIC.nlVar(x)
	class(x)<-c("summary.VAR", "VAR")
	return(x)
}



print.summary.VAR<-function(x,...){
	cat("#############\n###Model", x$model,"\n#############")
	cat("\nFull sample size:",x$T, "\tEnd sample size:", x$t) 
	cat("\nNumber of variables:", x$k,"\tNumber of estimated slope parameters", x$npar)
	cat("\nAIC",x$aic , "\tBIC", x$bic)
	if(x$model=="VECM")
		cat("\nCointegrating vector:", x$model.specific$betaLT)
	cat("\n\n")
	print(noquote(x$bigcoefficients))

}


toLatex.VAR<-function(object,..., digits=4){
	x<-object
	if(inherits(x,"summary.VAR")){
		coef<-x$bigcoefficients
		a<-as.numeric(sub("\\([[:print:]]*", "",coef)) #extract coef values
		a<-myformat(a,digits,toLatex=TRUE)
		b<-as.numeric(sub("\\).*", "",sub(".*\\(", "",coef))) #extract st dev
		b<-myformat(b,digits,toLatex=TRUE) #put the scientific notation in \text latex fromat
		if(getOption("show.signif.stars"))					d<-paste("^{",sub(".*\\)", "",coef),"}", sep="")#extract stars and add ^{}
		else
			d<-NULL
		coef<-matrix(paste(a,"(",b,")",d, sep=""),ncol=ncol(coef), nrow=nrow(coef))
		}
	else{
		coef<-myformat(x$coefficients, digits)}
	varNames<-rownames(x$coefficients)
	res<-character()
	res[1]<-"%insert in the preamble and uncomment the line you want for usual /medium /small matrix"
	res[2]<-"%\\usepackage{amsmath} \\newenvironment{smatrix}{\\begin{pmatrix}}{\\end{pmatrix}} %USUAL"
	res[3]<-"%\\usepackage{amsmath} \\newenvironment{smatrix}{\\left(\\begin{smallmatrix}}{\\end{smallmatrix}\\right)} %SMALL"
	res[4]<-"%\\usepackage{nccmath} \\newenvironment{smatrix}{\\left(\\begin{mmatrix}}{\\end{mmatrix}\\right)} %MEDIUM"
	res[5]<-"\\begin{equation}"
	res[6]<- "\\begin{smatrix} %explained vector"
	res[7]<-TeXVec(paste("X_{t}^{",seq(1, x$k),"}", sep=""))
	res[8]<- "\\end{smatrix}="
 	res<-include(x, res, coef)
	ninc<-switch(x$include, "const"=1, "trend"=1,"none"=0, "both"=2)
	if(x$model=="VECM"){
		len<-length(res)
		res[len+1]<-"+\\begin{smatrix}  %ECT"
		res[len+2]<-TeXVec(coef[,ninc+1])
		res[len+3]<-"\\end{smatrix}ECT_{-1}"
		ninc<-ninc+1}
	res<-LagTeX(res, x, coef, ninc)
	res[length(res)+1]<-"\\end{equation}"
	res<-gsub("slash", "\\", res, fixed=TRUE)
	res<-res[res!="blank"]
	
	return(structure(res, class="Latex"))
}



