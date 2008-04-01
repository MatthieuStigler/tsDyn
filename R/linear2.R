linear2<-function(data, lag, include = c( "const", "trend","none", "both")){
y <- as.matrix(data)
Torigin <- nrow(y) 	#Size of original sample
T <- nrow(y) 		#Size of start sample
p <- lag
t <- T-p 		#Size of end sample
k <- ncol(y) 		#Number of variables
t<-T-p			#Size of end sample
if(is.null(colnames(data)))
	colnames(data)<-paste("Var", c(1:k), sep="")
include<-match.arg(include)

Y <- y[(p+1):T,] #
Z <- embed(y, p+1)[, -seq_len(k)]	#Lags matrix

if(include=="const")
	Z<-cbind(1, Z)
else if(include=="trend")
	Z<-cbind(seq_len(t), Z)
else if(include=="both")
	Z<-cbind(rep(1,t),seq_len(t), Z)

npar <- ncol(Z)			#Number of parameters

##################
###Linear model
#################
 B<-t(Y)%*%Z%*%solve(t(Z)%*%Z)		#B: OLS parameters, dim 2 x npar


allpar<-ncol(B)*nrow(B)

rownames(B)<-paste("Equation",colnames(data))
LagNames<-c(paste(rep(colnames(data),p), -rep(seq_len(p), each=k)))

if(include=="const")
	Bnames<-c("Intercept",LagNames)
else if(include=="trend")
	Bnames<-c("Trend",LagNames)
else if(include=="both")
	Bnames<-c("Intercept","Trend",LagNames)
else 
	Bnames<-c(LagNames)
colnames(B)<-Bnames
fitted<-Z%*%t(B)
res<-Y-fitted


z<-list(residuals=res,  coefficients=B,  k=k, t=t, nparB=ncol(B), nthresh=0, type="linear", fitted.values=fitted)
class(z)<-"nlVar"
return(z)
}


if(FALSE) { #usage example
###Hansen Seo data
data(zeroyld)
data<-zeroyld

a<-linear2(data[1:100,], lag=2, include="none")
#lag2, 2 thresh, trim00.05: 561.46
class(a)
print(a)
logLik(a)
AIC(a)
deviance(a)
}
