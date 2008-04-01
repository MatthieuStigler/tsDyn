print.nlVar<-function(x){
	if(x$model.specific$nthresh==0) 
		cat("Linear VAR model\n")
	else
		cat("\n\nNon Linear Model\n")
}

logLik.nlVar<-function(x){
	res<-x$residuals
	k<-x$k
	t<-x$t
	Sigmabest<-matrix(1/t*crossprod(res),ncol=k)
	log(det(Sigmabest))
}

AIC.nlVar<-function(x, k=2){
	t<-x$t
	t*logLik.nlVar(x)+k*(x$nparB+x$model.specific$nthresh)
}

deviance.nlVar<-function(x){
	as.numeric(crossprod(c(x$residuals)))
}