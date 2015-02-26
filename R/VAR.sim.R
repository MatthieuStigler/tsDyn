VAR.gen <- function(B, n=200, lag=1, include = c("const", "trend","none", "both"),  
                    starting=NULL, innov, exogen=NULL,
                    show.parMat=FALSE, returnStarting=FALSE,
                    seed){
  
  include<-match.arg(include)
  
  ## Create some variables/parameters
  p <- lag
  ninc <- switch(include, "none"=0, "const"=1, "trend"=1, "both"=2)
  k <- nrow(B)
  T <- n   	#Size of start sample
  t <- T-p  #Size of end sample
  y <- matrix(0,ncol=k, nrow=n+p)
  trend<-c(rep(NA, p), 1:T) 
  
  ## exogen
  if(!is.null(exogen)){
    if(!is.matrix(exogen)) exogen <- as.matrix(exogen)
    nExogen <- ncol(exogen)
  } else {
    exogen <- matrix(0, ncol=1, nrow=n+p)
    nExogen <- 0
  }
  

  ## Check inputs
  npar <- p*k+ninc+nExogen
  if(ncol(B)!=npar){
    stop("bad specification of B. Expected: ", 
         ninc, " const/trend, ",
         p, " * ", k, " lags, ", nExogen, " exogens")
  }
  
  if(!is.null(starting)&&!all(dim(as.matrix(starting))==c(p,k))){
    stop("Bad specification of starting values. Expected dim ", p, "*", k)
  }
  
  if(!all(dim(as.matrix(innov))==c(n,k))){
    stop("Bad specification of innovations. Expected dim ", n, "*", k)
  }
  
  
  ## Augment B to include always constant/trend (eventually 0)
  addInc <- switch(include, "none"=1:2, "trend"=1, "const"=2, "both"=NULL)
  Bfull <- myInsertCol(B, c=addInc ,0)
  if(nExogen==0) Bfull <- cbind(Bfull, 0)
  
  
  ## Starting values
  if(!is.null(starting)){
    if(all(dim(as.matrix(starting))==c(p,k)))
      y[seq_len(p),]<- as.matrix(starting)
    else
      stop("Bad specification of starting values. Should have nrow = lag and ncol = number of variables. But is: ", dim(as.matrix(starting)), sep="")
  }

  ## innovations
  resb <- rbind(matrix(0,nrow=p, ncol=k),innov)	
  
  ## MAIN loop:  
  for(i in (p+1):(n+p)){
    Y <- matrix(t(y[i-c(1:p),, drop=FALSE]), ncol=1)
    Yexo <-  rbind(Y, t(exogen[i,]))
    y[i,]<-rowSums(cbind(Bfull[,1],  # intercept
                         Bfull[,2]*trend[i], #trend
                         Bfull[,-c(1,2)]%*%Yexo, #lags
                         resb[i,])) #residuals
  }
  
  

  
  if(show.parMat) print(Bmat)
  if(!returnStarting) y <- y[-c(1:p),] 
  return(y)
}

if(FALSE){
  
  B <- matrix(c(0.3, 0.2, 0.1, 0.3, 0.2, 0.4),nrow=2 )
  n <- 200
  inno <- matrix(rnorm(200*2), ncol=2)
  environment(VAR.gen) <- environment(TVECM)
  VAR.gen(B=B, include="const", lag=1, innov=inno)
  
}

VAR.sim <- function(B, n=200, lag=1, include = c("const", "trend","none", "both"),  
                    starting=NULL, innov=rmnorm(n, mean=0, varcov=varcov), 
                    varcov=diag(1,k), 
                    show.parMat=FALSE, returnInitial=FALSE, seed){
  
  VAR.gen(B=B, n=n, lag=lag, include = include,  
          starting=NULL, innov=innov, 
          show.parMat=FALSE, returnInitial=returnInitial, 
          seed)
}
    
VAR.boot <- function(VARobject, boot.scheme=c("resample", "wild1", "wild2", "check"),
                     seed){
  
  boot.scheme <- match.arg(boot.scheme)
  
  B <- coef(VARobject)
  t <- VARobject$t
  k <- VARobject$k
  lags <- VARobject$lag
  include <- VARobject$include
  
  yorig <- VARobject$model[,1:k]
  starts <- yorig[1:lags,, drop=FALSE]
  resids <- residuals(VARobject)

  ## boot it
  innov <- switch(boot.scheme, 
                  "resample"=  resids[sample(seq_len(t), replace=TRUE),], 
                  "wild1"=resids+rnorm(t), 
                  "wild2"=resids+sample(c(-1,1), size = t, replace=TRUE),
                  "check"=  resids)
  
  res <- VAR.gen(B=B, n=t, lag=lags, include = include,  
                 starting=starts, innov=innov, 
                 show.parMat=FALSE, returnStarting=TRUE)
  colnames(res) <- colnames(yorig)
  res
}


if(FALSE){
  barry_mat <- as.matrix(as.data.frame(barry))
  
  va <- lineVar(barry, lag=1)
  a<-VAR.boot(va, boot.scheme="check")
  all.equal(a, as.matrix(as.data.frame(barry)), check.attributes = FALSE)
  
  var_l2_const <-VAR.boot(lineVar(barry, lag=2), boot.scheme="check")
  all.equal(var_l2_const, barry_mat)
  
  var_l3_const <-VAR.boot(lineVar(barry, lag=3), boot.scheme="check")
  all.equal(var_l3_const, barry_mat)
  
  var_l2_none <-VAR.boot(lineVar(barry, lag=2, include="none"), boot.scheme="check")
  all.equal(var_l2_none, barry_mat)

  var_l2_trend <-VAR.boot(lineVar(barry, lag=2, include="trend"), boot.scheme="check")
  all.equal(var_l2_trend, barry_mat)
  
}
