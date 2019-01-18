###

model.gen <- function(model = c("setar", "TVAR", "TVECM"),
                      B, n=200, lag=1, include=c("const", 'trend', "none", "both"), 
                      nthresh=0, thDelay=0, Thresh, 
                      trendStart=1, 
                      starting=NULL,  innov, 
                      exo, # setar
                      round_digits = 10,
                      returnStarting = FALSE,
                      add.regime =FALSE,
                      show.parMat = FALSE, 
                      thVar=NULL, mTh=1, # TVAR
                      beta, #TVECM
                      ...){
 
  model <- match.arg(model) 
  if(model == "setar") {
    res <- setar.gen(B = B, n=n, lag=lag, include=include, 
                     nthresh=nthresh, Thresh=Thresh, 
                     thDelay=thDelay, ## SPECIFIC to: setar, TVAR
                     starting=as.numeric(starting), innov = innov, 
                     exo = exo, ## SPECIFIC
                     trendStart=trendStart, ## SPECIFIC to: setar, TVAR
                     round_digits = round_digits,
                     returnStarting = returnStarting,
                     add.regime =add.regime,
                     show.parMat = show.parMat, ...)
    
  } else if(model == "TVAR") {
    res <- TVAR.gen(B = B, n=n, lag=lag, include=include, 
                    nthresh=nthresh, Thresh=Thresh, 
                    thDelay=thDelay, ## SPECIFIC to: setar, TVAR
                    starting=starting, innov = innov, 
                    thVar=thVar, mTh=mTh,  ## SPECIFIC
                    trendStart=trendStart, ## SPECIFIC to: setar, TVAR
                    round_digits=round_digits,
                    returnStarting = returnStarting,
                    add.regime =add.regime,
                    show.parMat=show.parMat, 
                    ...)
  } else if(model == "TVECM") {
    ## not used: thDelay
    res <- TVECM.gen(B = B, n=n, lag=lag, include=include, 
                     beta = beta, ## SPECIFIC
                     nthresh=nthresh, Thresh=Thresh, 
                     starting=starting, innov = innov, 
                     round_digits=round_digits,
                     returnStarting = returnStarting,
                     show.parMat=show.parMat, 
                     ...)
      
  }
  as.data.frame(res)
}
  


  
  
    