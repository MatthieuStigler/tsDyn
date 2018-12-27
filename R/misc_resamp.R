resamp <- function(x, boot.scheme=c("resample", "wild1", "wild2", "check"), seed = NULL){
  
  boot.scheme <- match.arg(boot.scheme)
  
  ## everything is matrix here
  X <-  as.matrix(x)
  t <-  nrow(X)

  if(!is.null(seed)) set.seed(seed) 
  X_boot <- switch(boot.scheme, 
                  "resample" =  X[sample(seq_len(t), replace=TRUE),,drop= FALSE], 
                  "wild1" = X+rnorm(t), 
                  "wild2" = X+sample(c(-1,1), size = t, replace=TRUE),
                  "check" = X)
  if(!is.matrix(x)) X_boot <- X_boot[,1]
  X_boot
}

if(FALSE) {
  head(resamp(x=lynx, seed = 123))
  head(resamp(cbind(lynx, lynx), seed = 123))
}