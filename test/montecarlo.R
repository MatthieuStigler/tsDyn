library(tsDyn);

noIter <- 20

#ncstarModels <- array(list(), noIter)
ncstarModels <- array(0, noIter)

sum <- array(0, 5)
phi1_mean <- array(0, 9)
th_mean <- array(0, 2)
gamma_mean <- array(0, 2)
noRegimes_mean <- array(0, noIter)
threeRegimes <- array(0, noIter)

for(i in 1:noIter) {
  cat("\n-------------------------------------------------------------\n",
      "DATASET NUMBER ", i,
      "\n-------------------------------------------------------------\n")
  
  ncstarModels[i] <-
    try(list(ncstar(ncstar.data[,i], m=2,noRegimes=10, control="method=BFGS,hessian = TRUE, maxit=1000")));

#  if(! is.na(ncstarModels[[i]]) ) {
#    if(ncstarModels[[i]]$model.specific$noRegimes == 3) {
#      threeRegimes[i] <- 1;
#3    }
    
    noRegimes_mean[i] <- ncstarModels[[i]]$noRegimes
    cat("\n--- Mean of the number of regimes so far: ",
        mean(noRegimes_mean[1:i]))
    
#    if(ncstarModels[[i]]$model.specific$noRegimes == 3) { 
      
#      if(i == 1) {
#        phi1_mean <- ncstarModels[[i]]$model.specific$phi1;
#        th_mean <- ncstarModels[[i]]$model.specific$phi2omega[,2];
#        gamma_mean<- ncstarModels[[i]]$model.specific$phi2omega[,1];
#      }
#      else {
#        phi1_mean <- apply(rbind(phi1_mean,
#               as.vector(ncstarModels[[i]]$model.specific$phi1)), 2, "mean")
#        th_mean <- apply(rbind(th_mean,
#               as.vector(ncstarModels[[i]]$model.specific$phi2omega[,2])), 2, "mean")
#        gamma_mean <- apply(rbind(gamma_mean,
#               as.vector(ncstarModels[[i]]$model.specific$phi2omega[,1])), 2, "mean")
#      }
      
#      cat("\n--- Mean of the linear parameters:\n    ", phi1_mean)
#      cat("\n--- Mean of the gamma parameters:\n     ", gamma_mean)
#      cat("\n--- Mean of the th parameters:\n     ", th_mean)

#    }
#  }
}
