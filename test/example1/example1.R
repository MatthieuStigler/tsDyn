datasize <- 1000; # Número de instancias del modelo que serán generadas
T <- 500; # Tamaño de cada instancia
svIter <- 100; # Número de iteraciones de la búsqueda de p. iniciales

alg <- "BFGS"
cluster <- NULL
#cluster <- c("localhost","localhost")
#library(snow)

plot <- TRUE

# PARÁMETROS DEL MODELO ORIGINAL
noRegimes <- 1
m <- 2 
sigma <- 1; # Varianza del modelo

phi1 <- c(0.8, -0.5, 0.3)

########################################################
########################################################

library(tsDyn)
ncstarModels <- array(list(), 2000)
nR <- array(NA, datasize)
nR_acc <- NA

first <- TRUE
for(i in 1:datasize) {
  cat("\n\n-------------------------------------------------------------\n",
      "DATASET NUMBER ", i, 
      "\n-------------------------------------------------------------\n")

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # BUILD UP A TIME SERIES
  
  cat("Building time series with ", noRegimes, "regimes...\n")
  
  e <- sigma * rnorm(T+501)

  y <- rnorm(2)
  f1 <- NA
  f2 <- NA

  for(t in 3:(T+500)) {
    x <- c(1, y[(t-1):(t-2)]);

    y[t] <- phi1 %*% x + e[t];
  }

  y <- y[501:(T+500)]

  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # ESTIMATE THE MODEL
  
  ncstarModels[i] <-
    try(list(ncstar(y, m, alg = alg,
                     cluster = cluster, svIter = svIter, trace = TRUE)));

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # COMPUTE MEDIANS

  nR[i] <- ncstarModels[[i]]$noRegimes
  nR_acc[i] <- mean(nR[1:i])
  cat("\n*** Mean of the number of regimes so far: ", mean(nR[1:i]))

  if(ncstarModels[[i]]$model.specific$noRegimes == noRegimes) { 

    if(first == TRUE) {
      phi1_acc <- as.vector(ncstarModels[[i]]$model.specific$phi1)
      phi1_median_acc <- phi1_acc
      first <- FALSE
    }
    else {
      cat("\n*** Median of the linear parameters:\n    ")
      phi1_acc <- rbind(phi1_acc,
                   as.vector(ncstarModels[[i]]$model.specific$phi1))
      phi1_median <- apply(phi1_acc, 2, median)
      phi1_median_acc <- rbind(phi1_median_acc, phi1_median)
      print(phi1_median)
      
      if(plot) {
        par(mfrow=c(2, 2), pty = "m")

        plot(nR[1:i], type="l", ylab="No. Regimes");
        lines(nR_acc, lty=2)
        abline(h=noRegimes, col="red")
        
        plot(phi1_acc[,1], type="l", ylab="phi_1");
        lines(phi1_median_acc[, 1], lty=2) 
        abline(h=phi1[1], col="red")
        
        plot(phi1_acc[,2], type="l", ylab="phi_2");
        lines(phi1_median_acc[, 2], lty=2) 
        abline(h=phi1[2], col="red")
        
        plot(phi1_acc[,3], type="l", ylab="phi_3");
        lines(phi1_median_acc[, 3], lty=2) 
        abline(h=phi1[3], col="red")
        
      }
      
    }

  }
}
