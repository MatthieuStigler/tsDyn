datasize <- 1000; # Número de instancias del modelo que serán generadas
T <- 500; # Tamaño de cada instancia
svIter <- 1000; # Número de iteraciones de la búsqueda de p. iniciales

alg <- "BFGS"
#cluster <- NULL
cluster <- c("localhost","localhost")
library(snow)

plot <- TRUE

# PARÁMETROS DEL MODELO ORIGINAL
noRegimes <- 3 
m <- 2 
sigma <- 1; # Varianza del modelo

gamma <- c(8.49, 8.49) 
th <- c(-1.0607, 1.0607)
omega <- cbind(c(0.7071, -0.7071), c(0.7071,-0.7071))
phi <- rbind(c(0.5, 0.8, -0.2),
             c(1.5, -0.6, -0.3),
             c(-0.5, -1.2, 0.7))

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
    f1[t] <- 1 / (1 + exp(- gamma[1] *
                          (omega[1,1] * y[t-1] + omega[2,1] * y[t-2] - th[1])))
    f2[t] <- 1 / (1 + exp(- gamma[2] *
                          (omega[1,2] * y[t-1] + omega[2,2] * y[t-2] - th[2])))
 
    x <- c(1, y[(t-1):(t-2)]);

    y[t] <- phi[1,] %*% x +
        phi[2,] %*% x * f1[t] +
        phi[3,] %*% x * f2[t] + e[t];
  }

  y <- y[501:(T+500)]

  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # ESTIMATE THE MODEL
  
  ncstarModels[i] <-
    try(list(ncstar(y, m, noRegimes = noRegimes, alg = alg,
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
      phi2_acc <- as.vector(ncstarModels[[i]]$model.specific$phi2)
      phi2_median_acc <- phi2_acc
      first <- FALSE
    }
    else {
      cat("\n*** Median of the linear parameters:\n    ")
      phi1_acc <- rbind(phi1_acc,
                   as.vector(ncstarModels[[i]]$model.specific$phi1))
      phi1_median <- apply(phi1_acc, 2, median)
      dim(phi1_median) <- c(noRegimes, m+1)
      print(phi1_median)
      phi1_median_acc <- rbind(phi1_median_acc, phi1_median)
      
      cat("\n*** Median of the nonlinear parameters:\n")
      phi2_acc <- rbind(phi2_acc,
                   as.vector(ncstarModels[[i]]$model.specific$phi2))
      phi2_median <- apply(phi2_acc, 2, median) 
      phi2_median_acc <- rbind(phi2_median_acc, phi2_median)

      cat("gamma = ", phi2_median[1:(noRegimes - 1)])
      cat("; th = ", phi2_median[noRegimes:(2*(noRegimes - 1))])
      cat("\nomega=", phi2_median[(2 * (noRegimes - 1) + 1):length(phi2_median)])

      if(plot) {
        par(mfrow=c(2, 2), pty = "m")

        plot(phi2_acc[,1], type="l", ylab="Gamma_1");
        lines(phi2_median_acc[,1], lty=2)
        abline(h=gamma[1], col="red")
        
        plot(phi2_acc[,2], type="l", ylab="Gamma_2");
        lines(phi2_median_acc[,2], lty=2) 
        abline(h=gamma[2], col="red")
        
        plot(phi2_acc[,3], type="l", ylab="Th_1");
        lines(phi2_median_acc[,3], lty=2)
        abline(h=th[1], col="red")
        
        plot(phi2_acc[,4], type="l", ylab="Th_2")
        lines(phi2_median_acc[,4], lty=2)
        abline(h=th[2], col="red")

      }
      
    }

  }
}
