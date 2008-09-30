datasize <- 1000; # Número de instancias del modelo que serán generadas
T <- 500; # Tamaño de cada instancia
svIter <- 1000; # Número de iteraciones de la búsqueda de p. iniciales

alg <- "BFGS"
#cluster <- NULL
#cluster <- c("brain","brain")
cluster <- c("brain", "brain", "neuron1", "neuron1", "neuron2", "neuron2","neuron3", "neuron3","neuron4", "neuron4","neuron5", "neuron5","neuron6", "neuron6","neuron7", "neuron7")
library(snow)

plot <- FALSE

# PARAMETROS DEL MODELO ORIGINAL
noRegimes <- 3 
m <- 2 
sigma <- 0.5; # Varianza del modelo

gamma <- c(3.13, 2.12) 
th <- rbind(c(-0.5016, 0.5016),
            c(0.5152, -0.5152))
phi <- rbind(c(0.5, 0.8, -0.2),
             c(1.5, -0.6, -0.3),
             c(-0.5, -1.2, 0.7))

########################################################
########################################################

library(tsDyn)
ncstarModels <- array(list(), 2000)
noRegimes_mean <- array(NA, datasize)
phi1_acc <- array(NA, noRegimes * (m + 1))
phi2_acc <- array(NA, (m+1) * (noRegimes - 1))
phi2_acc_median <- array(NA, (m+1) * (noRegimes - 1))

first <- TRUE
for(i in 1:datasize) {
  cat("\n-------------------------------------------------------------\n",
      "DATASET NUMBER ", i, 
      "\n-------------------------------------------------------------\n")

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # BUILD UP A TIME SERIES
  
  e <- sigma * rnorm(T+501)

  y <- rnorm(2)
  f1 <- NA
  f2 <- NA

  for(t in 3:(T+500)) {
    f1[t] <- prod(exp(- gamma[1] * (y[(t-1):(t-2)] - th[1,])^2))
    f2[t] <- prod(exp(- gamma[2] * (y[(t-1):(t-2)] - th[2,])^2))
 
    x <- c(1, y[(t-1):(t-2)]);

    y[t] <- phi[1,] %*% x +
        phi[2,] %*% x * f1[t] +
        phi[3,] %*% x * f2[t] + e[t];
  }

  y <- y[501:(T+500)]

  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # ESTIMATE THE MODEL
  
  ncstarModels[i] <-
    try(list(ncgstar(y, m, alg = alg,
                     cluster = cluster, svIter = svIter, trace = TRUE)));

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # COMPUTE MEDIANS

  noRegimes_mean[i] <- ncstarModels[[i]]$noRegimes
  cat("\n*** Mean of the number of regimes so far: ",
      mean(noRegimes_mean[1:i]))

  if(ncstarModels[[i]]$model.specific$noRegimes == noRegimes) { 

    if(first == TRUE) {
      phi1_acc <- as.vector(ncstarModels[[i]]$model.specific$phi1)
      phi2_acc <- as.vector(ncstarModels[[i]]$model.specific$phi2)
      phi2_acc_median <- phi2_acc
      first <- FALSE
    }
    else {
      cat("\n*** Median of the linear parameters:\n    ")
      phi1_acc <- rbind(phi1_acc,
                   as.vector(ncstarModels[[i]]$model.specific$phi1))
      phi1_median <- apply(phi1_acc, 2, median)
      dim(phi1_median) <- c(3,3)
      print(phi1_median)
      
      cat("\n*** Median of the nonlinear parameters:\n     ")
      phi2_acc <- rbind(phi2_acc,
                   as.vector(ncstarModels[[i]]$model.specific$phi2))
      phi2_median <- apply(phi2_acc, 2, median) 
      phi2_acc_median <- rbind(phi2_acc_median, phi2_median)

      cat("gamma = ", phi2_median[1:(noRegimes - 1)])

      th_median <- phi2_median[noRegimes:length(phi2_median)]
      dim(th_median) <- c(m, noRegimes - 1)
      cat("\nth = ")
      print(th_median)

      if(plot) {
        par(mfrow=c(3, 2), pty = "m")

        plot(phi2_acc_median[,1], type="l", ylab="Gamma_1");
        lines(phi2_acc[,1], lty=2)
        abline(h=gamma[1], col="red")
        
        plot(phi2_acc_median[,2], type="l", ylab="Gamma_2");
        lines(phi2_acc[,2], lty=2) 
        abline(h=gamma[2], col="red")
        
        plot(phi2_acc_median[,3], type="l", ylab="Th_11");
        lines(phi2_acc[,3], lty=2)
        abline(h=th[1,1], col="red")
        
        plot(phi2_acc_median[,4], type="l", ylab="Th_12");
        lines(phi2_acc[,4], lty=2)
        abline(h=th[1,2], col="red")

        plot(phi2_acc_median[,5], type="l", ylab="Th_21");
        lines(phi2_acc[,5], lty=2)
        abline(h=th[2,1], col="red")
        
        plot(phi2_acc_median[,6], type="l", ylab="Th_22");
        lines(phi2_acc[,6], lty=2)
        abline(h=th[2,2], col="red")
      }
      
    }

  }
}
