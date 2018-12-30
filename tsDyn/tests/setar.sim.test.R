library(tsDyn)

################
### tets boot
################

setar.boot.check <-  function(object, n_digits = 10) {
  mod_boot <- setar.boot(object, boot.scheme = "check",n_digits = n_digits)  
  orig_series <- as.numeric(object$str$x)
  all.equal(mod_boot$serie, orig_series)
}

linear.boot.check <-  function(object, n_digits = 10) {
  mod_boot <- linear.boot(object, boot.scheme = "check",n_digits = n_digits)  
  orig_series <- as.numeric(object$str$x)
  all.equal(mod_boot$serie, orig_series)
}



## nthresh ==0
ar_1_noInc <- linear(log(lynx), m = 1, include = "none")
ar_2_noInc <- linear(log(lynx), m = 2, include = "none")
ar_1_const <- linear(log(lynx), m = 1, include = "const")
ar_2_const <- linear(log(lynx), m = 2, include = "const")


setar.boot.check(ar_1_noInc)
setar.boot.check(ar_2_noInc)
setar.boot.check(ar_1_const)
setar.boot.check(ar_2_const)

linear.boot.check(ar_1_noInc)
linear.boot.check(ar_2_noInc)
linear.boot.check(ar_1_const)
linear.boot.check(ar_2_const)



## nthresh ==1
set_1th_l1 <-  setar(lynx, nthresh=1, m=1)
set_1th_l2 <-  setar(lynx, nthresh=1, m=2)
set_1th_l1_tr <-  setar(lynx, nthresh=1, m=1, include = "trend")


setar.boot.check(set_1th_l1)
setar.boot.check(set_1th_l1, n_digits = 2)
setar.boot.check(set_1th_l2)
setar.boot.check(set_1th_l2, n_digits = 5)
setar.boot.check(set_1th_l1_tr)
setar.boot.check(set_1th_l1_tr, n_digits = 1)


## why difference?
if(FALSE) {
  library(tidyverse)
  getTh(set_1th_l2)
  filt_diff <-  function(x, minus=2, tol =1) {
    x2 <- x %>% 
      mutate(diff = x$boot - x$orig)
    first <- which(abs(x2$diff)>tol)[1]
    filter(x2, between(n_row, first -minus, first +minus))
  }
  set_1th_l2_b <- setar.boot(setarObject = set_1th_l2, boot.scheme = "check", n_digits = 7)
  
  df_comp <- data_frame(orig = lynx, boot = set_1th_l2_b$serie) %>% 
    mutate(n_row = 1:n(),
           th1 = getTh(set_1th_l1_tr)[1],
           th2 = getTh(set_1th_l1_tr)[2],
           reg = regime(set_1th_l1_tr)) 
  
  df_comp %>% 
    filt_diff(tol = 0.01)   
  df_comp %>% 
    qplot(x=n_row, y = as.numeric(orig), data =., geom = "line") +
    geom_point(aes(colour = as.numeric(reg) %>%  factor))
    geom_line(aes(y = boot), colour = "red")
}

## nthresh == 2

### boot
set_2th_l1 <-  setar(lynx, nthresh=2, m=1)
set_2th_l2 <-  setar(lynx, nthresh=2, m=2)
set_2th_l1_tr <-  setar(lynx, nthresh=2, m=1, include = "trend")


setar.boot.check(set_2th_l1)
setar.boot.check(set_2th_l2)
setar.boot.check(set_2th_l2, n_digits = 1)
setar.boot.check(set_2th_l1_tr)
setar.boot.check(set_2th_l1_tr, n_digits = 1)


################
### tets sim
################

## nthresh ==0
set.seed(123)
innov_1 <-  rnorm(200)
sim_nth0 <- setar.sim(B=0.5, lag=1, nthresh=0, 
                      include ="none",
                      starting= 0.4,
                      innov=innov_1,
                      show.parMat = TRUE)$serie

head(sim_nth0)

## nthresh ==1
Bvals <- c(2.9,-0.4,-0.1, 2, 0.2,0.3)
sim_new <- setar.sim(B=Bvals,lag=2, nthresh=1, Thresh=2, starting=c(2.8,2.2),
                     innov=innov_1, show.parMat = TRUE)$serie

head(sim_new)
