

#'Extract threshold(s) coefficient
#'
#'Extract threshold coefficient(s)
#'
#'
#'@aliases getTh
#'@param object object of class \code{setar}, \code{summary.setar},
#'\code{nlVar}
#'@param \dots additional arguments to \code{getTh}
#'@return Threshold value.
#'@author Matthieu Stigler
#'@keywords ts
#'@export
#'@examples
#'
#'set<-setar(lynx, m=3)
#'getTh(set)
#'getTh(summary(set))
#'
getTh<- function (object, ...)  
  UseMethod("getTh")

#' @rdname getTh
#' @export
getTh.default <- function(object, ...){
  # look first just in object
  allth<-object[grep("th",names(object))]
  # look then in coef(object)
  if(length(allth)==0){
    allth<-coef(object)[grep("th",names(coef(object)))]
  }
  if(length(allth)==0) allth <- NULL
  # remove thD if there
  if(any(grepl("thD",names(allth))))
    allth<-allth[-grep("thD",names(allth))]
  return(allth)
}

#' @export
getTh.setar<-function(object,...){
  object<-object$coefficients
  getTh.default(object)
}

#' @export
getTh.lstar<-function(object,...){
  object <-object$coef
  object["th"]
}

#' @export
getTh.summary.setar<-function(object,...){
  object$th
}

#' @export
getTh.nlVar<-function(object,...){
  nth <- object$model.specific$nthresh
  
  if(nth>0){
    th<-object$model.specific$Thresh
    if(length(th)==1){
      names(th)<-"th"
    } else{
      names(th)<-c("th1", "th2")
    }
  } else {
    th <- NULL
  }
  
  return(th)
}


######## Get lag  ########
get_lag <- function (object, ...)  UseMethod("get_lag")
get_lag.nlar <- function (object, ...) as.integer(object$str$m)
get_lag.nlVar <- function (object, ...) as.integer(object$lag)

######## Get nthresh  ########
get_nthresh <- function (object, ...)  UseMethod("get_nthresh")
get_nthresh.nlar <- function (object, ...) as.integer(object$model.specific$nthresh)
get_nthresh.nlVar <- function (object, ...) as.integer(object$model.specific$nthresh)
get_nthresh.lstar <- function (object, ...) 1L


######## Get nVar ########
get_nVar <- function (object, ...)  UseMethod("get_nVar")
get_nVar.nlar <- function (object, ...) 1L
get_nVar.nlVar <- function (object, ...) as.integer(object$k)


######## get series ########
get_series <- function(x) {
  if(inherits(x ,"nlVar")) {
    res <- colnames(x$model)[seq_len(x$k)]
  } else {
    res <- x$str$series
  }
  res
}

######## get orig data ########
get_data_orig <- function (object, as.df = FALSE, ...)  UseMethod("get_data_orig")

get_data_orig.nlar <- function (object, as.df = FALSE, ...) {
  res <- object$str$x
  if(as.df) res <- as.data.frame(res)
  res
}
get_data_orig.nlVar <- function (object, as.df = FALSE, ...) {
  res <- object$model[, seq_len(get_nVar(object))]
  if(as.df) res <- as.data.frame(res)
  res
}





if(FALSE) {
  library(tidyverse)
  library(tsDyn)
  
  ############################
  ### Load data
  ############################
  path_mod_uni <- system.file("inst/testdata/models_univariate.rds", package = "tsDyn")
  if(path_mod_uni=="") path_mod_uni <- system.file("testdata/models_univariate.rds", package = "tsDyn")
  
  path_mod_multi <- system.file("inst/testdata/models_multivariate.rds", package = "tsDyn")
  if(path_mod_multi=="") path_mod_multi <- system.file("testdata/models_multivariate.rds", package = "tsDyn")
  
  models_univariate <- readRDS(path_mod_uni) 
  models_multivariate <- readRDS(path_mod_multi)
  models_all <- bind_rows(models_univariate, 
                          models_multivariate)
  ############################
  ### Test
  ############################
  
  mod_uni_1 <-  models_all$object[[1]]
  mod_multi_1 <-  models_multivariate$object[[1]]
  
  get_data_orig(object=mod_uni_1, as.df = TRUE) %>%  head
  get_data_orig(object=mod_multi_1, as.df = TRUE)%>%  head
  
  
  
  ## on all
  sapply(models_all$object, getTh)
  sapply(models_all$object, get_lag)
  sapply(models_all$object, get_series)
  sapply(models_all$object, get_data_orig) %>% 
    sapply(head,1)

  checks <- models_all %>% 
    mutate(lag_out = map_int(object, get_lag),
           nthresh_out = map(object, get_nthresh),
           nVar_out = map_int(object, get_nVar)) %>% 
    select(-object, -object_vars) %>% 
    as.data.frame()
  
  checks
  
  checks %>% 
    filter(lag !=lag_out)
  
}
