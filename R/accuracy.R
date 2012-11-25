

##########
##### acc_stats: function to get accuracy for true/fitted
##########

accuracy_stat <- function(object, ...)
  UseMethod("accuracy_stat")

accuracy_stat.default <- function(object, true, ...) accuracy_stat_simple(fit=object, true=true)

accuracy_stat.pred_roll <- function(object, w, ...) {

  is_multiH <- "n.ahead" %in% colnames(object$pred)

  if(is_multiH){
    n.aheads <- unique(object$pred[,"n.ahead"])
    nvar <- ncol(object$true)
    li<-list()
    for(i in 1:length(n.aheads)){
      li[[i]] <- accuracy_stat_simple(fit=subset(object$pred, n.ahead==n.aheads[i])[,1:nvar], true=object$true)
    }
    res_raw <- simplify2df(li)
    res_raw <- data.frame(var=rownames(res_raw),res_raw, row.names=1:nrow(res_raw))

  ## add means: 
    if(missing(w)) w <- rep(1/length(n.aheads), length(n.aheads))
    means <- aggregate(res_raw[,2:6], list(res_raw$var), weighted.mean, w=w)
    colnames(means)[1] <- "var"
    res_withmeans <- rbind(res_raw, means)

  ## add horizont column:
    res_withmeans[,"n.ahead"] <- rep(c(n.aheads,"all"), each=nvar+1)
    res <- res_withmeans
    res <- res[order(res$var, res$n.ahead),]


  } else {
    res <- accuracy_stat_simple(fit=object$pred, true=object$true)
  }

return(res)

}
  



accuracy_stat_simple <- function(fit, true){

  isMulti <- !is.null(dim(true))&&ncol(true)>1
  isSameLength <- if(isMulti) !all(dim(true)==dim(fit)) else if(is.null(dim(true)))  length(true)!=length(fit) else nrow(true)!=nrow(fit)
  if(isSameLength) stop("'true' and 'fit' should be of same dimension")

  res<- true -fit
  pe <- res/true * 100 # Percentage error

  res_mat <- as.matrix(res)
  pe_mat <- as.matrix(pe)

  out <- c(mean(res_mat,na.rm=TRUE), sqrt(mean(res_mat^2,na.rm=TRUE)), mean(abs(res_mat),na.rm=TRUE), mean(pe_mat,na.rm=TRUE), mean(abs(pe_mat),na.rm=TRUE))
  names(out) <- c("ME","RMSE","MAE","MPE","MAPE")


  if(isMulti){
    out2 <- t(rbind(colMeans(res,na.rm=TRUE), sqrt(colMeans(res^2,na.rm=TRUE)), colMeans(abs(res),na.rm=TRUE), colMeans(pe,na.rm=TRUE), colMeans(abs(pe),na.rm=TRUE)))
    colnames(out2) <- c("ME","RMSE","MAE","MPE","MAPE")
    out <- rbind(out2, out)
    if(all(rownames(out)[-nrow(out)]=="")) rownames(out)[-nrow(out)] <- paste("Var", 1:(nrow(out)-1), sep="")
    rownames(out)[nrow(out)] <- "global"
  } else {
    out <- t(as.data.frame(out))
    rownames(out) <- ifelse(is.null(names(true)),  "Var1", names(true))
  }

  return(out)

}

simplify2df <- function(x, res=c("matrix", "df")) {
  res<- match.arg(res)

  out <- x[[1]]
  for(i in 2:length(x)){
    out <- rbind(out, x[[i]])
  }
  if(res=="df") out <- as.data.frame(out)
  out
}

if(FALSE){
library(tsDyn)


## univariate
mod_ar <- linear(lynx[1:100], m=1)
mod_ar_pred <- predict_rolling(mod_ar, newdata=lynx[101:114])
accuracy_stat(object=mod_ar_pred$pred, true=mod_ar_pred$true)
accuracy_stat(object=as.matrix(mod_ar_pred$pred), true=as.matrix(mod_ar_pred$true))


## multivariate
data(barry)
mod_var <- lineVar(barry, lag=1)

mod_var_pred <-predict_rolling(object=mod_var, nroll=10, n.ahead=1:3)
accuracy_stat(object=mod_var_pred)
accuracy_stat(object=mod_var_pred, w=c(0.7, 0.2, 0.1))


accuracy_stat(object=mod_var_pred$pred, true=mod_var_pred$true)

accuracy_stat(object=as.matrix(mod_var_pred$pred), true=as.matrix(mod_var_pred$true))
accuracy_stat(object=mod_var_pred$pred[,1], true=mod_var_pred$true[,1])

mod_var_pred_multih <-predict_rolling(object=mod_var, nroll=10, n.ahead=1:2)
accuracy_stat(object=mod_var_pred_multih)
accuracy_stat(object=mod_var_pred_multih, w=c(1,0))

}