#' US monthly industrial production from Hansen (1999)
#' 
#' This data, used as example in Hansen (1999), contains the US monthly
#' industrial production.
#' 
#' 
#' @name IIPUs
#' @docType data
#' @format A monthly time series of class ts starting in January 1960 and
#' ending in September 1997. Note that the series ends at 1997 and not 1998 as
#' in the paper of Hansen, even if the data was taken from hi site and the
#' graph is exactly the same.
#' @source Hansen (1999) Testing for linearity, Journal of Economic Surveys,
#' Volume 13, Number 5, December 1999 , pp. 551-576(26) available at:
#' \url{https://www.ssc.wisc.edu/~bhansen/papers/cv.htm}
#' @keywords datasets
#' @examples
#' 
#' data(IIPUs)
#' end(IIPUs) #not same date as in the paper
#' 
#' ## Figure 2 in paper (p. 555)
#' plot(IIPUs)
#' 
#' ## Table 1 in paper (p. 556)
#' ar_1 <- linear(IIPUs, m=16)
#' coef(summary(ar_1))[, 1:2]
#' deviance(ar_1)
#' 
#' ## Table 2 in paper (p. 559)
#' set_1 <- setar(IIPUs, m=16, thDelay=5, trim=0.1)
#' ## tsDyn finds another threshold, with a better SSR:
#' getTh(set_1)
#' deviance(set_1)
#' 
#' ## estimate with Hansen threshold:
#' set_1_han <- setar(IIPUs, m=16, thDelay=5, trim=0.1, th = 0.23)
#' deviance(set_1_han)
#' set_1_CO <- coef(summary(set_1_han))[, 1:2]
#' cbind(set_1_CO[1:17,], set_1_CO[18:34,])
#' 
#' 
#' ## Table 4 in paper (p. 561)
#' set_2 <- setar(IIPUs, m=16, thDelay=5, trim=0.1, nthresh =2)
#' getTh(set_2)
#' deviance(set_2)
#' 
#' ##most of the results agree, except constant in the low regime which has opposed signs!
#' set_2_CO <- coef(summary(set_2))[, 1:2]
#' cbind(set_2_CO[1:17,], set_1_CO[18:34,])
#' 
#' #this is obviously a error in Hansen, see:
#' XX<-embed(IIPUs, 17)
#' Y<-XX[,1]
#' X<-XX[,-1]
#' dummyDown<-ifelse(X[,6]<= -2.5, 1,0)
#' sum(dummyDown)
#' M<-cbind(1*dummyDown,X*dummyDown )
#' lm(Y~M-1)
#' 
#' ## setar test (takes long to compute, even with small nboot)
#' \dontrun{
#'   test_1 <- setarTest(IIPUs, m=16, thDelay=5, nboot=10)
#'   #because of the discrepency. test1vs2 does not correspond, test 1vs3 conforms
#'   test_1$Ftests ## compare with Table 2 (F12)  and Table 4 (F13, F23)
#'   summary(test_1)
#'   
#'   setarTest(IIPUs, m=16, thDelay=5, nboot=10, test="2vs3")
#'   #test 2vs3 is also different of the version in the article (27)
#' }
#' 
#' ## results from the test is stored in: setarTest_IIPUs_results
#' data(setarTest_IIPUs_results)
#' 
#' ## Table 5 and 6
#' test_1vs <- setarTest_IIPUs_results$test_1
#' test_1vs
#' 
#' ## Table 7
#' test_2vs <- setarTest_IIPUs_results$test_2
#' test_2vs
#' plot(test_2vs)
"IIPUs"


#' Results from the setarTest, applied on Hansen (1999) data
#' 
#' Saved objects
#' 
#' 
#' @name setarTest_IIPUs_results
#' @docType data
#' @format A list containing output from two tests
#' @keywords datasets
#' @seealso Example in \code{\link{IIPUs}}
#' @examples
#' library(tsDyn)
#' data(IIPUs)
#' B <-  1000
#' \dontrun{
#'   test_1 <- setarTest(IIPUs, m=16, thDelay=5, nboot=B)
#'   test_2 <- setarTest(IIPUs, m=16, thDelay=5, nboot=B, test = "2vs3")
#'}
"setarTest_IIPUs_results"