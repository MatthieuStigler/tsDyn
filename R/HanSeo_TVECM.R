HanSeo_TVECM <- function(dat, lag=1, gn=300, bn=300, trim=0.05, boot=1000, UserSpecified_beta=NULL, cov=1, p_ests=1, intercept=TRUE, UserSpecified_gamma=NULL, boot.type=c("FixedReg", "ResBoot")) {
warning("This function may be changed in a future version of tsDyn\n")
##de tar_ci.m je crois
# Lags in VAR beyond EC @
#gn=300;			# number of gridpoints for gamma @
#bn = 300;		# number of gridpoints for beta      @
#trim=0.05;		# trimming percentage for threshold @
#boot=1000;		% number of bootstrap replications 
#		        % set equal to zero to not do testing  @
#coint=1;		# set to 1 to estimate cointegrating vector
#			    % set to 0 to fix cointegrating vector at _cvalue @
#cvalue_=1;		% cointegrating vector, if coint_=0 @
#cov=1;		 	#covariance matrix estimation method
#		        % set to 1 for Eicker-White
#		        % set to 0 for conventional homoskedastic estimator @
#p_ests_=1;		% set to 1 to print estimates, else 0 @
#graph_=0;		% set to 1 to generate graph of nonlinear ECM, else 0 @
#graph_rotate_=0;	% set to 1 to generate rotated graphs 
#		        %(useful for some print jobs, but ackward for screen viewing) @


boot.type<-match.arg(boot.type)

### Organize Data
dat<-as.matrix(dat)
n<-nrow(dat) #n=length(dat(:,1))
k<-lag
y<-dat[(2+k):n,]-dat[(1+k):(n-1),] #y=dat(2+k:n,:)-dat(1+k:n-1,:);
t<-nrow(y) #t=length(y(:,1));
xlag<-dat[(1+k):(n-1),] # X(t-1)
cvalue<-UserSpecified_beta
if(ncol(dat)>2) {warning("Please no more than two equations")}
if(is.null(colnames(dat))==TRUE){colnames(dat)<-paste("Var", c(1:2), sep="")}
colnames(xlag)<-paste(colnames(dat), rep(-1,ncol(dat)))
x<-rep(1,t)		#x=ones(t,1);
#j<-1			#j=1;
#while (j<=k){
 #   x<-cbind(x,dat[(2+k-j):(n-j),]-dat[(1+k-j):(n-1-j),])		#x=[x,dat(2+k-j:n-j,:)-dat(1+k-j:n-1-j,:)];
  #  j<-j+1}
#x: intercept and lags matrix

##ma manière
T<-nrow(dat)
nequ<-ncol(dat)

DeltaX<-t(embed(diff(dat),k+1)[,-(1:nequ)])
if(intercept==TRUE){
DeltaX<-rbind(rep(1,T-k-1), DeltaX)}
x<-t(DeltaX)

########################
### Compute Linear Model
########################
xx<- solve(t(x)%*%x)		#xx=inv(x'*x); 
u<-y-x%*%xx%*%t(x)%*%y		#u=y-x*xx*(x'*y); #Residuals for the first auxiliary regression

###Estimation of the cointegrating parameter
if (is.null(UserSpecified_beta)==TRUE){
  v<-xlag-x%*%xx%*%(t(x)%*%xlag)		#Residuals for the first auxiliary regression  v=xlag-x*xx*(x'*xlag);
    S00<-t(u)%*%u				#S00 of size 2x2		#uu=u'*u;
    S11<-t(v)%*%v				#S11 of size 2x2
    m<-solve(S11)%*%(t(v)%*%u)%*%solve(S00)%*%(t(u)%*%v)	#m2x2 #m=inv(v'*v)*(v'*u)*inv(uu)*(u'*v);
    ve<-eigen(m)$vectors 			#eigenvectors 2x2		#[ve,va]=eig(m);
    va<-eigen(m)$values				#Eigenvalues 2
    maa<-which.max(va)				#Selectioon of the biggest eigenvalue	#   [temp,maa]=max(va);
    h<-ve[,maa]					#Eigenvector of the Biggest eigenvalue		#h=ve(:,maa);
    b0<- -h[2]/h[1]				#Eigenvalue ratio		#b0= -h(2)/h(1);
		}
else{ b0<-UserSpecified_beta}				#b0=cvalue_;


###Estimation of the parameters
w0<-xlag%*%c(1,-b0)				#ECT 	#w0=xlag*[1;-b0];
z0<-cbind(w0,x)					#ALL REGRESSORS	#z0=[w0,x];
kk<-length(z0[1,])				#kk=length(z0(1,:));
zz0<-solve(t(z0)%*%z0)			#	zz0=inv(z0'*z0);
zzz0<-z0%*%zz0				#	zzz0=z0*zz0;
beta0<-t(zzz0)%*%y			#ALL Parameters	beta0=zzz0'*y;
colnames(beta0)<-paste("Equation",colnames(dat))
if(intercept==TRUE){rownames(beta0)<-c("ECM","Intercept",c(paste(rep(colnames(dat),k), -rep(1:k, each=nequ))))}
if(intercept==FALSE) {rownames(beta0)<-c("ECM",c(paste(rep(colnames(dat),k), -rep(1:k, each=nequ))))}

e<-y-z0%*%beta0				#RESIDUALS of the VECM nx2	e=y-z0*beta0;
sige<-(t(e)%*%e)/t			#SIGMA	2x2 sige=e'*e/t;
nlike<-(t/2)*log(det(sige))		#	nlike=(t/2)*log(det(sige));
bic<-nlike+log10(t)*4*(1+k)		#bic #=nlike+log10(t)*4*(1+k); ###BIC
aic<-nlike+2*4*(1+k)			#aic=nlike+2*4*(1+k);		###AIC

###Standard error
b_like<-function(b){
	z<-cbind(xlag%*%c(1,-b),x)		# z=[xlag*[1;-b],x];
	sigma<-(t(y-z%*%(solve(t(z)%*%z)%*%t(z)%*%y))%*%(y-z%*%(solve(t(z)%*%z)%*%t(z)%*%y)))/t#sigma=((y-z*(y/z))'*(y-z*(y'/z')'))/t;
like<-(t/2)*log(det(sigma))}

nlike1<-b_like(b0+0.001) 		#nlike1=b_like(b0+.001);
nlike2<-b_like(b0-0.001)			#nlike2=b_like(b0-.001);
hp<-(nlike1+nlike2-2*nlike)/(0.001^2)		#hp=(nlike1+nlike2-2*nlike)/(.001^2);
seb<-1/sqrt(hp)				#seb=1/sqrt(hp);

###VARIANCe COVARIANCE
if (cov==1){
	lz0<-ncol(z0) 			#Number of regressors (ECT, constant and lags)
	xea1<-z0*(e[,1]%*%matrix(1,ncol=lz0))
	xea2<-z0*(e[,2]%*%matrix(1,ncol=lz0))
	xe<-cbind(xea1,xea2)	# xe=[(z0.*(e(:,1)*ones(1,length(z0(1,:))))),(z0.*(e(:,1)*ones(1,length(z0(1,:)))))];
	#lr<-nrow(zz0)  #=length(zz0[,1])
	#lc<-ncol(zz0)  #=length(zz0[1,])
	am0<-matrix(0, ncol=ncol(zz0), nrow=nrow(zz0))
	m0<-rbind(cbind(zz0, am0),cbind( am0, zz0)) 
		#m0=[zz0,zeros(length(zz0(:,1)),length(zz0(1,:)));zeros(length(zz0(:,1)),length(zz0(1,:))),zz0];
	v<- m0%*%(t(xe)%*%xe)%*%m0		#v=m0*(xe*xe)*m0; taille attendue 2*para(cst +ECM+k*2

	}
#print(paste("dim v", dim(v)))
#print(paste("dim sige", dim(sige)))
#print(paste("dim zz0", dim(zz0)))
if (cov==0){ v<-sige%x%solve(t(z0)%*%z0)}			#v=sige.*zz0; 


se<-sqrt(diag(v));
seprint<-matrix(se, nrow=2, byrow=TRUE)
rownames(seprint)<-paste("Equation",colnames(dat))
if(intercept==TRUE){colnames(seprint)<-c("ECM","Intercept",c(paste(rep(colnames(dat),k), -rep(1:k, each=nequ))))}
if(intercept==FALSE) {colnames(seprint)<-c("ECM",c(paste(rep(colnames(dat),k), -rep(1:k, each=nequ))))}

######Graphical output for the linear vecm
file<-""
if(cov==1)std<-"Eicker-White covariance matrix estimator"
if(cov==0) std<-"standard covariance matrix estimator"
cat("\n############################ \n###Linear VECM estimates (by Johansen MLE)\n############################\n \n",  file=file, append=T)
cat("Cointegrating vector", c(1, -b0))
cat("\nStandard error of the cointegrating value", seb)
cat("\nParameters \n")
print(t(beta0))
cat(paste("\n Standard errors with ", std), "\n")
print(seprint)
cat("\nNegative LL \t", nlike)
cat("\nAIC \t\t", aic)
cat("\nBIC \t\t", bic)

###Set up of the grid
qa<-unique(w0)
q<-sort(qa)
gamma1<-q[round(seq(from=trim, to=(gn-1)*((1-2*trim)/gn)+trim, by=((1-2*trim)/gn))*length(q))]
gamma2<-q[round(seq(from=1/(gn+1), by=1/(gn+1), to=(gn/(gn+1)))*length(q))] #original: 1/(gn+1) moi:1/t

if(is.null(UserSpecified_gamma)==FALSE) {gamma1<-UserSpecified_gamma; gamma2<-UserSpecified_gamma; gn<-1}

betas<-seq(from=b0 -6*seb, by=12*seb/(bn-1), to=b0 +6*seb)
if (is.null(UserSpecified_beta)==FALSE) {betas<-UserSpecified_beta;bn<-1}
################
####Grid search
################
store<-matrix(rep(0, gn*bn),ncol=bn)	#store=zeros(gn,bn)+nlike;

j<-1
while (j<=gn){
  gam<-gamma2[j]
  bj<-1
  while (bj<=bn){
    betab<-betas[bj]
    w<-xlag%*%c(1,-betab)	#ECT
    z<-cbind(w,x)		#All variables: ECT and lag
    d1<-ifelse(w<=gam, 1,0)	#Dummy vector 		#d1=(w<=gam);
    n1<-sum(d1)		#Number of elements of the ECT under the threshold
    if(is.na(n1)==TRUE){n1<-0}
    if (min(n1,t-n1)/t>trim) {
      zj<-cbind(z*c(d1), w)		#zj: Matrix with dummy on lags and intercept, but not on ECT #zj=[(z.*(d1*ones(1,length(z(1,:))))),w];		
      betaj<-solve(t(x)%*%x)%*%t(x)%*%zj        	    	#OLS regression zj=f(x)
      zzj<-zj-x%*%betaj					# Residuals zj=f(x) #zzj=zj-x*xx*(x'*zj);
            # warning off;
            #lastwarn(' ');
      bz<-solve(t(zzj)%*%zzj)%*%t(zzj)%*%u		#OLS regression u=f(zzj) where u: Residuals for the first auxiliary regression inlinear VECM # (u'/zzj')';
      u_zzj<-u-zzj%*%bz   				#Residuals   
		#[mw,idw] = lastwarn;
            #lastwarn(' ');
            #warning on;
            #if (1-mw==' ');
            #    bz=pinv(zzj'*zzj)*(zzj'*u);
            #	end;
	
      store[j,bj]<-(t/2)*log(det((t(u_zzj)%*%u_zzj)/t))#store(j,bj)=(t/2)*log(det((u-zzj*bz)*(u-zzj*bz)/t))
    }
    bj<-bj+1}
  j<-j+1}


m<-min(store)
c<-which(store==min(store), arr.ind=TRUE)[2]		#c<-min(diag(store(m,:)));
r<-which(store==min(store), arr.ind=TRUE)[1]		#r=m(c);

betaNLL<-apply(store,2,FUN=min)
gammaNLL<-apply(store,1,FUN=min)


if(is.null(UserSpecified_beta)==TRUE){
  layout(c(1,2))
  plot(gamma2,gammaNLL, type="l", xlab="Threshold parameter gamma", ylab="Negative Log-Likelihood", main="Grid Search")
  plot(betas,betaNLL, type="l", xlab="Cointegrating parameter beta", ylab="Negative Log-Likelihood", main="Grid Search")}
else plot(gamma2,gammaNLL, type="l", xlab="Threshold parameter gamma", ylab="Negative Log-Likelihood", main="Grid Search")

gammahat<-gamma2[r]
b1<-betas[c]

nlike<-store[r,c]
bic<-nlike+log10(t)*8*(1+k)
aic<-nlike+2*8*(1+k)
w<-xlag%*%c(1,-b1)				# w= ECT with beta estimated by the grid
z<-cbind(w,x)					# z= ECT and ALL regressor

###################################
###Subsampling estimation over/under thresh
##################################
z1<-matrix(z[w<=gammahat], ncol=ncol(z))	#observations under the threshold
y1<-matrix(y[w<=gammahat], ncol=ncol(y))

z2<-matrix(z[w>gammahat], ncol=ncol(z))		#observations over the trehshold
y2<-matrix(y[w>gammahat], ncol=ncol(y))

beta1<-solve(t(z1)%*%z1)%*%t(z1)%*%y1
beta2<-solve(t(z2)%*%z2)%*%t(z2)%*%y2

colnames(beta1)<-paste("Equation",colnames(dat))
if(intercept==TRUE){rownames(beta1)<-c("ECM","Intercept",c(paste(rep(colnames(dat),k), -rep(1:k, each=nequ))))}
if(intercept==FALSE) {rownames(beta1)<-c("ECM",c(paste(rep(colnames(dat),k), -rep(1:k, each=nequ))))}

colnames(beta2)<-paste("Equation",colnames(dat))
if(intercept==TRUE){rownames(beta2)<-c("ECM","Intercept",c(paste(rep(colnames(dat),k), -rep(1:k, each=nequ))))}
if(intercept==FALSE) {rownames(beta2)<-c("ECM",c(paste(rep(colnames(dat),k), -rep(1:k, each=nequ))))}

zz1<-solve(t(z1)%*%z1)
zz2<-solve(t(z2)%*%z2)
e1<-y1-z1%*%beta1
e2<-y2-z2%*%beta2

###Var Cov for the threshold model

if (cov==1){
  lz1<-ncol(z1) 			#Number of regressors (ECT, constant and lags)
  xe1a<-z1*(e1[,1]%*%matrix(c(rep(1,lz1)),ncol=lz1))
  xe1b<-z1*(e1[,2]%*%matrix(c(rep(1,lz1)),ncol=lz1))
  xe1<-cbind(xe1a,xe1b)	
  am1<-matrix(c(rep(0, nrow(zz1)*ncol(zz1))), ncol=ncol(zz1))
  m1<-rbind(cbind(zz1, am1),cbind( am1, zz1)) 
  v1<- m1%*%(t(xe1)%*%xe1)%*%m1		

  lz2<-ncol(z2) 			#Number of regressors (ECT, constant and lags)
  xe2a<-z2*(e2[,1]%*%matrix(c(rep(1,lz2)),ncol=lz2))
  xe2b<-z2*(e2[,2]%*%matrix(c(rep(1,lz2)),ncol=lz2))
  xe2<-cbind(xe2a,xe2b)	
  am2<-matrix(c(rep(0, nrow(zz2)*ncol(zz2))), ncol=ncol(zz2))
  m2<-rbind(cbind(zz2,am2),cbind(am2,zz2)) 
  v2<- m2%*%(t(xe2)%*%xe2)%*%m2	
} #end if cov==1

else{
  sig1<-t(e1)%*%e1/nrow(e1)
  indx<-0
  indy<-0
  for (ii in 1:nrow(sig1)){
    for (jj in 1:ncol(sig1)){
      if (indx==0){
        tempx<-zz1*sig1[ii,jj]
        indx<-1
      } #end if indx==0
      else tempx<-cbind(tempx,zz1*sig1[ii,jj])
    }#end fir jj
    indx<-0
    if (indy==0){
      v1<-tempx
      indy<-1}
    else v1<-rbind(v1,tempx)
  } #end for ii
  
  sig2<-t(e2)%*%e2/nrow(e2)
  indx<-0
  indy<-0
  for (ii in 1:nrow(sig2)){
    for (jj in 1:ncol(sig2)){
      if (indx==0){
        tempx<-zz2*sig2[ii,jj]
        indx<-2}
      else tempx<-cbind(tempx,zz2*sig2[ii,jj])
    } #end for jj
    indx<-0
    if (indy==0){
      v2<-tempx
      indy<-2}
    else v2<-rbind(v2,tempx)
  } #end for ii
} #end else, cov=0


se1<-sqrt(diag(v1))	#v1: residuals?
se2<-sqrt(diag(v2))
rb<-nrow(beta1) 		#length(beta1)
bb<-beta1-beta2		#Differences of the parameters

se1print<-matrix(se1, nrow=2, byrow=TRUE)
rownames(se1print)<-paste("Equation",colnames(dat))
if(intercept==TRUE){colnames(se1print)<-c("ECM","Intercept",c(paste(rep(colnames(dat),k), -rep(1:k, each=nequ))))}
if(intercept==FALSE) {colnames(se1print)<-c("ECM",c(paste(rep(colnames(dat),k), -rep(1:k, each=nequ))))}

se2print<-matrix(se2, nrow=2, byrow=TRUE)
rownames(se2print)<-paste("Equation",colnames(dat))
if(intercept==TRUE){colnames(se2print)<-c("ECM","Intercept",c(paste(rep(colnames(dat),k), -rep(1:k, each=nequ))))}
if(intercept==FALSE) {colnames(se2print)<-c("ECM",c(paste(rep(colnames(dat),k), -rep(1:k, each=nequ))))}

########################
###Graphical outuput
########################

cat("\n############################ \n###Threshold VECM estimates \n############################\n \n")
cat("Threshold estimate \t\t")
print(gammahat)
cat("Cointegrating Vector Estimate: \t", b1)
    #for i=1:length(b1)        fprintf(out,'%f   ',b1(i));
cat("\nNegative Log-Like: \t\t", nlike)
cat("\nAIC \t\t\t\t", aic,"\nBIC \t\t\t\t", bic)

cat("\n###First Regime \n \n Percentage of observations \t")
print(nrow(z1)/nrow(z))
cat("\n Parameters \n")
print(t(beta1))
cat(paste("\n Standard errors with ", std))
cat("\n")
print(se1print)    
#for i=1:kk        fprintf(out,'%f   %f\n',beta1(i,1),se1(i));
#  for i=kk+1:2*kk        fprintf(out,'%f   %f\n',beta1(i-kk,2),se1(i));

cat("\n###Second Regime \n \n Percentage of observations \t")
print(nrow(z2)/nrow(z))
cat("\n Parameters \n")
print(t(beta2))
cat(paste("\n Standard errors with ", std, collapse="\n\n"))
cat("\n")
print(se2print)    
#for i=1:kk        fprintf(out,'%f   %f\n',beta1(i,1),se1(i));
#  for i=kk+1:2*kk        fprintf(out,'%f   %f\n',beta1(i-kk,2),se1(i));

#   for i=1:kk       fprintf(out,'%f   %f\n',beta2(i,1),se2(i));
#   for i=kk+1:2*kk        fprintf(out,'%f   %f\n',beta2(i-kk,2),se2(i));

#############
###Wald tests
##############
# bb<-bb
# rb<-rb
if (k>0){
  bw<-matrix(c(bb[-c(1,2),]), ncol=1)		#Vectorization of the lags parameters
  vr<-matrix(c(c(3:rb),c((rb+3):(2*rb))), ncol=1)
  vv<-v1[vr,vr]+v2[vr,vr]
  ww<-t(bw)%*%solve(vv)%*%bw
} #end if k<0

bbx<-matrix(c(bb[1,]), ncol=2)
wecm<-bbx%*%solve(v1[c(1, rb+1), c(1,rb+1)]+v2[c(1, rb+1), c(1,rb+1)])%*%t(bbx)
###########
###Lm Test
###########
y<-y
gn<-gn
x<-x
w0<-w0
k<-k
trim<-0.05
beta0<-beta0
b0<-b0
t<-length(y[,1])
#gamma1<-gamma1
#gammas<-gamma1

lmtest01<-function(y,x,w0,gammas){
#y: y var, x: intercept and lags matrix, w0: ECT term, gammas: potential thresholds

z0<-cbind(w0,x) 		#z0: ECT and intercept and lags
#warning off;
#lastwarn(' ');
#z0zz=z0*inv(z0'*z0);
#[mw,idw] = lastwarn;
#lastwarn(' ');
#warning on;
#if (1-mw==' ')
 #   z0zz=z0*pinv(z0'*z0);
#end;
z0zz<-z0%*%solve(t(z0)%*%z0)
e<-y-z0zz%*%(t(z0)%*%y)		#residuals from the linear VECM given b0
e1<-e[,1]
e2<-e[,2]
store<-rep(0, gn)
j<-1
while (j<=gn){
  d1<-ifelse(w0<=gammas[j],1,0)	#d1: dummy variable		#=w0<=gammas(j);
  n1<-sum(d1)
  if (min(c(n1,(t-n1)))/t>trim){
    z1<-z0*(d1%*%matrix(c(rep(1,ncol(z0))), ncol=ncol(z0)))		#All regressors but only with obs where w0<g
    z11<-z1-z0zz%*%(t(z0)%*%z1)		#residuals from regression of z1 on z0; z11 of dim 479x6
    zea<-z11*(e1%*%matrix(c(rep(1, ncol(z11))), ncol=ncol(z11)))	#  479x6 (479x1 1x6)
    zeb<-z11*(e2%*%matrix(c(rep(1, ncol(z11))), ncol=ncol(z11)))
    ze<-cbind(zea,zeb) #	[(z11.*(e1*ones(1,length(z11(1,:))))),(z11.*(e2*ones(1,length(z11(1,:)))))];
    v<-t(ze)%*%ze
    z11y<-t(z11)%*%y

    s<-matrix(c(z11y), ncol=1)				#vectorization of the parameter matrix z11y

#warning off;
        #lastwarn(' ');
        #sv=(s'/v')';
        #[mw,idw] = lastwarn;
        #lastwarn(' ');
        #warning on;
        #if (1-mw==' ')
        #     sv=pinv(v)*s;
        #end;
    store[j]<-t(s)%*%solve(t(v)%*%v)%*%t(v)%*%s 	
  } #end of the if	
  j<-j+1
} #end of the whole loop
lm01<-max(store)
lm01
} #end of the function lmtest01

lm01<-lmtest01(y,x,w0,gamma1)
#if p_ests_==1
#    out=fopen('tar_ci.txt','at');
#    fprintf(out,'Lagrange Multipler Threshold Test\n\n Test Statistic:     %f\n',lm01);
#    fclose(out);
#end;

##################################
### Fixed Regressor Bootstrap %
##################################
if (boot>0){
  if(boot.type=="FixedReg"){
    boot<-boot
    e<-e
    fix01<-rep(0,boot);
    r<-1
    while (r<=boot){
      yr<-rnorm(n=t*2,0,1)*e
      fix01[r]<-lmtest01(yr,x,w0,gamma1)
      r<-r+1
    }#end while
    
    pfix01a<-ifelse(fix01>lm01,1,0)
    pfix01<-sum(pfix01a)/length(fix01)#pfix01<-mean(fix01>lm01)
    fix01<-sort(fix01)
    cfix01<-fix01[round(.95*boot)]
  
##################################
###Parametric Bootstrap %
##################################
} else{
    x0<-dat[1:(1+k),] #k first observations of the data 
    mu<-beta0[2,] #Intercept parameters of the linear VECM
    ab<-c(1,-b0)%*%beta0[1,] 		###a A VOIR * ou %*% ????
    boot01<-rep(0,boot)
    if (k>0){
      capgamma<-beta0[3:(2+2*k),] 	#Parameters of the lags in the linear VECM
      temp<-x0[k+1,]-x0[k,]		#First left side delta(x)
      if(k>1){
        for (ii in k:2){
          temp<-rbind(temp,x0[ii,]-x0[ii-1,])
        } #end for		#Next left side delta(x) if k>1
      }#end if k>1
        #ind<-0 for (ii in 1:nrow(temp)){ for (jj in 1:ncol(temp)){if (ind==0){ dx0<-temp[ii,jj] ind<-1} else dx0<-rbind(dx0,temp[ii,jj])}}

      dx0<-matrix(c(t(as.matrix(temp))), ncol=1, byrow=TRUE) #vectorization of k the first values
    } #End if k>0
    r<-1
    while (r<=boot){
      datb<-x0
      if (k>0){dx<-dx0}
      datbi<-datb[k+1,]
      ei<-ceiling(runif(n=t,min=0,max=1)*t)
      eb<-e[ei,]			#Bootstraping the original residuals of linear VECM
      i<-k+2
      while (i<=n){
        print(i)
        browser()
        u<-mu+eb[i-k-1,]+datbi*ab
        if (k>0){u<-u+t(dx)%*%capgamma}
        datbi<-datbi+u
        datb<-rbind(datb,datbi)
        if (k==1){dx<-t(u)}
        if (k>1){dx<-rbind(t(u),matrix(dx[1:2*k-2], ncol=1))}
        i<-i+1
      } #end while i<n
      yr<-datb[2+k:n,]-datb[1+k:n-1,]
      xlagr<-datb[1+k:n-1,]
      xr<-rep(1,t)
      j<-1
      while (j<=k){
        xr<-cbind(xr,(datb[2+k-j:n-j,]-datb[1+k-j:n-1-j,]))
        j<-j+1
      }#end while j<k
      if (is.null(UserSpecified_beta)==TRUE){
        xxr<-solve(t(xr)%*%xr)
        u<-yr-xr%*%xxr%*%(t(xr)%*%xlagr)
        v<-xlagr-xr%*%xxr%*%(t(xr)%*%xlagr)
        m<-solve(t(v)%*%v)%*%(t(v)%*%u)%*%solve(t(u)%*%u)%*%(t(u)%*%v)
        ve<-eigen(m)$vectors
        va<-eigen(m)$values			#[ve,va]=eig(m);
        
        maxindva<-which(va==max, arr.ind=TRUE)		#[temp,maxindva]=max(va);
        
        h<-ve[,maxindva]
        wr<-xlagr%*%c(1,(h[2]/h[1]))
      }#end if coint==1
      else{ wr<-xlagr%*%c(1-UserSpecified_beta)}
      boot01[r]<-lmtest01(yr,xr,wr,gamma1)
      
#cat("wr")
#print(wr)
#cat("xr")
#print(xr)

#cat("yr")
#print(yr)}
    }#end of the while

    boot01<-sort(boot01)
    cb01<-boot01[round(.95*boot)]
    pb01a<-ifelse(boot01>lm01,1,0)				#pb01=mean(boot01>lm01);
    pb01<-sum(pb01a)/length(boot01)
  }#end if boot= ResBoot
}#end if boot>0
###########
###Output test

if (k>0){cat("\nWald Test for Equality of Dynamic Coefs: \t")
         print(round(c(ww,1-pchisq(ww,(rb-2)*2)),6))}
cat("\nWald Test for Equality of ECM Coef: \t")
print(round(c(wecm,1-pchisq(wecm,2)),6))

cat("\n###Lm Test", lm01)

if(boot>0){
  if(boot.type=="FixedReg"){
    cat("\n Fixed Regressor (Asymptotic) .05 Critical Value\t", cfix01)
    cat("\n Fixed Regressor (Asymptotic) P-Value:\t", pfix01, "\n")
  }else{
    cat("\n Bootstrap .05 Critical Value:\t", cb01)
    cat("\n Bootstrap P-Value: \t", pb01, "\n")
  }
} #end if boot
}#End of the whole function

###########################
### ENd of function
###########################
if(FALSE) {#usage example
###Test

#dat<-zeroyld[,7:62]
#dat<-dat[,c(30,13)]
#colnames(datax)<-c("short-run", "long-run")
#dat1<-zeroyld[,c(23, 49)]
#datax<-zeroyld[,c(36,19)]
#write.table(datax, file="/home/mat/zeroyld.txt")
#rs<-matrix(c(c(0:18), c(21,24,30), c(seq(from=63, to=120, by=12))))
#rs
#short<-12
#long<-120
#short_i<-13
#long_1<-30
#dat<-dat[,c(30,13)]
#dat
#dim(dat)
data(zeroyld)
data<-zeroyld


HanSeo_TVECM(data, lag=1, intercept=TRUE, UserSpecified_beta=0.98, cov=1, UserSpecified_gamma=-0.63, bn=30, gn=30, boot=0)

HanSeo_TVECM(data, lag=1, intercept=TRUE, cov=1, bn=30, gn=300, boot=0)


}

