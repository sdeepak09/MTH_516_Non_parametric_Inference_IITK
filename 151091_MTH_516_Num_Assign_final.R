## Written By Deepak Singh Student IITK India
## Date 05/04/2017
## Numerical Assignemnt MTH 516 2017 submitted to Dr. Subhajit Dutta

## Special thing is in this entire code for lopp is not used 
## Starting Lines 
## Q-1  16
## Q-2  350
## Q-3  504
## Q-4  597
## Q-5  644
## Q-6  843



## Q-1

## Grid for 25 dimension
grid.unif=function(n,dim){
  matrix(runif(n*dim,0,1),ncol = dim)
}

grid.gauss=function(n,dim){
  seq=seq(-3,3,length.out = 5000)
  matrix(sample(seq,n*dim,replace = T),ncol = dim)
}  # Here we may get same point but probability is very Low But again it will not creat problem

grid.exp=function(n,dim){
  seq=seq(0,3,length.out = 5000) ## As we know that 95% point in exp(1) is 2.995732
  matrix(sample(seq,n*dim,replace = T),ncol = dim)
}  # Here we may get same point but probability is very Low But again it will not creat problem

grid.cauchy=function(n,dim){
  seq=seq(-7,7,length.out = 50000) ## As we know that 95% point in cauchy(0,1) is 6.3137
  matrix(sample(seq,n*dim,replace = T),ncol = dim)
}  # Here we may get same point but probability is very Low But again it will not creat problem


## kernel Function
uni.gaussian=function(x){
  (1/sqrt(2*pi))*exp(-(x^2)/2)
}

mult.gaussian=function(x,trn.dt,d,h,H_inv){
  (1/(nrow(trn.dt)*(2*pi)^(d/2)*(h)^(d/2)))*sum(sapply(c(1:nrow(trn.dt)),function(k){exp(-(t(x-trn.dt[k,])%*%H_inv%*%(x-trn.dt[k,]))/2)}))
}

uni.unif.kernel=function(x){
  return(ifelse(x<=1 && x>=0,1,0))
}

prod.ker.unif=function(x){
  prod(sapply(c(1:length(x)),function(j){uni.unif.kernel(x[j])}))
}


dnst.compare=function(est.dnst,true.dnst){
  sapply(c(1:10),function(j){
    mean((est.dnst[[j]]-true.dnst)^2)
  })
}

plot.est.den=function(tst.dt,est.dnst,true.dnst){
  d=as.matrix(data.frame(tst.dt,est.dnst[[1]],est.dnst[[2]],est.dnst[[3]],est.dnst[[4]],est.dnst[[5]],est.dnst[[6]],est.dnst[[7]],est.dnst[[8]],est.dnst[[9]],est.dnst[[10]]))
  d=d[order(d[,1]),]
  plot(d[,1],true.dnst,type="b",col="black",xlab = "x",ylab = "Estimated Density",main = "Plot of estimated densitites w.r.t different h, Gaussian Kernel")
  lines(d[,1],d[,2],type="b",col="red")
  lines(d[,1],d[,3],type="b",col="green")
  lines(d[,1],d[,4],type="b",col="darkred")
  lines(d[,1],d[,5],type="b",col="purple")
  lines(d[,1],d[,6],type="b",col="aquamarine")
  lines(d[,1],d[,7],type="b",col="coral1")
  lines(d[,1],d[,8],type="b",col="burlywood4")
  lines(d[,1],d[,9],type="b",col="darkmagenta")
  lines(d[,1],d[,10],type="b",col="gold")
  lines(d[,1],d[,11],type="b",col="blue")
}


## Univariate Case Using Gaussian KErnel
uni.denst.est.gauss.kernl=function(j,trn.dt,tst.dt){
  h=2^(j)
  dnst=sapply(c(1:length(tst.dt)),function(j){
    sum(sapply(((tst.dt[j]-trn.dt)/h),uni.gaussian))/(length(trn.dt)*h)
  })
  return(dnst)
}
tst.dt.1.1.1=grid.gauss(500,1)
trn.dt.1.1.1=rnorm(100,0,1)
est.dnst.1.1.1=lapply(c(-4:5),uni.denst.est.gauss.kernl,trn.dt=trn.dt.1.1.1,tst.dt=tst.dt.1.1.1)
z=dnst.compare(est.dnst.1.1.1,dnorm((tst.dt.1.1.1)))
plot(dnst.compare(est.dnst.1.1.1,dnorm((tst.dt.1.1.1))))


tst.dt.1.1.2=grid.exp(500,1)
trn.dt.1.1.2=rexp(100,1)
est.dnst.1.1.2=lapply(c(-4:5),uni.denst.est.gauss.kernl,trn.dt=trn.dt.1.1.2,tst.dt=tst.dt.1.1.2)
plot.est.den(tst.dt.1.1.2,est.dnst.1.1.2,dexp(sort(tst.dt.1.1.2)))
z=dnst.compare(est.dnst.1.1.2,dexp((tst.dt.1.1.2)))
plot(dnst.compare(est.dnst.1.1.2,dexp((tst.dt.1.1.2))))

tst.dt.1.1.3=grid.unif(500,1)
trn.dt.1.1.3=runif(100,0,1)
est.dnst.1.1.3=lapply(c(-4:5),uni.denst.est.gauss.kernl,trn.dt=trn.dt.1.1.3,tst.dt=tst.dt.1.1.3)
plot.est.den=function(tst.dt,est.dnst,true.dnst){
  d=as.matrix(data.frame(tst.dt,est.dnst[[1]],est.dnst[[2]],est.dnst[[3]],est.dnst[[4]],est.dnst[[5]],est.dnst[[6]],est.dnst[[7]],est.dnst[[8]],est.dnst[[9]],est.dnst[[10]]))
  d=d[order(d[,1]),]
  plot(d[,1],d[,2],type="b",col="red",xlab = "x",ylab = "Estimated Density",main = "Plot of estimated densitites w.r.t different h, Gaussian Kernel")
  lines(d[,1],d[,3],type="b",col="green")
  lines(d[,1],d[,4],type="b",col="darkred")
  lines(d[,1],d[,5],type="b",col="purple")
  lines(d[,1],d[,6],type="b",col="aquamarine")
  lines(d[,1],d[,7],type="b",col="coral1")
  lines(d[,1],d[,8],type="b",col="burlywood4")
  lines(d[,1],d[,9],type="b",col="darkmagenta")
  lines(d[,1],d[,10],type="b",col="gold")
  lines(d[,1],d[,11],type="b",col="blue")
  lines(d[,1],true.dnst,type="b",col="blue")
  
}
## Problem in plotting is happening
plot.est.den(tst.dt.1.1.3,est.dnst.1.1.3,dunif(sort(tst.dt.1.1.3)))
dnst.compare(est.dnst.1.1.3,dunif((tst.dt.1.1.3)))
plot(dnst.compare(est.dnst.1.1.3,dunif((tst.dt.1.1.3))))


tst.dt.1.1.4=grid.cauchy(500,1)
trn.dt.1.1.4=rcauchy(100,0,1)
est.dnst.1.1.4=lapply(c(-4:5),uni.denst.est.gauss.kernl,trn.dt=trn.dt.1.1.4,tst.dt=tst.dt.1.1.4)
plot.est.den=function(tst.dt,est.dnst,true.dnst){
  d=as.matrix(data.frame(tst.dt,est.dnst[[1]],est.dnst[[2]],est.dnst[[3]],est.dnst[[4]],est.dnst[[5]],est.dnst[[6]],est.dnst[[7]],est.dnst[[8]],est.dnst[[9]],est.dnst[[10]]))
  d=d[order(d[,1]),]
  plot(d[,1],d[,2],type="b",col="red",xlab = "x",ylab = "Estimated Density",main = "Plot of estimated densitites w.r.t different h, Gaussian Kernel")
  lines(d[,1],d[,3],type="b",col="green")
  lines(d[,1],d[,4],type="b",col="darkred")
  lines(d[,1],d[,5],type="b",col="purple")
  lines(d[,1],d[,6],type="b",col="aquamarine")
  lines(d[,1],d[,7],type="b",col="coral1")
  lines(d[,1],d[,8],type="b",col="burlywood4")
  lines(d[,1],d[,9],type="b",col="darkmagenta")
  lines(d[,1],d[,10],type="b",col="gold")
  lines(d[,1],d[,11],type="b",col="blue")
  lines(d[,1],true.dnst,type="b",col="black")
  
}
plot.est.den(tst.dt.1.1.4,est.dnst.1.1.4,dcauchy(sort(tst.dt.1.1.4)))
dnst.compare(est.dnst.1.1.4,dcauchy((tst.dt.1.1.4)))
plot(dnst.compare(est.dnst.1.1.4,dcauchy((tst.dt.1.1.4))))

## Univariate Case Using Uniform KErnel

uni.denst.est.unif.kernl=function(j,trn.dt,tst.dt){
  h=2^(j)
  dnst=sapply(c(1:length(tst.dt)),function(j){
    sum(sapply(((tst.dt[j]-trn.dt)/h),uni.unif.kernel))/(length(trn.dt)*h)
  })
  return(dnst)
}

tst.dt.1.2.0=grid.gauss(500,1)
trn.dt.1.2.0=rnorm(100,0,1)
est.dnst.1.2.0=lapply(c(-4:5),uni.denst.est.unif.kernl,trn.dt=trn.dt.1.2.0,tst.dt=tst.dt.1.2.0)
plot.est.den(tst.dt.1.2.0,est.dnst.1.2.0,dnorm(sort(tst.dt.1.2.0)))
dnst.compare(est.dnst.1.2.0,dnorm((tst.dt.1.2.0)))
plot(dnst.compare(est.dnst.1.2.0,dnorm((tst.dt.1.2.0))))

tst.dt.1.2.1=grid.unif(500,1)
trn.dt.1.2.1=runif(100,0,1)
est.dnst.1.2.1=lapply(c(-4:5),uni.denst.est.unif.kernl,trn.dt=trn.dt.1.2.1,tst.dt=tst.dt.1.2.1)
plot.est.den(tst.dt.1.2.1,est.dnst.1.2.1,dunif(sort(tst.dt.1.2.1)))
dnst.compare(est.dnst.1.2.1,dunif((tst.dt.1.2.1)))
plot(dnst.compare(est.dnst.1.2.1,dunif((tst.dt.1.2.1))))

tst.dt.1.2.2=grid.exp(500,1)
trn.dt.1.2.2=rexp(100,1)
est.dnst.1.2.2=lapply(c(-4:5),uni.denst.est.unif.kernl,trn.dt=trn.dt.1.2.2,tst.dt=tst.dt.1.2.2)
plot.est.den(tst.dt.1.2.2,est.dnst.1.2.2,dexp(sort(tst.dt.1.2.2)))
dnst.compare(est.dnst.1.2.2,dexp((tst.dt.1.2.2)))
plot(dnst.compare(est.dnst.1.2.2,dexp((tst.dt.1.2.2))))


tst.dt.1.2.3=grid.cauchy(500,1)
trn.dt.1.2.3=rcauchy(100,0,1)
est.dnst.1.2.3=lapply(c(-4:5),uni.denst.est.unif.kernl,trn.dt=trn.dt.1.2.3,tst.dt=tst.dt.1.2.3)
plot.est.den(tst.dt.1.2.3,est.dnst.1.2.3,dcauchy(sort(tst.dt.1.2.3)))
dnst.compare(est.dnst.1.2.3,dcauchy((tst.dt.1.2.3)))
plot(dnst.compare(est.dnst.1.2.3,dcauchy((tst.dt.1.2.3))))


## Multidimensional Case USing Multidimensional Gaussian case

mult.denst.est.gauss.kernl=function(j,d,tst.dt,trn.dt){
  h=2^(j)
  H_inv=diag(1/h,d)
  dnst=sapply(c(1:nrow(tst.dt)),function(i){mult.gaussian(tst.dt[i,],trn.dt,d,h,H_inv)})
  return(dnst)
}

library(mvtnorm)
dim=2
tst.dt.1.2.1.2=grid.gauss(500,dim)
trn.dt.1.2.1.2=rmvnorm(100,mean =c(rep(0,dim)),sigma = diag(dim))
est.dnst.1.2.1.2=lapply(c(-6:5),mult.denst.est.gauss.kernl,d=dim,trn.dt=trn.dt.1.2.1.2,tst.dt=tst.dt.1.2.1.2)
dnst.compare(est.dnst.1.2.1.2,dmvnorm(tst.dt.1.2.1.2,mean =c(rep(0,dim)),sigma = diag(1,dim)))
plot(dnst.compare(est.dnst.1.2.1.2,dmvnorm(tst.dt.1.2.1.2,mean =c(rep(0,dim)),sigma = diag(dim))))


dim=5
tst.dt.1.2.1.5=grid.gauss(500,dim)
trn.dt.1.2.1.5=rmvnorm(100,mean =c(rep(0,dim)),sigma = diag(dim))
est.dnst.1.2.1.5=lapply(c(-6:5),mult.denst.est.gauss.kernl,d=dim,trn.dt=trn.dt.1.2.1.5,tst.dt=tst.dt.1.2.1.5)
z=dnst.compare(est.dnst.1.2.1.5,dmvnorm(tst.dt.1.2.1.5,mean =c(rep(0,dim)),sigma = diag(1,dim)))
plot(dnst.compare(est.dnst.1.2.1.5,dmvnorm(tst.dt.1.2.1.5,mean =c(rep(0,dim)),sigma = diag(dim))))

dim=10
tst.dt.1.2.1.10=grid.gauss(500,dim)
trn.dt.1.2.1.10=rmvnorm(100,mean =c(rep(0,dim)),sigma = diag(dim))
est.dnst.1.2.1.10=lapply(c(-6:5),mult.denst.est.gauss.kernl,d=dim,trn.dt=trn.dt.1.2.1.10,tst.dt=tst.dt.1.2.1.10)
dnst.compare(est.dnst.1.2.1.10,dmvnorm(tst.dt.1.2.1.10,mean =c(rep(0,dim)),sigma = diag(1,dim)))
plot(dnst.compare(est.dnst.1.2.1.10,dmvnorm(tst.dt.1.2.1.10,mean =c(rep(0,dim)),sigma = diag(dim))))

dim=25
tst.dt.1.2.1.25=grid.gauss(500,dim)
trn.dt.1.2.1.25=rmvnorm(100,mean =c(rep(0,dim)),sigma = diag(dim))
est.dnst.1.2.1.25=lapply(c(-6:5),mult.denst.est.gauss.kernl,d=dim,trn.dt=trn.dt.1.2.1.25,tst.dt=tst.dt.1.2.1.25)
z=dnst.compare(est.dnst.1.2.1.25,dmvnorm(tst.dt.1.2.1.25,mean =c(rep(0,dim)),sigma = diag(1,dim)))
plot(dnst.compare(est.dnst.1.2.1.25,dmvnorm(tst.dt.1.2.1.25,mean =c(rep(0,dim)),sigma = diag(dim))))


## MUltivariate Case data from Multivariate Cauchy using Multivariate Gaussian
install.packages("LaplacesDemon")
library(LaplacesDemon)
dim=2
tst.dt.2.2.2.2=grid.cauchy(500,dim)
trn.dt.2.2.2.2=rmvc(100,mu=c(rep(0,dim)),S=diag(dim))
est.dnst.2.2.2.2=lapply(c(-6:5),mult.denst.est.gauss.kernl,d=dim,trn.dt=trn.dt.2.2.2.2,tst.dt=tst.dt.2.2.2.2)
z=dnst.compare(est.dnst.2.2.2.2,dmvc(tst.dt.2.2.2.2,mu=c(rep(0,dim)),S=diag(1,dim)))
plot(dnst.compare(est.dnst.2.2.2.2,dmvc(tst.dt.2.2.2.2,mu=c(rep(0,dim)),S=diag(dim))))


dim=5
tst.dt.2.2.2.5=grid.cauchy(500,dim)
trn.dt.2.2.2.5=rmvc(100,mu=c(rep(0,dim)),S=diag(dim))
est.dnst.2.2.2.5=lapply(c(-6:5),mult.denst.est.gauss.kernl,d=dim,trn.dt=trn.dt.2.2.2.5,tst.dt=tst.dt.2.2.2.5)
dnst.compare(est.dnst.2.2.2.5,dmvc(tst.dt.2.2.2.5,mu=c(rep(0,dim)),S=diag(1,dim)))
plot(dnst.compare(est.dnst.2.2.2.5,dmvc(tst.dt.2.2.2.5,mu=c(rep(0,dim)),S=diag(dim))))

dim=10
tst.dt.2.2.2.10=grid.cauchy(500,dim)
trn.dt.2.2.2.10=rmvc(100,mu=c(rep(0,dim)),S=diag(dim))
est.dnst.2.2.2.10=lapply(c(-6:5),mult.denst.est.gauss.kernl,d=dim,trn.dt=trn.dt.2.2.2.10,tst.dt=tst.dt.2.2.2.10)
z=dnst.compare(est.dnst.2.2.2.10,dmvc(tst.dt.2.2.2.10,mu=c(rep(0,dim)),S=diag(1,dim)))
plot(dnst.compare(est.dnst.2.2.2.10,dmvc(tst.dt.2.2.2.10,mu=c(rep(0,dim)),S=diag(dim))))

dim=25
tst.dt.2.2.2.25=grid.cauchy(500,dim)
trn.dt.2.2.2.25=rmvc(100,mu=c(rep(0,dim)),S=diag(dim))
est.dnst.2.2.2.25=lapply(c(-6:5),mult.denst.est.gauss.kernl,d=dim,trn.dt=trn.dt.2.2.2.25,tst.dt=tst.dt.2.2.2.25)
z=dnst.compare(est.dnst.2.2.2.25,dmvc(tst.dt.2.2.2.25,mu=c(rep(0,dim)),S=diag(1,dim)))
plot(dnst.compare(est.dnst.2.2.2.25,dmvc(tst.dt.2.2.2.25,mu=c(rep(0,dim)),S=diag(dim))))


## MUltivariate Case data from Multivariate exponential using Multivariate Gaussian kernel
rmvexp=function(n,dim){
  sapply(1:dim,function(i,n){rexp(n)},n)
}
dmvexp=function(X){
  sapply(1:nrow(X),function(j){
    return(prod(dexp(X[j,])))
  })
}

dim=2
tst.dt.2.2.3.2=grid.exp(500,dim)
trn.dt.2.2.3.2=rmvexp(100,dim)
est.dnst.2.2.3.2=lapply(c(-6:5),mult.denst.est.gauss.kernl,d=dim,trn.dt=trn.dt.2.2.3.2,tst.dt=tst.dt.2.2.3.2)
z=dnst.compare(est.dnst.2.2.3.2,dmvexp(tst.dt.2.2.3.2))
plot(dnst.compare(est.dnst.2.2.3.2,dmvexp(tst.dt.2.2.3.2)))


dim=5
tst.dt.2.2.3.5=grid.exp(500,dim)
trn.dt.2.2.3.5=rmvexp(100,dim)
est.dnst.2.2.3.5=lapply(c(-6:5),mult.denst.est.gauss.kernl,d=dim,trn.dt=trn.dt.2.2.3.5,tst.dt=tst.dt.2.2.3.5)
dnst.compare(est.dnst.2.2.3.5,dmvexp(tst.dt.2.2.3.5))
plot(dnst.compare(est.dnst.2.2.3.5,dmvexp(tst.dt.2.2.3.5)))


dim=10
tst.dt.2.2.3.10=grid.exp(500,dim)
trn.dt.2.2.3.10=rmvexp(100,dim)
est.dnst.2.2.3.10=lapply(c(-6:5),mult.denst.est.gauss.kernl,d=dim,trn.dt=trn.dt.2.2.3.10,tst.dt=tst.dt.2.2.3.10)
dnst.compare(est.dnst.2.2.3.10,dmvexp(tst.dt.2.2.3.10))
plot(dnst.compare(est.dnst.2.2.3.10,dmvexp(tst.dt.2.2.3.10)))

dim=25
tst.dt.2.2.3.25=grid.exp(500,dim)
trn.dt.2.2.3.25=rmvexp(100,dim)
est.dnst.2.2.3.25=lapply(c(-6:5),mult.denst.est.gauss.kernl,d=dim,trn.dt=trn.dt.2.2.3.25,tst.dt=tst.dt.2.2.3.25)
dnst.compare(est.dnst.2.2.3.25,dmvexp(tst.dt.2.2.3.25))
plot(dnst.compare(est.dnst.2.2.3.25,dmvexp(tst.dt.2.2.3.25)))


rmvunif=function(n,dim){
  sapply(1:dim,function(i,n){runif(n)},n)
}
dmvunif=function(X){
  sapply(1:nrow(X),function(j){
    prod(dunif(X[j,]))
  })
}

dim=2
tst.dt.2.2.4.2=grid.unif(500,dim)
trn.dt.2.2.4.2=rmvunif(100,dim)
est.dnst.2.2.4.2=lapply(c(-6:5),mult.denst.est.gauss.kernl,d=dim,trn.dt=trn.dt.2.2.4.2,tst.dt=tst.dt.2.2.4.2)
dnst.compare(est.dnst.2.2.4.2,dmvunif(tst.dt.2.2.4.2))
plot(dnst.compare(est.dnst.2.2.4.2,dmvunif(tst.dt.2.2.4.2)))


dim=5
tst.dt.2.2.4.5=grid.unif(500,dim)
trn.dt.2.2.4.5=rmvunif(100,dim)
est.dnst.2.2.4.5=lapply(c(-6:5),mult.denst.est.gauss.kernl,d=dim,trn.dt=trn.dt.2.2.4.5,tst.dt=tst.dt.2.2.4.5)
dnst.compare(est.dnst.2.2.4.5,dmvunif(tst.dt.2.2.4.5))
plot(dnst.compare(est.dnst.2.2.4.5,dmvunif(tst.dt.2.2.4.5)))


dim=10
tst.dt.2.2.4.10=grid.unif(500,dim)
trn.dt.2.2.4.10=rmvunif(100,dim)
est.dnst.2.2.4.10=lapply(c(-6:5),mult.denst.est.gauss.kernl,d=dim,trn.dt=trn.dt.2.2.4.10,tst.dt=tst.dt.2.2.4.10)
z=dnst.compare(est.dnst.2.2.4.10,dmvunif(tst.dt.2.2.4.10))
plot(dnst.compare(est.dnst.2.2.4.10,dmvunif(tst.dt.2.2.4.10)))


dim=25
tst.dt.2.2.4.25=grid.unif(500,dim)
trn.dt.2.2.4.25=rmvunif(100,dim)
est.dnst.2.2.4.25=lapply(c(-6:5),mult.denst.est.gauss.kernl,d=dim,trn.dt=trn.dt.2.2.4.25,tst.dt=tst.dt.2.2.4.25)
dnst.compare(est.dnst.2.2.4.25,dmvunif(tst.dt.2.2.4.25))
plot(dnst.compare(est.dnst.2.2.4.25,dmvunif(tst.dt.2.2.4.25)))



############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
## Question 2
## NW estimator
NWE_gauss_2D.est=function(X,Y,x,h){
  est.dnst=sapply(c(1:length(X)),function(i){uni.gaussian((X[i]-x)/h)})
  NW_weights=est.dnst/sum(est.dnst)
  est.NWE=Y%*%NW_weights
  return(est.NWE)
}
NWE_gauss_2D.weights=function(X,x,h){
  est.dnst=sapply(c(1:length(X)),function(i){uni.gaussian((X[i]-x)/h)})
  NW_weights=est.dnst/sum(est.dnst)
  return(NW_weights)
}

epsilon=rnorm(200,0,0.1)
X=rnorm(200,0,0.4)
Y_1=1+X+epsilon
Y_2=0.3*X+X^(23)+epsilon
Y_3=exp(3+X)+epsilon
Y_4=sin(X)+epsilon

Y_1_NWE_gauss_2D=NWE_gauss_2D.est(X,Y_1,0.4,0.3)
Y_1_NWE_gauss_2D=sapply(c(1:length(X)),function(j){NWE_gauss_2D.est(X,Y_1,X[j],0.3)})
MSE_NWE_gauss_2D=mean((Y_1-est.Y_1)^2) ## sum of squares of Residuals for NW test
plot(X,Y_1_NWE_gauss_2D,col='Red')
## Simple Linear Regression
## Y=a+b*X   a=mean(Y)-b*mean(X)   b=cov(x,y)/V(x)
b_SLR=((Y_1-mean(Y_1))%*%(X-mean(X)))/sum((X-mean(X))^2)
a_SLR=mean(Y_1)-b_SLR*mean(X)
Y_1_SLR=a_SLR+b_SLR*X
MSE_SLR=mean((Y_1-Y_1_SLR)^2)

## Locally Linear NWE 
LLNWE_gauss_2D=function(X,Y,x,h,b.intl){
  eps=10^-2
  n=1
  wghts=NWE_gauss_2D.weights(X,x,h)
  NWE=sum(Y*wghts)
  a.intl=NWE-b.intl*(sum(wghts*(X-x)))
  repeat{
    b=(sum(wghts*Y*(X-x))-a.intl*(sum(wghts*(X-x))))/(sum(wghts*((X-x)^2)))
    a=NWE-b*(wghts%*%(X-x))
    if(abs(a.intl-a)<eps && abs(b.intl-b)<eps) break
    a.intl=a
    n=n+1
    b.intl=b
  }
  return(list(a=a,b=b,n=n))
}

LLNWE_gauss_2D_new=function(X,Y,x,h){
  X.x=matrix(c(rep(1,length(X)),X-x),nrow = length(X), byrow = F)
  W.x=diag(NWE_gauss_2D.weights(X,x,h))
  l.x.t=rbind(c(1,0))%*%solve(t(X.x)%*%W.x%*%X.x)%*%(t(X.x)%*%W.x)
  r.n.x=sum(l.x.t*Y)
  return(r.n.x)
}
est.Y_1=sapply(c(1:length(X)),function(i){LLNWE_gauss_2D_new(X,Y_1,X[i],0.345)})
mat=cbind(X,Y_1,Y_1_NWE_gauss_2D,est.Y_1)
mat=mat[order(mat[,1]),]
plot(mat[,1],mat[,2],col='black',type = 'b',main="Plots of trues data points and predicted values",ylab = "Y",xlab = "x")
lines(mat[,1],mat[,3],col='red',type = 'b')
lines(mat[,1],mat[,4],col='blue',type = 'b')
line.plot=function(x){
  a_SLR+b_SLR*x
}
curve(line.plot,from = -1,to=1,add = T,col='green')
legend(-1.0, 2.0, legend=c("True data point","NW estimate", "LLNW estimate","SLR"),
       col=c("black","red","blue","green"), lty=c(rep(1,11)),lwd = c(rep(3,11)), cex=0.60)

Y_2=0.3*X+X^(23)+epsilon
Y_2_NWE_gauss_2D=sapply(c(1:length(X)),function(j){NWE_gauss_2D.est(X,Y_2,X[j],0.345)})
b_SLR_Y_2=((Y_2-mean(Y_2))%*%(X-mean(X)))/sum((X-mean(X))^2)
a_SLR_Y_2=mean(Y_2)-b_SLR_Y_2*mean(X)
est.Y_2_LLNWE_gauss_2D=sapply(c(1:length(X)),function(i){LLNWE_gauss_2D_new(X,Y_1,X[i],0.345)})
mat=cbind(X,Y_2,Y_2_NWE_gauss_2D,est.Y_2_LLNWE_gauss_2D)
mat=mat[order(mat[,1]),]
line.plot=function(x){
  a_SLR_Y_2+b_SLR_Y_2*x
}

plot(mat[,1],mat[,2],ylim = c(min(c(mat[,1],mat[,2],mat[,3],line.plot(seq(-1,1,length.out = 50)))),max(c(mat[,1],mat[,2],mat[,3],line.plot(seq(-1,1,length.out = 50))))),col='black',type = 'b',main="Plots of trues data points and predicted values",ylab = "Y",xlab = "x")
lines(mat[,1],mat[,3],col='red',type = 'b')
lines(mat[,1],mat[,4],col='blue',type = 'b')
curve(line.plot,from = -1,to=1,add = T,col='green')
legend(-1.0, 2.0, legend=c("True data point","NW estimate", "LLNW estimate","SLR"),
       col=c("black","red","blue","green"), lty=c(rep(1,11)),lwd = c(rep(3,11)), cex=0.60)



Y_3=exp(3+X)+epsilon
Y_3_NWE_gauss_2D=sapply(c(1:length(X)),function(j){NWE_gauss_2D.est(X,Y_3,X[j],0.3)})
b_SLR_Y_3=((Y_3-mean(Y_3))%*%(X-mean(X)))/sum((X-mean(X))^2)
a_SLR_Y_3=mean(Y_3)-b_SLR_Y_3*mean(X)
est.Y_3_LLNWE_gauss_2D=sapply(c(1:length(X)),function(i){LLNWE_gauss_2D_new(X,Y_1,X[i],0.345)})
mat=cbind(X,Y_3,Y_3_NWE_gauss_2D,est.Y_3_LLNWE_gauss_2D)
line.plot=function(x){
  a_SLR_Y_3+b_SLR_Y_3*x
}
mat=mat[order(mat[,1]),]
plot(mat[,1],mat[,2],ylim = c(min(c(mat[,1],mat[,2],mat[,3],line.plot(seq(-1,1,length.out = 50)))),max(c(mat[,1],mat[,2],mat[,3],line.plot(seq(-1,1,length.out = 50))))),col='black',type = 'b',main="Plots of trues data points and predicted values",ylab = "Y",xlab = "x")
lines(mat[,1],mat[,3],col='red',type = 'b')
lines(mat[,1],mat[,4],col='blue',type = 'b')
curve(line.plot,from = -1,to=1,add = T,col='green')
legend(-1.0, 65, legend=c("True data point","NW estimate", "LLNW estimate","SLR"),
       col=c("black","red","blue","green"), lty=c(rep(1,11)),lwd = c(rep(3,11)), cex=0.60)


Y_4=sin(X)+epsilon
Y_4_NWE_gauss_2D=sapply(c(1:length(X)),function(j){NWE_gauss_2D.est(X,Y_4,X[j],0.3)})
b_SLR_Y_4=((Y_4-mean(Y_4))%*%(X-mean(X)))/sum((X-mean(X))^2)
a_SLR_Y_4=mean(Y_4)-b_SLR_Y_4*mean(X)
est.Y_4_LLNWE_gauss_2D=sapply(c(1:length(X)),function(i){LLNWE_gauss_2D_new(X,Y_1,X[i],0.345)})
mat=cbind(X,Y_4,Y_4_NWE_gauss_2D,est.Y_4_LLNWE_gauss_2D)
mat=mat[order(mat[,1]),]
line.plot=function(x){
  a_SLR_Y_4+b_SLR_Y_4*x
}
plot(mat[,1],mat[,2],ylim = c(min(c(mat[,2],mat[,3],mat[,4],line.plot(seq(-1,1,length.out = 50)))),max(c(mat[,1],mat[,2],mat[,3],line.plot(seq(-1,1,length.out = 50))))),col='black',type = 'b',main="Plots of trues data points and predicted values",ylab = "Y",xlab = "x")
lines(mat[,1],mat[,3],col='red',type = 'b')
lines(mat[,1],mat[,4],col='blue',type = 'b')
curve(line.plot,from = -1,to=1,add = T,col='green')
legend(-1.0, 1.3, legend=c("True data point","NW estimate", "LLNW estimate","SLR"),
       col=c("black","red","blue","green"), lty=c(rep(1,11)),lwd = c(rep(3,11)), cex=0.60)




### NWE for 5.4 data
dt=read.table("http://www.stat.cmu.edu/~larry/all-of-nonpar/=data/motor.dat")
View(dt)
dt=dt[-1,]
dt=as.data.frame(dt)
names(dt)=c("times","accel","strata","v")
dt=dt[,-c(3,4)]
dt[,1]=as.numeric(as.vector(dt[,1]))
dt[,2]=as.numeric(as.vector(dt[,2]))
summary(dt)
X=dt[,1]
Y=dt[,2]
pred.accel=sapply(c(1:length(dt[,1])),function(j){NWE_gauss_2D.est(X,Y,X[j],(95)^(-1/5))})
plot(X,Y,col="red",main = "plot of true acceleration and predicted acc at given time points",ylab = "Accelaration",xlab = "time")
lines(X,pred.accel,col="black",type = "b")
legend(3, -85, legend=c("True data points","estimated accel"),
       col=c("red","black"), lty=c(rep(1,11)),lwd = c(rep(3,11)), cex=0.60)

X_1=seq(0,60,length.out = 10000)
pred.accel=sapply(c(1:length(X_1)),function(j){NWE_gauss_2D.est(X,Y,X_1[j],(95)^(-1/5))})
plot(X_1,pred.accel,col="blue",main = "Plot of predicted values of acceleration w.r.t.time",ylab = "Predicted Acceleration", xlab = "time")

############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
## Q-3
Hodge.lehman.estimator=function(x){
  median(unlist(sapply(c(2:(length(x)+1)),function(i){sapply(c(1:(i-1)),function(j){(x[i]+x[j])/2})})),na.rm = T)
}

sample.mode=function(x){
  brk=seq(min(x),max(x),length.out=50)
  hist=hist(x,brk,plot=F);
  frq=hist$counts
  m=hist$mids
  mode=m[which(frq==max(frq))]
  return(mode[1])
}

## What to do in case of bimodal distribution
## To find out the efficiency Ihave used the information given in the below link
##http://www.emathzone.com/tutorials/basic-statistics/efficiency-of-an-estimator.html

est.normal=lapply(c(1:500),function(i){
  x=rnorm(100,0,1)
  return(list(mean=mean(x),median=median(x),mode=sample.mode(x),HLE=Hodge.lehman.estimator(x)))
})
est.mean.normal=as.vector(unlist(sapply(est.normal,'[','mean')))
est.median.normal=as.vector(unlist(sapply(est.normal,'[','median')))
est.mode.normal=as.vector(unlist(sapply(est.normal,'[','mode')))
est.HLE.normal=as.vector(unlist(sapply(est.normal,'[','HLE')))
mean(est.mean.normal^2)
mean(est.median.normal^2)
mean(est.mode.normal^2)
mean(est.HLE.normal^2)


est.unif=lapply(c(1:500),function(i){
  x=runif(100,-1,1)
  return(list(mean=mean(x),median=median(x),mode=sample.mode(x),HLE=Hodge.lehman.estimator(x)))
})
est.mean.unif=as.vector(unlist(sapply(est.unif,'[','mean')))
est.median.unif=as.vector(unlist(sapply(est.unif,'[','median')))
est.mode.unif=as.vector(unlist(sapply(est.unif,'[','mode')))
est.HLE.unif=as.vector(unlist(sapply(est.unif,'[','HLE')))
mean(est.mean.unif^2)
mean(est.median.unif^2)
mean(est.mode.unif^2)
mean(est.HLE.unif^2)

install.packages("rmutil")
library(rmutil)

est.laplace=lapply(c(1:500),function(i){
  x=rlaplace(100,0,1)
  return(list(mean=mean(x),median=median(x),mode=sample.mode(x),HLE=Hodge.lehman.estimator(x)))
})
est.mean.laplace=as.vector(unlist(sapply(est.laplace,'[','mean')))
est.median.laplace=as.vector(unlist(sapply(est.laplace,'[','median')))
est.mode.laplace=as.vector(unlist(sapply(est.laplace,'[','mode')))
est.HLE.laplace=as.vector(unlist(sapply(est.laplace,'[','HLE')))
mean(est.mean.laplace^2)
mean(est.median.laplace^2)
mean(est.mode.laplace^2)
mean(est.HLE.laplace^2)


est.logistic=lapply(c(1:500),function(i){
  x=rlogis(100,0,1)
  return(list(mean=mean(x),median=median(x),mode=sample.mode(x),HLE=Hodge.lehman.estimator(x)))
})
est.mean.logistic=as.vector(unlist(sapply(est.logistic,'[','mean')))
est.median.logistic=as.vector(unlist(sapply(est.logistic,'[','median')))
est.mode.logistic=as.vector(unlist(sapply(est.logistic,'[','mode')))
est.HLE.logistic=as.vector(unlist(sapply(est.logistic,'[','HLE')))
mean(est.mean.logistic^2)
mean(est.median.logistic^2)
mean(est.mode.logistic^2)
mean(est.HLE.logistic^2)


est.cauchy=lapply(c(1:500),function(i){
  x=rcauchy(100,0,1)
  return(list(mean=mean(x),median=median(x),mode=sample.mode(x),HLE=Hodge.lehman.estimator(x)))
})
est.mean.cauchy=as.vector(unlist(sapply(est.cauchy,'[','mean')))
est.median.cauchy=as.vector(unlist(sapply(est.cauchy,'[','median')))
est.mode.cauchy=as.vector(unlist(sapply(est.cauchy,'[','mode')))
est.HLE.cauchy=as.vector(unlist(sapply(est.cauchy,'[','HLE')))
mean(est.mean.cauchy^2)
mean(est.median.cauchy^2)
mean(est.mode.cauchy^2)
mean(est.HLE.cauchy^2)

############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
## Q-4
install.packages("randtests")
library(randtests)
seq.0.1=function(x){
  sapply(1:length(x),function(j){ifelse((x[j]<median(x)),0,1)})
}
run.count=function(x){
  return(c(1,sapply(2:length(x),function(j){ifelse((x[j]!=x[j-1]),1,0)})))
}
n.C.x=function(n,x){return(factorial(n) / (factorial(x) * factorial(n-x)))}
P.Rn.fun=function(N,n){
  P.Rn.1=ifelse((n==0 | n==N),1,0)
  P.Rn=sapply(2:N,function(j,N,n){ifelse(((j%%2)!=0),{(n.C.x(n-1,(j-1)/2)*n.C.x(N-n-1,((j-1)/2)-1)+n.C.x(N-n-1,(j-1)/2)*n.C.x(n-1,((j-1)/2)-1))/n.C.x(N,n)},{(n.C.x(n-1,(j/2)-1)*n.C.x(N-n-1,(j/2)-1)+n.C.x(N-n-1,(j/2)-1)*n.C.x(n-1,(j/2)-1))/n.C.x(N,n)})},N=N,n=n)
  P.Rn=c(P.Rn.1,P.Rn)
  return(P.Rn)
}
exact.run.test=function(R_N,alpha,p){
  r.1=cumsum(p)
  up=length(which(r.1>(alpha/2)))
  r.2=cumsum(rev(p))
  lo=length(which(r.2>(alpha/2)))
  return(ifelse((R_N>=up | R_N<=lo),"Null Hypothesis not accepted","Null Hypothesis accepted"))
}

data(sweetpotato)
dt=sweetpotato$yield
seq.0.1.s=seq.0.1(dt)
run.cnt=run.count(seq.0.1.s)
R_N=sum(run.cnt)
N=length(dt)
n=length(which(seq.0.1.s == 1))
p=P.Rn.fun(N,n)
exact.run.test(R_N,0.05,p)

#### Various Run tests from randtests Library
test.1=runs.test(dt,alternative="two.sided")
test.2=bartels.rank.test(dt, alternative="two.sided")
test.3=rank.test(dt, alternative="two.sided") 
test.4=turning.point.test(dt, alternative="two.sided") 
test.5=cox.stuart.test(dt, alternative="two.sided")

############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################

## Q-5 

emp.cdf=function(x){
  sapply(x,function(j){sum(x<=j)})/length(x)
}
Kol.simn=function(x){
  ks.stat=max(abs(emp.cdf(x)-pnorm(x,0,1)))
  tab.ks.stat=sqrt((-0.5)*log(0.05/2,exp(1))/100)
  ifelse(ks.stat>=tab.ks.stat,1,0) ## 1 means rejection of H0 and 0 means accept H0
}
## Cramer Van Misses test
cramer.v.mises.stat=function(x){
  x=sort(x)
  cvm=function(y,n){
    (pnorm(x[y],0,1)-(2*y-1)/(2*n))^2
  }
  cvm.stat=(1/(12*length(x)))+sum(sapply(c(1:length(x)),function(i,n){cvm(i,n)},n=length(x)))
  return(cvm.stat)
}
chi.sqr=function(x){
  h.decomp=hist(x,6,plot=F)
  cnts=h.decomp$counts
  brk=h.decomp$breaks
  if(sum(cnts==0)>=1){
    z=which(cnts==0)
    cnts=cnts[-z]
    brk=brk[-c(z+1)]
  }
  chi.sqr.stat=sum(sapply(c(1:length(cnts)),function(k){((length(x)*(pnorm(brk[k+1])-pnorm(brk[k]))-cnts[k])^2)/cnts[k]}))
  chi.sqr.tabl=qchisq(0.95,length(cnts)-1)
  ifelse(chi.sqr.stat>=chi.sqr.tabl,1,0)
}

## Power=P(reject H0|H1 is true)
## In this case power will be mean
mean(sapply(c(rep(100,500)),function(n){x=rnorm(n,0.05,1)
Kol.simn(x)}))
power_0.05_1=sapply(c(1:100),function(i){mean(sapply(c(rep(100,500)),function(n){x=rnorm(n,0.05,1)
Kol.simn(x)}))})
mean(power_0.05_1)
sd(power_0.05_1)

mean(sapply(c(rep(100,500)),function(n){x=rnorm(n,0,sqrt(2))
Kol.simn(x)}))
power_0_sqrt_2=sapply(c(1:100),function(i){mean(sapply(c(rep(100,500)),function(n){x=rnorm(n,0,sqrt(2))
Kol.simn(x)}))})
mean(power_0_sqrt_2)
sd(power_0_sqrt_2)

mean(sapply(c(rep(100,500)),function(n){x=rnorm(n,0.05,sqrt(2))
Kol.simn(x)}))
power_0.05_sqrt_2=sapply(c(1:100),function(i){mean(sapply(c(rep(100,500)),function(n){x=rnorm(n,0.05,sqrt(2))
Kol.simn(x)}))})
mean(power_0.05_sqrt_2)
sd(power_0.05_sqrt_2)


## Critical Values for Cramer Van Mises statistics using Simulation technique at 95% confidence interval
cvm.pval.simulation=as.numeric(quantile(sapply(c(rep(100,5000)),cramer.v.mises.stat),0.95))

z=sapply(c(rep(100,500)),cramer.v.mises.stat)
power.cvmn=function(){
  z=sapply(c(1:500),function(i){x=rnorm(100,0.05,1)
  cramer.v.mises.stat(x)})
  return(sum(z>=0.46)/length(z))
  
}
power.cvmn()
power.cvmn_0.05_1=sapply(c(1:100),function(j){power.cvmn()})
mean(power.cvmn_0.05_1)
sd(power.cvmn_0.05_1)
power.cvmn=function(){
  z=sapply(c(1:500),function(i){x=rnorm(100,0,sqrt(2))
  cramer.v.mises.stat(x)})
  return(sum(z>=0.46)/length(z))
  
}
power.cvmn()
power.cvmn_0_sqrt_2=sapply(c(1:100),function(j){power.cvmn()})
mean(power.cvmn_0_sqrt_2)
sd(power.cvmn_0_sqrt_2)
power.cvmn=function(){
  z=sapply(c(1:500),function(i){x=rnorm(100,0.05,sqrt(2))
  cramer.v.mises.stat(x)})
  return(sum(z>=0.46)/length(z))
  
}
power.cvmn()
power.cvmn_0.05_sqrt_2=sapply(c(1:100),function(j){power.cvmn()})
mean(power.cvmn_0.05_sqrt_2)
sd(power.cvmn_0.05_sqrt_2)




mean(sapply(c(1:500),function(i){x=rnorm(100,0.05,1)
chi.sqr(x)}))
power.chi.sqr_0.05_1=sapply(c(1:100),function(a){mean(sapply(c(1:500),function(i){x=rnorm(100,0.05,1)
chi.sqr(x)}))})
mean(power.chi.sqr_0.05_1)
sd(power.chi.sqr_0.05_1)
mean(sapply(c(1:500),function(i){x=rnorm(100,0,sqrt(2))
chi.sqr(x)}))
power.chi.sqr_0_sqrt_2=sapply(c(1:100),function(a){mean(sapply(c(1:500),function(i){x=rnorm(100,0,sqrt(2))
chi.sqr(x)}))})
mean(power.chi.sqr_0_sqrt_2)
sd(power.chi.sqr_0_sqrt_2)


mean(sapply(c(1:500),function(i){x=rnorm(100,0.05,sqrt(2))
chi.sqr(x)}))
power.chi.sqr_0.05_sqrt_2=sapply(c(1:100),function(a){mean(sapply(c(1:500),function(i){x=rnorm(100,0.05,sqrt(2))
chi.sqr(x)}))})
mean(power.chi.sqr_0.05_sqrt_2)
sd(power.chi.sqr_0.05_sqrt_2)



## when alternate distribution is DE(0,1)
install.packages("smoothmest")
library(smoothmest)
rdoublex(n,mu=0,lambda=1)

power_DE_0_1=sapply(c(1:100),function(i){mean(sapply(c(rep(100,500)),function(n){x=rdoublex(n,0,1)
Kol.simn(x)}))})
mean(power_DE_0_1)
sd(power_DE_0_1)

power.cvmn=function(){
  z=sapply(c(1:500),function(i){x=rdoublex(100,0,1)
  cramer.v.mises.stat(x)})
  return(sum(z>=0.46)/length(z))
  
}
power.cvmn_DE_0_1=sapply(c(1:100),function(j){power.cvmn()})
mean(power.cvmn_DE_0_1)
sd(power.cvmn_DE_0_1)


power.chi.sqr_DE_0_1=sapply(c(1:100),function(a){mean(sapply(c(1:500),function(i){x=rdoublex(100,0,1)
chi.sqr(x)}))})
mean(power.chi.sqr_DE_0_1)
sd(power.chi.sqr_DE_0_1)




#### when alternate distribution is C(0,1)
power_C_0_1=sapply(c(1:100),function(i){mean(sapply(c(rep(100,500)),function(n){x=rcauchy(n,0,1)
Kol.simn(x)}))})
mean(power_C_0_1)
sd(power_C_0_1)

power.cvmn=function(){
  z=sapply(c(1:500),function(i){x=rcauchy(100,0,1)
  cramer.v.mises.stat(x)})
  return(sum(z>=0.46)/length(z))
  
}
power.cvmn_C_0_1=sapply(c(1:100),function(j){power.cvmn()})
mean(power.cvmn_C_0_1)
sd(power.cvmn_C_0_1)


power.chi.sqr_C_0_1=sapply(c(1:100),function(a){mean(sapply(c(1:500),function(i){x=rcauchy(100,0,1)
chi.sqr(x)}))})
mean(power.chi.sqr_C_0_1)
sd(power.chi.sqr_C_0_1)



#### when alternate distribution is EXP(0,1)==exp(1)
power_exp_0_1=sapply(c(1:100),function(i){mean(sapply(c(rep(100,500)),function(n){x=rexp(n,1)
Kol.simn(x)}))})
mean(power_exp_0_1)
sd(power_exp_0_1)

power.cvmn=function(){
  z=sapply(c(1:500),function(i){x=rexp(100,1)
  cramer.v.mises.stat(x)})
  return(sum(z>=0.46)/length(z))
  
}
power.cvmn_exp_0_1=sapply(c(1:100),function(j){power.cvmn()})
mean(power.cvmn_exp_0_1)
sd(power.cvmn_exp_0_1)


power.chi.sqr_exp_0_1=sapply(c(1:100),function(a){mean(sapply(c(1:500),function(i){x=r2exp(n, rate = 1, shift = 0)
chi.sqr(x)}))})
mean(power.chi.sqr_exp_0_1)
sd(power.chi.sqr_exp_0_1)

############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################

## Q-6
emp.cdf.x.y=function(x,y){
  sapply(y,function(j){sum(x<=j)})/length(x)
}
g.n.t=emp.cdf.x.y(runif(50,0,1),t)

y.n.t=function(x,t){
  emp.cdf.x.y(x,t)-t
}


curve(sqrt(50)*(y.n.t(runif(50,0,1),x)),from = 0,to=1, main="Plot of 'sqrt(n)*Y_t' beetween 0 and 1 for n=50")

curve(sqrt(100)*(y.n.t(runif(100,0,1),x)),from = 0,to=1, main="Plot of 'sqrt(n)*Y_t' beetween 0 and 1 for n=100")

curve(sqrt(500)*(y.n.t(runif(500,0,1),x)),from = 0,to=1, main="Plot of 'sqrt(n)*Y_t' beetween 0 and 1 for n=500")

curve(sqrt(1000)*(y.n.t(runif(1000,0,1),x)),from = 0,to=1, main="Plot of 'sqrt(n)*Y_t' beetween 0 and 1 for n=1000")








