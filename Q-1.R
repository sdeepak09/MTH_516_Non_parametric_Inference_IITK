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
