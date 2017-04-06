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
