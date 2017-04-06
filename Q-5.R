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

library(tolerance)
  
  
