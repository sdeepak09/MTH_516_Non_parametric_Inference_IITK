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





