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
