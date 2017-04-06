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

