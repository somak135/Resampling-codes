#install.packages("beepr")
library(beepr)

generate<-function(n = 1000, g1 = 8, g2 = 0.25, b1 = 0.1, b2 = -0.4, sigma = 0.005) {
  x2 = rbinom(n,1,0.3) #urbanization category (1: urban, 0: rural)
  par = (3*x2) + (5*(1-x2))
  x1  = sapply(par, function(x){1+rpois(1,x)}) #family size
  mean_v = (log(200000)*x2) + (log(120000)*(1-x2))
  y = ceiling(x1/3) * sapply(mean_v, function(x){exp(rnorm(1,mean = x,sd = 0.2))}) #total income
  error = rnorm(n, 0, sigma)
  z = exp(g1 + g2*log(y) + b1*x1 + b2*x2 + error)
  data = data.frame(z, y, x1, x2)
  #summary(data); which(z>y); mean(z)/mean(y)
  return(data)
}

f<-function(x, g0, g2, c) {
  return(log(g0*x+0.2) - g2*log(x) - c)
}

f_prime<-function(x, g0, g2) {
  return((g0/(g0*x+0.2)) - (g2/x))
}

p_line = function(dat) {
  fit = lm(log(z) ~ log(y) + x1 + as.factor(x2), data = dat)
  est = fit$coefficients
  g0 = mean(dat[,1])/mean(dat[,2])
  c = est[1] + (est[3]*4) + (est[4]*1)
  theta = mean(dat[,2])
  repeat
  {
    theta1 = theta - f(theta, g0, est[2], c)/(f_prime(theta, g0, est[2]))
    if(abs(theta1 - theta) < 0.01)
    {
      break
    }
    else
    {
      theta=theta1
    }
  }
  return(theta1)
}

vjack = function(dat)
{
  n=nrow(dat)
  theta=rep(0,n)
  for(i in 1:n)
  {
    dat1=dat[-i,]
    theta[i]=p_line(dat1)
  }
  v=((n-1)^2/n)*var(theta)
  return(v)
}

vboot<-function(dat,B=1000)
{
  n=nrow(dat)
  theta=rep(0,B)
  for(i in 1:B)
  {
    s=sample(1:n,n,replace=TRUE)
    dat1=dat[s,]
    theta[i]=p_line(dat1)
  }
  v=((B-1)/B)*var(theta)
  return(v)
}

N=50
v1=rep(0,N)
v2=rep(0,N)
for(i in 1:N) {
  dat=generate(n=1000,g1=8,g2=0.25,b1=0.1,b2=-0.4,sigma=0.005)
  v1[i]=sqrt(vjack(dat))
  v2[i]=sqrt(vboot(dat))
}


print(paste0("ME-jack: ", mean(v1), ", SD-jack: ", sd(v1), 
             ", ME-boot: ", mean(v2), ", SD-boot: ", sd(v2)))
beep(4)

