#1
library(sandwich)
#read in the data
dat = read.csv("xy.csv")
x = dat$x
y = dat$y
n = length(x)



##############
#Population slope and intercept for simulations
############
beta0 = 1
beta1 = 5


##########
#Simulation for d-h
#########

sigma.error = sqrt(.75)
SD.b1.skew = sigma.error/(sqrt(n-1)*sd(x))

nsim = 10000
b1.skew = rep(0, nsim)
se.skew= rep(0, nsim)
se.skew.hc = rep(0, nsim)
Epsilon.skew = matrix(0, n, nsim)

for(i in 1:nsim)
{
  errors = exp(rnorm(n, 0, sqrt(log(1.5)))) - exp(log(1.5)/2)
  Y = beta0 + beta1*x + errors
  lm.temp = lm(Y~x)
  b1.skew[i] = lm.temp$coef[2]
  se.skew[i] = summary(lm.temp)$coef[2,2]
  se.skew.hc[i] = sqrt(diag(vcovHC(lm.temp, type = "HC2")))[2]
  Epsilon.skew[,i] = errors
}





###########
#Simulation for i-m
###############


Sigma = diag((3.6*sqrt(.75)*abs(x - 1/2)))^2
Xmat = cbind(rep(1,n), x)
SD.b1.het = sqrt((solve(t(Xmat)%*%Xmat)%*%t(Xmat)%*%Sigma%*%Xmat%*%solve(t(Xmat)%*%Xmat))[2,2])

nsim = 10000
b1.het = rep(0, nsim)
se.het = rep(0, nsim)
se.het.hc = rep(0, nsim)
Epsilon.het = matrix(0, n, nsim)
for(i in 1:nsim)
{
  errors = rnorm(n, 0, 3.6*sqrt(.75)*abs(x - 1/2))
  Y = beta0 + beta1*x + errors
  lm.temp = lm(Y~x)
  b1.het[i] = lm.temp$coef[2]
  se.het[i] = summary(lm.temp)$coef[2,2]
  se.het.hc[i] = sqrt(diag(vcovHC(lm.temp, type = "HC2")))[2]
  Epsilon.het[,i] = errors
}










