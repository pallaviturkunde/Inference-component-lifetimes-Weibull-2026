
# Code for comparison with existing estimators mentioned in the paper
# detail description has given in the paper.
m = 50
mu = 0.5
sig = 1.5
beta = exp(mu)
theta = 1/sig
n = 10
k = 5
#n = 5
#k = 2
#prob = c(0,0.2,0.6,0.2,0) ### bridge structure
prob = c(rep(0,4),1/42,4/42,9/42,14/42,14/42,0) # Con/5/10
v1 = function(s)
{
  value1 = 0
  r = 0.00001
  while(r<5)
  {
    value1 = value1 + log(r)*exp(-(n-s+1)*r)*((1-exp(-r))^(s-1))*0.001
    r = r + 0.001
  }
  return(value1)
}

v2 = function(s)
{
  value2 = 0
  r = 0.00001
  while(r<5)
  {
    value2 = value2 + (log(r))^2*exp(-(n-s+1)*r)*((1-exp(-r))^(s-1))*0.001
    r = r + 0.001
  }
  return(value2)
}

cdf = function(t,a,b)
{
  return(1-exp(-exp((t-a)/b)))
}

compdf = function(t,i,a,b)
{
  term1 = n*choose(n-1,i-1)
  term2 = cdf(t,a,b)^(i-1)*(1-cdf(t,a,b))^(n-i+1)
  term3 = exp((t-a)/b)/b
  return(term1*term2*term3)
}

syspdf = function(t,a,b)
{
  sum = 0
  for(i in k:n)
    sum = sum + compdf(t,i,a,b)*prob[i]
  return(sum)
}


### Likelihood function
LogLik = function(par)
{
  sum = 0
  for(k in 1:m)
    sum = sum + log(syspdf(u[k],par[1],par[2]))
  return(-sum)
}

CondLogLik = function(par)
{
  sum = 0
  for(i in 1:m)
  {
    term1 = (w[i]-1)*log(1-exp(-exp((u[i]-par[1])/par[2])))
    term2 = -(n-w[i]+1)*exp((u[i]-par[1])/par[2])
    term3 = (u[i]-par[1])/par[2]-log(par[2])
    sum = sum + term1 + term2 + term3
  }
  return(-sum)
}

#Simulation of the system lifetimes
estmu = c()
estsig = c()
sim = 5000
mle = matrix(nrow=sim,ncol=2)
con = c() 
condmle = matrix(nrow=sim,ncol=2)
conv = c()

for(simno in 1:sim)
{
  set.seed(simno)
  u = c()
  w = sample(1:n,size=m,prob=prob,replace=T)
  for(i in 1:m)
  {
      z = log(rweibull(n,theta,beta))
      a = sort(z)
      u[i] = a[w[i]]
  }
  sexp = rep(0,m)
  svar = rep(0,m)
  uw = sort(unique(w))
  lw = length(uw)
  cw = tapply(uw,INDEX=1:lw,FUN=v1)*n*choose(n-1,uw-1)
  dw = tapply(uw,INDEX=1:lw,FUN=v2)*n*choose(n-1,uw-1)-cw^2
  for(i in 1:lw)
  {
    a = which(w==uw[i])
    sexp[a] = cw[i]
    svar[a] = dw[i]
  }
  sum1 = sum(u/svar)
  sum2 = sum(sexp/svar) 
  sum3 = sum(u*sexp/svar) 
  sum4 = sum(1/svar)
  sum5 = sum(sexp^2/svar)
  estsig[simno] = (sum1*sum2-sum3*sum4)/(sum2^2-sum5*sum4)
  if(estsig[simno]<0) estsig[simno] = 0.01
  estmu[simno] = (sum1-estsig[simno]*sum2)/sum4
  # initial = c(estmu[simno],estsig[simno])
  # if(initial[2]<=0) initial[2] = 1
  ty = ifelse(estsig[simno]<0.1,0.1,estsig[simno])
  B = constrOptim(c(estmu[simno],ty),CondLogLik,grad=NULL,ui=c(0,1),ci=0)
  condmle[simno,] = B$par
  conv[simno] = B$convergence
  initial1 = c(condmle[simno,1],condmle[simno,2])
  A = constrOptim(initial1,LogLik,grad=NULL,ui=c(0,1),ci=0,outer.iterations=10000)
  mle[simno,] = A$par
  con[simno] = A$convergence
}
l = which(estsig==0.01)
d = length(l)
m11 = mean(estmu)
m21 = mean(condmle[,1])
m31 = mean(mle[,1])
m12 = mean(estsig)
m22 = mean(condmle[,2])
m32 = mean(mle[,2])
mse11 = sum((estmu-mu)^2)/sim
mse21 = sum((condmle[,1]-mu)^2)/sim
mse31 = sum((mle[,1]-mu)^2)/sim
mse12 = sum((estsig-sig)^2)/sim
mse22 = sum((condmle[,2]-sig)^2)/sim
mse32 = sum((mle[,2]-sig)^2)/sim
cat(m11,"",m12,"",m21,"",m22,"",m31,"",m32,"\n")
cat(mse11,"",mse12,"",mse21,"",mse22,"",mse31,"",mse32,"\n")
estbeta = exp(estmu)
esttheta = 1/estsig
NBalabeta = exp(mle[,1])
NBalatheta = 1/mle[,2]
CondMLEbeta = exp(condmle[,1])
CondMLEtheta = 1/condmle[,2]
m41 = mean(estbeta);m42 = mean(esttheta)
mse41 = sum((estbeta-beta)^2)/(sim)
mse42 = sum((esttheta-theta)^2)/(sim)
m51 = mean(CondMLEbeta);m52 = mean(CondMLEtheta)
mse51 = sum((CondMLEbeta-beta)^2)/sim
mse52 = sum((CondMLEtheta-theta)^2)/(sim)
m61 = mean(NBalabeta); m62 = mean(NBalatheta)
mse61 = sum((NBalabeta-beta)^2)/sim
mse62 = sum((NBalatheta-theta)^2)/(sim)
cat(m41,"",m42,"",m51,"",m52,"",m61,"",m62,"\n")
cat(mse41,"",mse42,"",mse51,"",mse52,"",mse61,"",mse62,"\n")

############
