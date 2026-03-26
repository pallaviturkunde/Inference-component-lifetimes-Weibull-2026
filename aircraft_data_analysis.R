####### Data Set Description
#Times between successive failures of air condition
#(A/C) equipment in a Boeing 720 aircraft
#This data set is first discussed in Proschan 1963,Technometrics then in 
# Pal, Jin and Lim (2005, p. 175) and finally in Bhaumik et al Technometrics 2009
# Load data from repository (relative path)

data <- read.csv("../data/aircraft_data.csv") # import data 

# Extract Time Between Failures 
TBF <- data$TimeBetweenFailure  

summary(TBF)

# Compute failure times (cumulative sum)
x <- cumsum(TBF)
summary(x)

# Number of observations
m <- length(x)
summary(x) 
m = length(x)
m ### size of data
## Log transformation
y = log(x)
summary(y)
### Step-1: Computation of empirical cdf
p <- seq(1,10,0.01)
EmpSysCDF <- rep(0,length(p))
ys <- sort(unique(y))
um <- length(ys)
for(i in 1:um)
{
  if(i==1) S = which(p<=ys[1])
  if(i>1 & i<um) 
  {
    S <- which(p>=ys[i]&p<ys[i+1])
  }
  if(i==um) S <- which(p>=ys[um])
  j <- length(which(y<=ys[i]))
  EmpSysCDF[S] <- j/m
}
summary(EmpSysCDF) ### Empirical CDF 
### Step-2: k-out-of-n: F Model
n = 4; k = 2
pw = c(0,1/4,1/2,1/4) ## Signature of system
### Step-3: Generating values of W
qx <- quantile(x,c(0.25,0.5,0.75))
w = rep(3,m)
k1 = which(x<=qx[1]); k2 = which(x>qx[3])
w[k1] = 2; w[k2] = 4
table(w)/m
### Step-4: Computing MLEs of parameterss
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
sum1 = sum(y/svar)
sum2 = sum(sexp/svar) 
sum3 = sum(y*sexp/svar) 
sum4 = sum(1/svar)
sum5 = sum(sexp^2/svar)
estsig = (sum1*sum2-sum3*sum4)/(sum2^2-sum5*sum4) ### theta.est
estmu = (sum1-estsig*sum2)/sum4  ### zi.est
estmu;estsig
LogLik = function(par)
{
  sum = 0
  for(i in 1:m)
  {
    term1 = (w[i]-1)*log(1-exp(-exp((y[i]-par[1])/par[2])))
    term2 = -(n-w[i]+1)*exp((y[i]-par[1])/par[2])
    term3 = (y[i]-par[1])/par[2]-log(par[2])
    sum = sum + term1 + term2 + term3
  }
  return(-sum)
}
A = constrOptim(c(estmu,estsig),LogLik,grad=NULL,ui=c(0,1),ci=0)
A$par ### MLEs
A$convergence
mumle = A$par[1]
sigmle = A$par[2]
mumle;sigmle
### Step-5: Validation of Model
cmpcdf <- function(t1)
{
  z1 <- (t1-mumle)/sigmle
  term1 <- exp(-exp(z1))
  return(1-term1)
}
Sys.Surv = function(t1)
{
  s1 = 0
  for(i in 1:n)
  {
    s2 = 0
    for(j in 0:(i-1))
    {
      s2 = s2 + choose(n,j)*(cmpcdf(t1))^(j)*(1-cmpcdf(t1))^(n-j)
    }
    s1 = s1 + s2*pw[i]
  }
  return(s1)
}
ssurv <- Sys.Surv(p)
summary(1-ssurv)
plot(p,1-ssurv,xlab = "t",ylab="CDF of System Lifetime",ylim=c(0,1),type="l",lwd=1.2)
lines(p,EmpSysCDF,col="red")
KScal <- max(abs(EmpSysCDF-(1-ssurv)))
KScal
KStab <- 0.24 ## at 5% level of significance
KScal < KStab
