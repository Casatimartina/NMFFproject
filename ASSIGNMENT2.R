# ASSIGNEMENT 2 FOR "NUMERICAL METHODS FOR FINANCE" computed by
# Arash Sadid
# Martina Casati
# Andrea Frattini

# EXERCISE 1

library(readxl)
dataset <- read_excel("/Users/andreafrattini/Desktop/University/Numerical Methods for Finance/Assignments/Dataset.xls")
# add here your own path to the dataset.xls file
df <- data.frame(dataset)

Maturity <- 51
TimeToMat <- Maturity/365

mertConst <- function(PriceOpt, S, K, r, TimeToMat, TypeOpt){
  dimens <-length(PriceOpt)
  resCheck <- logical(dimens)
  for(i in c(1:length(PriceOpt))){
    if(TypeOpt[i]=='call'){
      condLow <- PriceOpt[i] >= max(c(S[i] - K[i] * exp(-r[i] * (TimeToMat)), 0))
      condUpp <- PriceOpt[i] <= S[i]
      resCheck[i] <- condLow & condUpp
    }else{
      condLow <- PriceOpt[i] >= max(c(K[i] * exp(-r[i] * TimeToMat)-S[i], 0))
      condUpp <- PriceOpt[i] <= K[i]*exp(-r[i] * TimeToMat)
      resCheck[i] <- condLow & condUpp
    }
  }
  return(resCheck)
}

check <- mertConst(PriceOpt=df$Price.Option, S=df$Underlyng, K=df$Strike, r=df$RiskFreeRate_ccr, 
                   TimeToMat=TimeToMat, TypeOpt = df$Type.Option)
View(check) #All the options satisfy the Merton's constraint, so I can move on

OutOfMoney=numeric(nrow(dataset))
for (i in 1:nrow(dataset)) {
  if (df$Type.Option[i]=="call") {
    OutOfMoney[i]= df$Strike[i]>df$Underlyng[i]
  } else {
    OutOfMoney[i]=df$Strike[i]<df$Underlyng[i]
  }
}
df2=df[which(OutOfMoney==1), ]

# Black & Scholes formulas
TypeOption="put"
if (TypeOption=='call'){
  fun <- function(S0,K,sigma,r,PriceOpt,TimeToMat){
    d1 <- (log(S0/K)+(r+0.5*sigma^2)*(TimeToMat))/(sigma*sqrt(TimeToMat))
    d2 <- d1-sigma*sqrt(TimeToMat)
    error <- S0*pnorm(d1)-K*exp(-r*TimeToMat)*pnorm(d2)- PriceOpt
    return(error^2)
  }
  pos= which(df2[,"Type.Option"]==TypeOption)
  K=df2[pos,"Strike"]
  Price=df2[pos,"Price.Option"]
  implied.vol=numeric(length(pos))
  convergence=numeric(length(pos))
}else{
  fun <- function(S0,K,sigma,r,PriceOpt,TimeToMat){
    d1 <- (log(S0/K)+(r+0.5*sigma^2)*(TimeToMat))/(sigma*sqrt(TimeToMat))
    d2 <- d1-sigma*sqrt(TimeToMat)
    error <- K*exp(-r*TimeToMat)*pnorm(-d2) - S0*pnorm(-d1)-PriceOpt
    return(error^2)
  }
  pos= which(df2[,"Type.Option"]==TypeOption)
  K=df2[pos,"Strike"]
  Price=df2[pos,"Price.Option"]
  implied.vol=numeric(length(pos))
  convergence=numeric(length(pos))
}

S0=52.67
for(i in c(1:length(pos))){
  implied.vol[i] <- optim(par=0.1, fn=fun, lower=0, method ="L-BFGS-B",
                          S0=S0 ,K=df2$Strike[i] , TimeToMat=TimeToMat, 
                          r=df2$RiskFreeRate_ccr[i] , PriceOpt=df2$Price.Option[i])$par
}

plot(K, implied.vol, main = paste0("Volatility Surface ", TypeOption))

if (TypeOption=="call") {
  fun <- function(S, K, r, sigma, dt, Market.Price) {   # CALL BS FORMULA  
    d1=(log(S/K)+(r+0.5*sigma^2)*(dt))/(sigma*sqrt(dt))
    d2=d1-sigma*sqrt(dt)
    error=S*pnorm(d1)-K*exp(1)^(-r*dt)*pnorm(d2)-Market.Price
    return(mean(objective^2))
  }
  
  pos=which(MyTable.OutMoney[,"type"]==TypeOption)
  K=MyTable.OutMoney[pos,"K"]
  Price=MyTable.OutMoney[pos,"Price"]
  
} else {
  objective.fun=function(S, K, r, sigma, dt, Market.Price) {   # PUT BS FORMULA  
    d1=(log(S/K)+(r+0.5*sigma^2)*(dt))/(sigma*sqrt(dt))
    d2=d1-sigma*sqrt(dt)
    objective=K*exp(1)^(-r*dt)*pnorm(-d2)-S*pnorm(-d1)-Market.Price
    return(mean(objective^2))
  }
  
  pos=which(MyTable.OutMoney[,"type"]==TypeOption)
  K=MyTable.OutMoney[pos,"K"]
  Price=MyTable.OutMoney[pos,"Price"]
  
}



calibration=optim(par=0.1, objective.fun, S=S, K=K, dt=residualDays/252, r=r, Market.Price=Price, method = "L-BFGS-B", lower = 0, upper = 100)$par
calibration

# EXERCISE 2
#Inputs
X0=0
alpha=0.2
mu=0.05
sigma=0.2
N=1000 
M=10000
t0=0
FinalT=1
# Point 1
a=0
b=0.7
standDevOfXT<-sigma/sqrt(2*alpha)*sqrt(1-exp(-2*alpha*FinalT-t0))
d_bar1 <- (a-X0*exp(-alpha*(FinalT-t0))-mu*(1-exp(-alpha*(FinalT-t0))))/standDevOfXT
d_bar2 <- (b-X0*exp(-alpha*(FinalT-t0))-mu*(1-exp(-alpha*(FinalT-t0))))/standDevOfXT
ExactProb <- 1+pnorm(d_bar1)-pnorm(d_bar2)
ExactProb
# Point 2
# Euler simulation scheme
set.seed(1)
samplePath<-matrix(0,nrow=M,ncol=N+1)
gridtime<-seq(t0,FinalT, length.out=N+1)
Delta<- (FinalT-t0)/(N)
samplePath[,1]<-X0
for(pos in c(2: (N+1))){
  samplePath[,pos] <-samplePath[,pos-1]+alpha*(mu-samplePath[,pos-1])*Delta+sigma*sqrt(Delta)*rnorm(M) 
}
colnames(samplePath)<-paste0("t = ",gridtime)
#Inputs
X0=0
alpha=0.2
mu=0.05
sigma=0.2
N=10
M=10000
t0=0
FinalT=1

# Exact simulation scheme
set.seed(1)
samplePathExact<-matrix(0,nrow=M,ncol=N+1)
gridtime<-seq(t0,FinalT, length.out=N+1)
samplePathExact[,1]<-X0
Delta<- (FinalT-t0)/(N)
for(pos in c(2:(N+1))){
  samplePathExact[,pos]<- samplePathExact[,pos-1]*exp(-alpha*(Delta))+mu*(1-exp(-alpha*(Delta)))+sigma/sqrt(2*alpha)*sqrt(1-exp(-2*alpha*Delta))*rnorm(M)
}
colnames(samplePathExact)<-paste0("t = ",gridtime)

# MC with Euler simulation Scheme
EulerXT<-samplePath[,dim(samplePath)[2]]
MCProbEuler1 <- mean(EulerXT<a)
MCProbEuler2 <- mean(EulerXT>b)
MCProbEuler <- MCProbEuler1+MCProbEuler2
#CONFIDENCE INTERVAL
UBa <- MCProbEuler1 + 1.96*standDevOfXT/sqrt(N+1)
UBb <- MCProbEuler2 + 1.96*standDevOfXT/sqrt(N+1)
UB <- UBa+UBb

LBa <- MCProbEuler1 - 1.96*standDevOfXT/sqrt(N+1)
LBb <- MCProbEuler2 - 1.96*standDevOfXT/sqrt(N+1)
LB <- LBa+LBb 

CI <- c(LB,MCProbEuler,UB)
CI
# MC with Exact simulation Scheme
ExactXT <- samplePathExact[,dim(samplePathExact)[2]]
MCProbExact1 <- mean(ExactXT<a)
MCProbExact2 <- mean(ExactXT>b)
MCProbExact <- MCProbExact1+MCProbExact2

ExactProb
MCProbEuler
MCProbExact

#CONFIDENCE INTERVAL 
UBExacta <- MCProbExact1 + 1.96*standDevOfXT/sqrt(N+1)
UBExactb <- MCProbExact2 + 1.96*standDevOfXT/sqrt(N+1)
UBExact <- UBExacta+UBExactb 

LBExacta <- MCProbExact1 - 1.96*standDevOfXT/sqrt(N+1)
LBExactb <- MCProbExact2 - 1.96*standDevOfXT/sqrt(N+1)
LBExact <- LBExacta+LBExactb 

CIExact <- c(LBExact,MCProbExact,ub_exact)
CIExact
# Inputs
K=100
S0=100
sigma=0.2
r=0.04 # on yearly basis
TimeToMat <- 90/252 # 3 months on yearly basis
N <- 1000000
seed<-1

BS <- function(S,K, TimeToMat, Sigma, Rfree){
  d1<- (log(S/K)+(Rfree+0.5*Sigma^2)*(TimeToMat))/(Sigma*sqrt(TimeToMat))
  d2 <- d1-Sigma*sqrt(TimeToMat)
  res<- S*pnorm(d1)-K*exp(-Rfree*TimeToMat)*pnorm(d2)
  return(res)
}

MCprice <-function(K, S0, sig, r, TimeToMat, N, seed){
  if(!is.null(seed)){
    set.seed(seed)
  }
  Z<-rnorm(n=N)
  S_T<-S0*exp((r-0.5*sig^2)*TimeToMat+sig*sqrt(TimeToMat)*Z)
  
  DiscountedFinalPayoff <- exp(-r*(TimeToMat))*pmax(S_T-K,0)
  g_n = mean(DiscountedFinalPayoff)
  TrueValue<- BS(S0,K, TimeToMat, Sigma=sigma, Rfree=r)
  S2<- var(DiscountedFinalPayoff)
  UB<- TrueValue+1.96*sqrt(S2)/sqrt(N)
  LB<- TrueValue-1.96*sqrt(S2)/sqrt(N)
  res <- list("N"=N,"Discounted Final Payoff"= g_n,"Upper Bound"= UB,"Lower Bound" = LB,
              "Difference between Ub and LB" = UB-LB,"True Value"= TrueValue)
  return(res)
}

result <- MCprice(K=K, S0=S0, sig=sigma, r=r, TimeToMat=TimeToMat, N=N, seed=seed)
result

# EXERCISE 3
# POINT 1
library(quantmod)
getSymbols("^VIX",  from="2019-01-02", to="2020-12-31")
Price <- na.omit(VIX[,4]) 
daysInOneYear<-365
timeOfObervation <-as.numeric(index(Price))
deltati<-diff(timeOfObervation)/daysInOneYear
length(deltati)
logprice <- log(as.numeric(Price)) 
logret <- diff(logprice)
logprice_old <- log(as.numeric(Price))[-length(as.numeric(Price))]
length(logret)
length(logprice_old)


minuslogLik<- function(par,logprice_old,logret,deltati){
  alpha<-par[1]
  mu<-par[2]
  sig<- par[3]
  vecMean <-  logprice_old*exp(-alpha*deltati)+mu*(1-exp(-alpha*deltati))-logprice_old
  vecsd <- sig/sqrt(2*alpha)*sqrt(1-exp(-2*alpha*deltati))
  -sum(dnorm(logret,mean=vecMean,sd=vecsd,log=TRUE))
}

minuslogLik(par=c(0.1,0.1,0.1),logprice_old=logprice_old,logret=logret,deltati=deltati)

res <- optim(par=c(0.1,0.1,0.1),fn = minuslogLik,method ="L-BFGS-B", 
             lower=c(0.00001, -Inf, 0.0000001), logprice_old=logprice_old,
             logret=logret,deltati=deltati)
res$par
alpha <- res$par[1]
mu <- res$par[2]
sig <- res$par[3]
-1*res$value
res$convergence

# POINT 2
deltati<-diff(timeOfObervation)/daysInOneYear
logret <- diff(log(as.numeric(Price)))
logprice_old <- log(as.numeric(Price))[-length(as.numeric(Price))]

minuslogLik1<- function(par,logprice,logprice_old,deltati){
  alpha<-par[1]
  mu<-par[2]
  sig<- par[3]
  vecMean <-  logprice_old*exp(-alpha*deltati)+mu*(1-exp(-alpha*deltati))
  vecsd <- sqrt((sig^2)/(2*alpha)*(1-exp(-2*alpha*deltati)))
  -sum(dnorm(logprice,mean=vecMean,sd=vecsd,log=TRUE))
}

minuslogLik1(par=c(1,0.5,0.2),logprice=logprice,logprice_old=logprice_old,deltati=deltati)

res <- optim(par=c(1,0.5,0.2),fn = minuslogLik1,method ="L-BFGS-B",
             logprice=logprice, logprice_old=logprice_old,deltati=deltati)
res$par
-1*res$value
res$convergence

# POINT 3

minusQloglik_Vasicek <- function(par,logret,logprice_old,deltati){
  alpha <- par[1]
  mu <- par[2]
  sig <- par[3]
  logdens <- dnorm(logret, mean = alpha*(mu-logprice_old)*deltati, sd= sig*sqrt(deltati),log=T )
  -sum(logdens)
}

minusQloglik_Vasicek(par=c(0.7,0.5,0.2),logret=logret,
                     logprice_old=logprice_old,deltati=deltati)

res <- optim(par=c(0.7,0.5,0.2), fn =minusQloglik_Vasicek,
             logret=logret, logprice_old=logprice_old, 
             deltati=deltati, lower=c(-Inf,-Inf,0,0), method ="L-BFGS-B")

res$par
res$value
res$convergence


# EXERCISE 4 
# POINT 1
# Estimate the volatility (on yearly basis) using the MLE approach.
library(quantmod)
getSymbols("^GSPC",  from="2019-01-02", to="2020-12-31")
tail(GSPC)
Price <- na.omit(GSPC[,4]) 
daysInOneYear<- 365 
timeOfObervation <-as.numeric(index(Price))
deltati<-diff(timeOfObervation)/daysInOneYear
Xt_i <- diff(log(as.numeric(Price)))

minuslogLik<- function(par,Xt_i,deltati){
  mu<-par[1]
  sig<-par[2]
  vecMean <-  (mu-0.5*sig^2)*deltati
  vecsd <-sig*sqrt(deltati)
  -sum(dnorm(Xt_i,mean=vecMean,sd=vecsd,log=TRUE))
}

minuslogLik(par=c(0.5,0.2),Xt_i=Xt_i,deltati=deltati)

res <- optim(par=c(0.5,0.2),fn = minuslogLik, lower=c(-Inf,0),method ="L-BFGS-B", 
             Xt_i=Xt_i,deltati=deltati)
res$par
-1*res$value
res$convergence
volatility <- res$par[2]
volatility

# POINT 2

#Inputs
#t0 (today) is December 30th 2019
#Maturity of the option is February 27th 2020
sigma=volatility
S0=as.numeric(Price[504,])
k=S0
r=0.015
TimeToMat=60 
daysInOneYear <- 365
deltati = TimeToMat/daysInOneYear 
d1 <- (log(S0/k)+(r+0.5*sigma^2)*(deltati))/(sigma*sqrt(deltati))
d2 <- d1-sigma*sqrt(deltati)

BSput=k*exp(-r*(deltati))*pnorm(-d2)- S0*pnorm(-d1)
 

