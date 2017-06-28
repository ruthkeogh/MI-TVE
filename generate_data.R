
#---------------------
#set seed
#---------------------

set.seed(1911)

#---------------------
#sample size
#---------------------

n=2000

#---------------------
#Generate covariates - illustrated here for binary covariates
#---------------------

x1=rbinom(n,1,0.2)

expit=function(x){exp(x)/(1+exp(x))} #expit funcion, used below

x2=rbinom(n,1,expit(x1))

#for continuous covariates, this can be replaced with the following
#library(MASS)
#x=mvrnorm(n,mu=c(0,0),Sigma=matrix(c(1,0.5,0.5,1),nrow=2,ncol=2))
#x1=x[,1]
#x2=x[,2]

#---------------------
#Generate event times
#---------------------

delta=0.01
t.grid=seq(0,10-delta,delta) #fine grid of times used to generate event times in a piecewise fashion
lambda=0.01 #baseline hazard
lambda.c=0.076 #hazard for censoring

t.event.pw=matrix(nrow=n,ncol=length(t.grid))

for(i in 1:length(t.grid)){
  u=runif(n,0,1)
  t.event.pw.a=t.grid[i]-log(u)/(lambda*exp((0.32+1.42*exp(-t.grid[i])-0.02*(t.grid[i]^0.7))*x1+0.5*x2))
  t.event.pw[,i]=ifelse(t.event.pw.a>=t.grid[i] & t.event.pw.a<(t.grid[i]+delta),t.event.pw.a,NA)
}

t.event=apply(t.event.pw,1,function(x) min(x,na.rm=T))
t.event=ifelse(t.event==Inf,100,t.event)

#---------------------
#Generate censoring times
#---------------------

u.cens=runif(n,0,1)
t.cens=-log(u.cens)/lambda.c

#---------------------
#Generate event/censoring indicator
#1=event,2=random censored,3=administrative censoring
#---------------------

d.a=ifelse(t.event<t.cens,1,2)
d.a=ifelse(pmin(t.event,t.cens)>10,3,d.a)

table(d.a)/n

#---------------------
#Generate time and indicator used in analyses
#---------------------

d=ifelse(d.a==1,1,0)
t=pmin(t.event,t.cens)
t=ifelse(t>10,10,t)

#---------------------
#Generate missing data in x1 and x2
#---------------------

#missing data parameters, used below

x1.m.a=0.4
x1.m.b=0.5
x2.m.a=0.4
x2.m.b=0.5
x1x2.m=0.3

#divide the data into 3 approximately equal sections

u.miss=runif(n,0,1)
section=cut(u.miss,c(0,0.33,0.66,1),labels=F)

#in section 1 X2 is fully observed and x1 is missing depending on x2

x1.m=rbinom(n,1,expit(x1.m.a+x1.m.b*x2))

#in section 2 x1 is fully observed and x2 is missing depending on x1

x2.m=rbinom(n,1,expit(x2.m.a+x2.m.b*x1))

#in section 3 both x1 and x2 are missing completely at random

x1x2.m=rbinom(n,1,x1x2.m)

#generate missing data in x1 and x2

x1.miss=ifelse(section==1 & x1.m==1,NA,x1)
x1.miss=ifelse(section==3 & x1x2.m==1,NA,x1.miss)

x2.miss=ifelse(section==2 & x2.m==1,NA,x2)
x2.miss=ifelse(section==3 & x1x2.m==1,NA,x2.miss)

#---------------------
#create data frames
#---------------------

complete.data=data.frame(id=1:n,t,d,x1,x2)

data.with.missingness=data.frame(id=1:n,t,d,x1=x1.miss,x2=x2.miss)
