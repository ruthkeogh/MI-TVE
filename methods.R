
#---------------------
#packages
#---------------------

library(survival)
library(timereg)
library(mice)
library(smcfcs)

#-----------------------------------------
#-----------------------------------------
#complete-data analysis
#-----------------------------------------
#-----------------------------------------

#---
#split survival data on event times

complete.data.split=survSplit(Surv(t,d)~.,complete.data,cut=complete.data$t[complete.data$d==1])

#---
#create restricted cubic spline functions based on 5 knots (note: this could also be achieved using the rms package)

v=function(t,km1,kmax,k){(pmax(t-k,0))^3-((kmax-k)/(kmax-km1))*(pmax(t-km1,0))^3-((k-km1)/(kmax-km1))*(pmax(t-kmax,0))^3} #this is the spline basis function

knots.5=quantile(complete.data$t[complete.data$d==1],probs=c(0.05,0.25,0.5,0.75,0.95)) #knot positions

complete.data.split$t.sp5.1=sapply(complete.data.split$t,function(t)v(t,knots.5[4],knots.5[5],knots.5[1]))
complete.data.split$t.sp5.2=sapply(complete.data.split$t,function(t)v(t,knots.5[4],knots.5[5],knots.5[2]))
complete.data.split$t.sp5.3=sapply(complete.data.split$t,function(t)v(t,knots.5[4],knots.5[5],knots.5[3]))

complete.data.split$x1t=complete.data.split$x1*complete.data.split$t
complete.data.split$x2t=complete.data.split$x2*complete.data.split$t

complete.data.split$x1t.sp1=complete.data.split$x1*complete.data.split$t.sp5.1
complete.data.split$x2t.sp1=complete.data.split$x2*complete.data.split$t.sp5.1

complete.data.split$x1t.sp2=complete.data.split$x1*complete.data.split$t.sp5.2
complete.data.split$x2t.sp2=complete.data.split$x2*complete.data.split$t.sp5.2

complete.data.split$x1t.sp3=complete.data.split$x1*complete.data.split$t.sp5.3
complete.data.split$x2t.sp3=complete.data.split$x2*complete.data.split$t.sp5.3

#---
#fit the Cox model with TVEs to the complete data

comp.data.mod=coxph(Surv(tstart,t,d)~x1+x1t+x1t.sp1+x1t.sp2+x1t.sp3+
                      x2+x2t+x2t.sp1+x2t.sp2+x2t.sp3,
                    data=complete.data.split)

#---
#tests of the proportional hazards assumption for x1 and x2, using joint wald tests

wald.test(coef=comp.data.mod$coef[2:5],Sigma=comp.data.mod$var[2:5,2:5],print.coef=F)$p.value
wald.test(coef=comp.data.mod$coef[7:10],Sigma=comp.data.mod$var[7:10,7:10],print.coef=F)$p.value

#-----------------------------------------
#-----------------------------------------
#complete-case analysis
#-----------------------------------------
#-----------------------------------------

#---
#split survival data on event times

data.with.missingness.split=survSplit(Surv(t,d)~.,data.with.missingness,cut=data.with.missingness$t[data.with.missingness$d==1])

#---
#create restricted cubic spline functions based on 5 knots (note: this could also be achieved using the rms package)

v=function(t,km1,kmax,k){(pmax(t-k,0))^3-((kmax-k)/(kmax-km1))*(pmax(t-km1,0))^3-((k-km1)/(kmax-km1))*(pmax(t-kmax,0))^3} #this is the spline basis function

knots.5=quantile(data.with.missingness$t[data.with.missingness$d==1],probs=c(0.05,0.25,0.5,0.75,0.95)) #knot positions

data.with.missingness.split$t.sp5.1=sapply(data.with.missingness.split$t,function(t)v(t,knots.5[4],knots.5[5],knots.5[1]))
data.with.missingness.split$t.sp5.2=sapply(data.with.missingness.split$t,function(t)v(t,knots.5[4],knots.5[5],knots.5[2]))
data.with.missingness.split$t.sp5.3=sapply(data.with.missingness.split$t,function(t)v(t,knots.5[4],knots.5[5],knots.5[3]))

data.with.missingness.split$x1t=data.with.missingness.split$x1*data.with.missingness.split$t
data.with.missingness.split$x2t=data.with.missingness.split$x2*data.with.missingness.split$t

data.with.missingness.split$x1t.sp1=data.with.missingness.split$x1*data.with.missingness.split$t.sp5.1
data.with.missingness.split$x2t.sp1=data.with.missingness.split$x2*data.with.missingness.split$t.sp5.1

data.with.missingness.split$x1t.sp2=data.with.missingness.split$x1*data.with.missingness.split$t.sp5.2
data.with.missingness.split$x2t.sp2=data.with.missingness.split$x2*data.with.missingness.split$t.sp5.2

data.with.missingness.split$x1t.sp3=data.with.missingness.split$x1*data.with.missingness.split$t.sp5.3
data.with.missingness.split$x2t.sp3=data.with.missingness.split$x2*data.with.missingness.split$t.sp5.3

#---
#fit the Cox model with TVEs to the complete cases

comp.case.mod=coxph(Surv(tstart,t,d)~x1+x1t+x1t.sp1+x1t.sp2+x1t.sp3+
                      x2+x2t+x2t.sp1+x2t.sp2+x2t.sp3,
                    data=data.with.missingness.split)

#---
#tests of the proportional hazards assumption for x1 and x2, using joint wald tests

wald.test(coef=comp.case.mod$coef[2:5],Sigma=comp.case.mod$var[2:5,2:5],print.coef=F)$p.value
wald.test(coef=comp.case.mod$coef[7:10],Sigma=comp.case.mod$var[7:10,7:10],print.coef=F)$p.value

#-----------------------------------------
#-----------------------------------------
#MI-Approx
#-----------------------------------------
#-----------------------------------------

#---
#obtain Nelson-Aalen estimate of the cumulative hazard

data.with.missingness$H0=nelsonaalen(data.with.missingness,t,d)

#---
#perform the imputation

predmat=matrix(0,nrow=dim(data.with.missingness)[2],ncol=dim(data.with.missingness)[2])
rownames(predmat)=names(data.with.missingness)
colnames(predmat)=names(data.with.missingness)
predmat["x1",c("d","H0","x2")]=1
predmat["x2",c("d","H0","x1")]=1
methods=matrix("",nrow=dim(data.with.missingness)[2],ncol=1)
rownames(methods)=names(data.with.missingness)
methods[c("x1","x2"),]="logreg" #for continuous covariates "logreg" is replaced by "norm"

mi_approx=mice(data.with.missingness,m=10,method=as.vector(methods),predictorMatrix = predmat)

#---
#perform Cox regression with TVEs for x and x2, using the imputed data sets

mi_approx.coef=matrix(nrow=10,ncol=10)
mi_approx.varcov=array(dim=c(10,10,10))
for(imp in 1:10){
  temp=survSplit(Surv(t,d)~.,complete(mi_approx,imp),cut=data.with.missingness$t[data.with.missingness$d==1])
  temp$t.sp5.1=data.with.missingness.split$t.sp5.1
  temp$t.sp5.2=data.with.missingness.split$t.sp5.2
  temp$t.sp5.3=data.with.missingness.split$t.sp5.3
  mi_approx.mod=coxph(Surv(tstart,t,d)~x1+x1:t+x1:t.sp5.1+x1:t.sp5.2+x1:t.sp5.3+
                        x2+x2:t+x2:t.sp5.1+x2:t.sp5.2+x2:t.sp5.3,data=temp)
  mi_approx.coef[imp,]=mi_approx.mod$coef
  mi_approx.varcov[imp,,]=mi_approx.mod$var}

#---
#apply Rubin's rules

mi_approx.pool=MIcombine(list(mi_approx.coef[1,],mi_approx.coef[2,],mi_approx.coef[3,],mi_approx.coef[4,],mi_approx.coef[5,],
                                  mi_approx.coef[6,],mi_approx.coef[7,],mi_approx.coef[8,],mi_approx.coef[9,],mi_approx.coef[10,]),
                          list(mi_approx.varcov[1,,],mi_approx.varcov[2,,],mi_approx.varcov[3,,],mi_approx.varcov[4,,],mi_approx.varcov[5,,],
                               mi_approx.varcov[6,,],mi_approx.varcov[7,,],mi_approx.varcov[8,,],mi_approx.varcov[9,,],mi_approx.varcov[10,,]))

#---
#tests of the proportional hazards assumption for x1 and x2, using joint wald tests

wald.test(coef=mi_approx.pool$coef[2:5],Sigma=mi_approx.pool$var[2:5,2:5],print.coef = F)$p.value
wald.test(coef=mi_approx.pool$coef[7:10],Sigma=mi_approx.pool$var[7:10,7:10],print.coef = F)$p.value

#-----------------------------------------
#-----------------------------------------
#MI-TVE-Approx
#-----------------------------------------
#-----------------------------------------

#---
#obtain Nelson-Aalen estimate of the cumulative hazard

data.with.missingness$H0=nelsonaalen(data.with.missingness,t,d)

#---
#generate interactions between D (event indicator) and the TVE function (a restricted cubic spline with 5 knots)

data.with.missingness$t.sp5.1=sapply(data.with.missingness$t,function(t)v(t,knots.5[4],knots.5[5],knots.5[1]))
data.with.missingness$t.sp5.2=sapply(data.with.missingness$t,function(t)v(t,knots.5[4],knots.5[5],knots.5[2]))
data.with.missingness$t.sp5.3=sapply(data.with.missingness$t,function(t)v(t,knots.5[4],knots.5[5],knots.5[3]))

data.with.missingness$dt=data.with.missingness$d*data.with.missingness$t
data.with.missingness$dt.sp5.1=data.with.missingness$d*data.with.missingness$t.sp5.1
data.with.missingness$dt.sp5.2=data.with.missingness$d*data.with.missingness$t.sp5.2
data.with.missingness$dt.sp5.3=data.with.missingness$d*data.with.missingness$t.sp5.3

#---
#perform the imputation

predmat=matrix(0,nrow=dim(data.with.missingness)[2],ncol=dim(data.with.missingness)[2])
rownames(predmat)=names(data.with.missingness)
colnames(predmat)=names(data.with.missingness)
predmat["x1",c("d","H0","x2","t","t.sp5.1","t.sp5.2","t.sp5.3","dt","dt.sp5.1","dt.sp5.2","dt.sp5.3")]=1
predmat["x2",c("d","H0","x1","t","t.sp5.1","t.sp5.2","t.sp5.3","dt","dt.sp5.1","dt.sp5.2","dt.sp5.3")]=1
methods=matrix("",nrow=dim(data.with.missingness)[2],ncol=1)
rownames(methods)=names(data.with.missingness)
methods[c("x1","x2"),]="logreg" #for continuous covariates "logreg" is replaced by "norm"

mi_tve_approx=mice(data.with.missingness,m=10,method=as.vector(methods),predictorMatrix = predmat)

#---
#perform Cox regression with TVEs for x and x2, using the imputed data sets

mi_tve_approx.coef=matrix(nrow=10,ncol=10)
mi_tve_approx.varcov=array(dim=c(10,10,10))
for(imp in 1:10){
  temp=survSplit(Surv(t,d)~.,complete(mi_tve_approx,imp),cut=data.with.missingness$t[data.with.missingness$d==1])
  temp$t.sp5.1=data.with.missingness.split$t.sp5.1
  temp$t.sp5.2=data.with.missingness.split$t.sp5.2
  temp$t.sp5.3=data.with.missingness.split$t.sp5.3
  mi_tve_approx.mod=coxph(Surv(tstart,t,d)~x1+x1:t+x1:t.sp5.1+x1:t.sp5.2+x1:t.sp5.3+
                        x2+x2:t+x2:t.sp5.1+x2:t.sp5.2+x2:t.sp5.3,data=temp)
  mi_tve_approx.coef[imp,]=mi_tve_approx.mod$coef
  mi_tve_approx.varcov[imp,,]=mi_tve_approx.mod$var}

#---
#apply Rubin's rules

mi_tve_approx.pool=MIcombine(list(mi_tve_approx.coef[1,],mi_tve_approx.coef[2,],mi_tve_approx.coef[3,],mi_tve_approx.coef[4,],mi_tve_approx.coef[5,],
                              mi_tve_approx.coef[6,],mi_tve_approx.coef[7,],mi_tve_approx.coef[8,],mi_tve_approx.coef[9,],mi_tve_approx.coef[10,]),
                         list(mi_tve_approx.varcov[1,,],mi_tve_approx.varcov[2,,],mi_tve_approx.varcov[3,,],mi_tve_approx.varcov[4,,],mi_tve_approx.varcov[5,,],
                              mi_tve_approx.varcov[6,,],mi_tve_approx.varcov[7,,],mi_tve_approx.varcov[8,,],mi_tve_approx.varcov[9,,],mi_tve_approx.varcov[10,,]))

#---
#tests of the proportional hazards assumption for x1 and x2, using joint wald tests

wald.test(coef=mi_tve_approx.pool$coef[2:5],Sigma=mi_tve_approx.pool$var[2:5,2:5],print.coef = F)$p.value
wald.test(coef=mi_tve_approx.pool$coef[7:10],Sigma=mi_tve_approx.pool$var[7:10,7:10],print.coef = F)$p.value

#-----------------------------------------
#-----------------------------------------
#MI-SMC
#-----------------------------------------
#-----------------------------------------

#---
#perform the imputation

method=as.matrix(rep("",dim(data.with.missingness)[2]))
rownames(method)=names(data.with.missingness)
method[c("x1","x2"),]="logreg" #for continuous covariates "logreg" is replaced by "norm"
pred.mat=matrix(0,nrow=length(method),ncol=length(method))
colnames(pred.mat)=names(data.with.missingness)
rownames(pred.mat)=names(data.with.missingness)
pred.mat["x1","x2"]=1
pred.mat["x2","x1"]=1

mi_smc=smcfcs(originaldata=data.with.missingness,smtype="coxph",smformula="Surv(t, d) ~ x1+x2",
                      method=method,predictorMatrix = pred.mat,m=10,numit=10,rjlimit=5000,noisy=F)

#---
#perform Cox regression with TVEs for x and x2, using the imputed data sets

mi_smc.coef=matrix(nrow=10,ncol=10)
mi_smc.varcov=array(dim=c(10,10,10))
for(imp in 1:10){
  temp=survSplit(Surv(t,d)~.,mi_smc$impDatasets[[imp]],cut=data.with.missingness$t[data.with.missingness$d==1])
  temp$t.sp5.1=data.with.missingness.split$t.sp5.1
  temp$t.sp5.2=data.with.missingness.split$t.sp5.2
  temp$t.sp5.3=data.with.missingness.split$t.sp5.3
  mi_smc.mod=coxph(Surv(tstart,t,d)~x1+x1:t+x1:t.sp5.1+x1:t.sp5.2+x1:t.sp5.3+
                            x2+x2:t+x2:t.sp5.1+x2:t.sp5.2+x2:t.sp5.3,data=temp)
  mi_smc.coef[imp,]=mi_smc.mod$coef
  mi_smc.varcov[imp,,]=mi_smc.mod$var}

#---
#apply Rubin's rules

mi_smc.pool=MIcombine(list(mi_smc.coef[1,],mi_smc.coef[2,],mi_smc.coef[3,],mi_smc.coef[4,],mi_smc.coef[5,],
                                  mi_smc.coef[6,],mi_smc.coef[7,],mi_smc.coef[8,],mi_smc.coef[9,],mi_smc.coef[10,]),
                             list(mi_smc.varcov[1,,],mi_smc.varcov[2,,],mi_smc.varcov[3,,],mi_smc.varcov[4,,],mi_smc.varcov[5,,],
                                  mi_smc.varcov[6,,],mi_smc.varcov[7,,],mi_smc.varcov[8,,],mi_smc.varcov[9,,],mi_smc.varcov[10,,]))

#---
#tests of the proportional hazards assumption for x1 and x2, using joint wald tests

wald.test(coef=mi_smc.pool$coef[2:5],Sigma=mi_smc.pool$var[2:5,2:5],print.coef = F)$p.value
wald.test(coef=mi_smc.pool$coef[7:10],Sigma=mi_smc.pool$var[7:10,7:10],print.coef = F)$p.value

#-----------------------------------------
#-----------------------------------------
#MI-TVE-SMC
#-----------------------------------------
#-----------------------------------------

#---
#perform the imputation

method=as.matrix(rep("",dim(data.with.missingness)[2]))
rownames(method)=names(data.with.missingness)
method[c("x1","x2"),]="logreg" #for continuous covariates "logreg" is replaced by "norm"
pred.mat=matrix(0,nrow=length(method),ncol=length(method))
colnames(pred.mat)=names(data.with.missingness)
rownames(pred.mat)=names(data.with.missingness)
pred.mat["x1","x2"]=1
pred.mat["x2","x1"]=1

source("smc_tve") #you will need to use the appropriate file path here
mi_tve_smc=smcfcs.tve(originaldata=data.with.missingness,smformula.timefixed="Surv(t, d) ~ x1+x2",
              method=method,tvar=c("x1","x2"), tvarfunc=NULL,
             spline=T,knots=NULL,nknots=5,predictorMatrix = pred.mat,m=10,numit=10,rjlimit=5000,noisy=F,quick=F)

#---
#perform Cox regression with TVEs for x and x2, using the imputed data sets

mi_tve_smc.coef=matrix(nrow=10,ncol=10)
mi_tve_smc.varcov=array(dim=c(10,10,10))
for(imp in 1:10){
  temp=survSplit(Surv(t,d)~.,mi_tve_smc$impDatasets[[imp]],cut=data.with.missingness$t[data.with.missingness$d==1])
  temp$t.sp5.1=data.with.missingness.split$t.sp5.1
  temp$t.sp5.2=data.with.missingness.split$t.sp5.2
  temp$t.sp5.3=data.with.missingness.split$t.sp5.3
  mi_tve_smc.mod=coxph(Surv(tstart,t,d)~x1+x1:t+x1:t.sp5.1+x1:t.sp5.2+x1:t.sp5.3+
                     x2+x2:t+x2:t.sp5.1+x2:t.sp5.2+x2:t.sp5.3,data=temp)
  mi_tve_smc.coef[imp,]=mi_tve_smc.mod$coef
  mi_tve_smc.varcov[imp,,]=mi_tve_smc.mod$var}

#---
#apply Rubin's rules

mi_tve_smc.pool=MIcombine(list(mi_tve_smc.coef[1,],mi_tve_smc.coef[2,],mi_tve_smc.coef[3,],mi_tve_smc.coef[4,],mi_tve_smc.coef[5,],
                           mi_tve_smc.coef[6,],mi_tve_smc.coef[7,],mi_tve_smc.coef[8,],mi_tve_smc.coef[9,],mi_tve_smc.coef[10,]),
                      list(mi_tve_smc.varcov[1,,],mi_tve_smc.varcov[2,,],mi_tve_smc.varcov[3,,],mi_tve_smc.varcov[4,,],mi_tve_smc.varcov[5,,],
                           mi_tve_smc.varcov[6,,],mi_tve_smc.varcov[7,,],mi_tve_smc.varcov[8,,],mi_tve_smc.varcov[9,,],mi_tve_smc.varcov[10,,]))

#---
#tests of the proportional hazards assumption for x1 and x2, using joint wald tests

wald.test(coef=mi_tve_smc.pool$coef[2:5],Sigma=mi_tve_smc.pool$var[2:5,2:5],print.coef = F)$p.value
wald.test(coef=mi_tve_smc.pool$coef[7:10],Sigma=mi_tve_smc.pool$var[7:10,7:10],print.coef = F)$p.value



