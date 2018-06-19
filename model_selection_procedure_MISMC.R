

for(var in vars.remaining){
  eval(parse(text=paste0("mi.",var,".k5.coef=matrix(nrow=nimp,ncol=4)")))
  eval(parse(text=paste0("mi.",var,".k4.coef=matrix(nrow=nimp,ncol=3)")))
  eval(parse(text=paste0("mi.",var,".k3.coef=matrix(nrow=nimp,ncol=2)")))
  eval(parse(text=paste0("mi.",var,".lin.coef=matrix(nrow=nimp,ncol=1)")))
  
  eval(parse(text=paste0("mi.",var,".k5.varcov=array(dim=c(nimp,4,4))")))
  eval(parse(text=paste0("mi.",var,".k4.varcov=array(dim=c(nimp,3,3))")))
  eval(parse(text=paste0("mi.",var,".k3.varcov=array(dim=c(nimp,2,2))")))
  eval(parse(text=paste0("mi.",var,".lin.varcov=array(dim=c(nimp,1,1))")))
}

for(imp in 1:nimp){
  
  print(imp)
  temp=survSplit(Surv(t,d)~.,mi.dat$impDatasets[[imp]],cut=cut.points)
  temp=temp[temp$t>temp$tstart,]
  temp$t.k5.1=mydata.split$t.k5.1
  temp$t.k5.2=mydata.split$t.k5.2
  temp$t.k5.3=mydata.split$t.k5.3
  temp$t.k4.1=mydata.split$t.k4.1
  temp$t.k4.2=mydata.split$t.k4.2
  temp$t.k3=mydata.split$t.k3

#-------------
#age

if("age"%in%vars.remaining){
  formula.tve.sp=as.formula(paste(formula.base,"+age:t+age:t.k5.1+age:t.k5.2+age:t.k5.3"))
  mi.age.k5=coxph(formula.tve.sp,data=temp)
  
  formula.tve.sp=as.formula(paste(formula.base,"+age:t+age:t.k4.1+age:t.k4.2"))
  mi.age.k4=coxph(formula.tve.sp,data=temp)

  formula.tve.sp=as.formula(paste(formula.base,"+age:t+age:t.k3"))
  mi.age.k3=coxph(formula.tve.sp,data=temp)

  formula.tve.sp=as.formula(paste(formula.base,"+age:t"))
  mi.age.lin=coxph(formula.tve.sp,data=temp)
  
  mi.age.k5.coef[imp,]=mi.age.k5$coef[c("age:t","age:t.k5.1","age:t.k5.2","age:t.k5.3")]
  mi.age.k4.coef[imp,]=mi.age.k4$coef[c("age:t","age:t.k4.1","age:t.k4.2")]
  mi.age.k3.coef[imp,]=mi.age.k3$coef[c("age:t","age:t.k3")]
  mi.age.lin.coef[imp,]=mi.age.lin$coef[c("age:t")]
  
  mi.age.k5.varcov[imp,,]=vcov(mi.age.k5)[c("age:t","age:t.k5.1","age:t.k5.2","age:t.k5.3"),c("age:t","age:t.k5.1","age:t.k5.2","age:t.k5.3")]
  mi.age.k4.varcov[imp,,]=vcov(mi.age.k4)[c("age:t","age:t.k4.1","age:t.k4.2"),c("age:t","age:t.k4.1","age:t.k4.2")]
  mi.age.k3.varcov[imp,,]=vcov(mi.age.k3)[c("age:t","age:t.k3"),c("age:t","age:t.k3")]
  mi.age.lin.varcov[imp,,]=vcov(mi.age.lin)[c("age:t"),c("age:t")]

}

#-------------
#sized1

if("sized1"%in%vars.remaining){
  formula.tve.sp=as.formula(paste(formula.base,"+sized1:t+sized1:t.k5.1+sized1:t.k5.2+sized1:t.k5.3"))
  mi.sized1.k5=coxph(formula.tve.sp,data=temp)
  
  formula.tve.sp=as.formula(paste(formula.base,"+sized1:t+sized1:t.k4.1+sized1:t.k4.2"))
  mi.sized1.k4=coxph(formula.tve.sp,data=temp)
  
  formula.tve.sp=as.formula(paste(formula.base,"+sized1:t+sized1:t.k3"))
  mi.sized1.k3=coxph(formula.tve.sp,data=temp)
  
  formula.tve.sp=as.formula(paste(formula.base,"+sized1:t"))
  mi.sized1.lin=coxph(formula.tve.sp,data=temp)
  
  mi.sized1.k5.coef[imp,]=mi.sized1.k5$coef[c("sized1:t","sized1:t.k5.1","sized1:t.k5.2","sized1:t.k5.3")]
  mi.sized1.k4.coef[imp,]=mi.sized1.k4$coef[c("sized1:t","sized1:t.k4.1","sized1:t.k4.2")]
  mi.sized1.k3.coef[imp,]=mi.sized1.k3$coef[c("sized1:t","sized1:t.k3")]
  mi.sized1.lin.coef[imp,]=mi.sized1.lin$coef[c("sized1:t")]
  
  mi.sized1.k5.varcov[imp,,]=vcov(mi.sized1.k5)[c("sized1:t","sized1:t.k5.1","sized1:t.k5.2","sized1:t.k5.3"),c("sized1:t","sized1:t.k5.1","sized1:t.k5.2","sized1:t.k5.3")]
  mi.sized1.k4.varcov[imp,,]=vcov(mi.sized1.k4)[c("sized1:t","sized1:t.k4.1","sized1:t.k4.2"),c("sized1:t","sized1:t.k4.1","sized1:t.k4.2")]
  mi.sized1.k3.varcov[imp,,]=vcov(mi.sized1.k3)[c("sized1:t","sized1:t.k3"),c("sized1:t","sized1:t.k3")]
  mi.sized1.lin.varcov[imp,,]=vcov(mi.sized1.lin)[c("sized1:t"),c("sized1:t")]
  
}
#-------------
#sized2

if("sized2"%in%vars.remaining){
  formula.tve.sp=as.formula(paste(formula.base,"+sized2:t+sized2:t.k5.1+sized2:t.k5.2+sized2:t.k5.3"))
  mi.sized2.k5=coxph(formula.tve.sp,data=temp)
  
  formula.tve.sp=as.formula(paste(formula.base,"+sized2:t+sized2:t.k4.1+sized2:t.k4.2"))
  mi.sized2.k4=coxph(formula.tve.sp,data=temp)
  
  formula.tve.sp=as.formula(paste(formula.base,"+sized2:t+sized2:t.k3"))
  mi.sized2.k3=coxph(formula.tve.sp,data=temp)
  
  formula.tve.sp=as.formula(paste(formula.base,"+sized2:t"))
  mi.sized2.lin=coxph(formula.tve.sp,data=temp)
  
  mi.sized2.k5.coef[imp,]=mi.sized2.k5$coef[c("sized2:t","sized2:t.k5.1","sized2:t.k5.2","sized2:t.k5.3")]
  mi.sized2.k4.coef[imp,]=mi.sized2.k4$coef[c("sized2:t","sized2:t.k4.1","sized2:t.k4.2")]
  mi.sized2.k3.coef[imp,]=mi.sized2.k3$coef[c("sized2:t","sized2:t.k3")]
  mi.sized2.lin.coef[imp,]=mi.sized2.lin$coef[c("sized2:t")]
  
  mi.sized2.k5.varcov[imp,,]=vcov(mi.sized2.k5)[c("sized2:t","sized2:t.k5.1","sized2:t.k5.2","sized2:t.k5.3"),c("sized2:t","sized2:t.k5.1","sized2:t.k5.2","sized2:t.k5.3")]
  mi.sized2.k4.varcov[imp,,]=vcov(mi.sized2.k4)[c("sized2:t","sized2:t.k4.1","sized2:t.k4.2"),c("sized2:t","sized2:t.k4.1","sized2:t.k4.2")]
  mi.sized2.k3.varcov[imp,,]=vcov(mi.sized2.k3)[c("sized2:t","sized2:t.k3"),c("sized2:t","sized2:t.k3")]
  mi.sized2.lin.varcov[imp,,]=vcov(mi.sized2.lin)[c("sized2:t"),c("sized2:t")]
  
}

#-------------
#grade

if("grade"%in%vars.remaining){
  formula.tve.sp=as.formula(paste(formula.base,"+grade:t+grade:t.k5.1+grade:t.k5.2+grade:t.k5.3"))
  mi.grade.k5=coxph(formula.tve.sp,data=temp)
  
  formula.tve.sp=as.formula(paste(formula.base,"+grade:t+grade:t.k4.1+grade:t.k4.2"))
  mi.grade.k4=coxph(formula.tve.sp,data=temp)
  
  formula.tve.sp=as.formula(paste(formula.base,"+grade:t+grade:t.k3"))
  mi.grade.k3=coxph(formula.tve.sp,data=temp)
  
  formula.tve.sp=as.formula(paste(formula.base,"+grade:t"))
  mi.grade.lin=coxph(formula.tve.sp,data=temp)
  
  mi.grade.k5.coef[imp,]=mi.grade.k5$coef[c("grade:t","grade:t.k5.1","grade:t.k5.2","grade:t.k5.3")]
  mi.grade.k4.coef[imp,]=mi.grade.k4$coef[c("grade:t","grade:t.k4.1","grade:t.k4.2")]
  mi.grade.k3.coef[imp,]=mi.grade.k3$coef[c("grade:t","grade:t.k3")]
  mi.grade.lin.coef[imp,]=mi.grade.lin$coef[c("grade:t")]
  
  mi.grade.k5.varcov[imp,,]=vcov(mi.grade.k5)[c("grade:t","grade:t.k5.1","grade:t.k5.2","grade:t.k5.3"),c("grade:t","grade:t.k5.1","grade:t.k5.2","grade:t.k5.3")]
  mi.grade.k4.varcov[imp,,]=vcov(mi.grade.k4)[c("grade:t","grade:t.k4.1","grade:t.k4.2"),c("grade:t","grade:t.k4.1","grade:t.k4.2")]
  mi.grade.k3.varcov[imp,,]=vcov(mi.grade.k3)[c("grade:t","grade:t.k3"),c("grade:t","grade:t.k3")]
  mi.grade.lin.varcov[imp,,]=vcov(mi.grade.lin)[c("grade:t"),c("grade:t")]
  
}
#-------------
#enodessq

if("enodessq"%in%vars.remaining){
  formula.tve.sp=as.formula(paste(formula.base,"+enodessq:t+enodessq:t.k5.1+enodessq:t.k5.2+enodessq:t.k5.3"))
  mi.enodessq.k5=coxph(formula.tve.sp,data=temp)
  
  formula.tve.sp=as.formula(paste(formula.base,"+enodessq:t+enodessq:t.k4.1+enodessq:t.k4.2"))
  mi.enodessq.k4=coxph(formula.tve.sp,data=temp)
  
  formula.tve.sp=as.formula(paste(formula.base,"+enodessq:t+enodessq:t.k3"))
  mi.enodessq.k3=coxph(formula.tve.sp,data=temp)
  
  formula.tve.sp=as.formula(paste(formula.base,"+enodessq:t"))
  mi.enodessq.lin=coxph(formula.tve.sp,data=temp)
  
  mi.enodessq.k5.coef[imp,]=mi.enodessq.k5$coef[c("enodessq:t","enodessq:t.k5.1","enodessq:t.k5.2","enodessq:t.k5.3")]
  mi.enodessq.k4.coef[imp,]=mi.enodessq.k4$coef[c("enodessq:t","enodessq:t.k4.1","enodessq:t.k4.2")]
  mi.enodessq.k3.coef[imp,]=mi.enodessq.k3$coef[c("enodessq:t","enodessq:t.k3")]
  mi.enodessq.lin.coef[imp,]=mi.enodessq.lin$coef[c("enodessq:t")]
  
  mi.enodessq.k5.varcov[imp,,]=vcov(mi.enodessq.k5)[c("enodessq:t","enodessq:t.k5.1","enodessq:t.k5.2","enodessq:t.k5.3"),c("enodessq:t","enodessq:t.k5.1","enodessq:t.k5.2","enodessq:t.k5.3")]
  mi.enodessq.k4.varcov[imp,,]=vcov(mi.enodessq.k4)[c("enodessq:t","enodessq:t.k4.1","enodessq:t.k4.2"),c("enodessq:t","enodessq:t.k4.1","enodessq:t.k4.2")]
  mi.enodessq.k3.varcov[imp,,]=vcov(mi.enodessq.k3)[c("enodessq:t","enodessq:t.k3"),c("enodessq:t","enodessq:t.k3")]
  mi.enodessq.lin.varcov[imp,,]=vcov(mi.enodessq.lin)[c("enodessq:t"),c("enodessq:t")]
  
}
#-------------
#hormon

if("hormon"%in%vars.remaining){
  formula.tve.sp=as.formula(paste(formula.base,"+hormon:t+hormon:t.k5.1+hormon:t.k5.2+hormon:t.k5.3"))
  mi.hormon.k5=coxph(formula.tve.sp,data=temp)
  
  formula.tve.sp=as.formula(paste(formula.base,"+hormon:t+hormon:t.k4.1+hormon:t.k4.2"))
  mi.hormon.k4=coxph(formula.tve.sp,data=temp)
  
  formula.tve.sp=as.formula(paste(formula.base,"+hormon:t+hormon:t.k3"))
  mi.hormon.k3=coxph(formula.tve.sp,data=temp)
  
  formula.tve.sp=as.formula(paste(formula.base,"+hormon:t"))
  mi.hormon.lin=coxph(formula.tve.sp,data=temp)
  
  mi.hormon.k5.coef[imp,]=mi.hormon.k5$coef[c("hormon:t","hormon:t.k5.1","hormon:t.k5.2","hormon:t.k5.3")]
  mi.hormon.k4.coef[imp,]=mi.hormon.k4$coef[c("hormon:t","hormon:t.k4.1","hormon:t.k4.2")]
  mi.hormon.k3.coef[imp,]=mi.hormon.k3$coef[c("hormon:t","hormon:t.k3")]
  mi.hormon.lin.coef[imp,]=mi.hormon.lin$coef[c("hormon:t")]
  
  mi.hormon.k5.varcov[imp,,]=vcov(mi.hormon.k5)[c("hormon:t","hormon:t.k5.1","hormon:t.k5.2","hormon:t.k5.3"),c("hormon:t","hormon:t.k5.1","hormon:t.k5.2","hormon:t.k5.3")]
  mi.hormon.k4.varcov[imp,,]=vcov(mi.hormon.k4)[c("hormon:t","hormon:t.k4.1","hormon:t.k4.2"),c("hormon:t","hormon:t.k4.1","hormon:t.k4.2")]
  mi.hormon.k3.varcov[imp,,]=vcov(mi.hormon.k3)[c("hormon:t","hormon:t.k3"),c("hormon:t","hormon:t.k3")]
  mi.hormon.lin.varcov[imp,,]=vcov(mi.hormon.lin)[c("hormon:t"),c("hormon:t")]
  
}

#-------------
#chemo

if("chemo"%in%vars.remaining){
  formula.tve.sp=as.formula(paste(formula.base,"+chemo:t+chemo:t.k5.1+chemo:t.k5.2+chemo:t.k5.3"))
  mi.chemo.k5=coxph(formula.tve.sp,data=temp)
  
  formula.tve.sp=as.formula(paste(formula.base,"+chemo:t+chemo:t.k4.1+chemo:t.k4.2"))
  mi.chemo.k4=coxph(formula.tve.sp,data=temp)
  
  formula.tve.sp=as.formula(paste(formula.base,"+chemo:t+chemo:t.k3"))
  mi.chemo.k3=coxph(formula.tve.sp,data=temp)
  
  formula.tve.sp=as.formula(paste(formula.base,"+chemo:t"))
  mi.chemo.lin=coxph(formula.tve.sp,data=temp)
  
  mi.chemo.k5.coef[imp,]=mi.chemo.k5$coef[c("chemo:t","chemo:t.k5.1","chemo:t.k5.2","chemo:t.k5.3")]
  mi.chemo.k4.coef[imp,]=mi.chemo.k4$coef[c("chemo:t","chemo:t.k4.1","chemo:t.k4.2")]
  mi.chemo.k3.coef[imp,]=mi.chemo.k3$coef[c("chemo:t","chemo:t.k3")]
  mi.chemo.lin.coef[imp,]=mi.chemo.lin$coef[c("chemo:t")]
  
  mi.chemo.k5.varcov[imp,,]=vcov(mi.chemo.k5)[c("chemo:t","chemo:t.k5.1","chemo:t.k5.2","chemo:t.k5.3"),c("chemo:t","chemo:t.k5.1","chemo:t.k5.2","chemo:t.k5.3")]
  mi.chemo.k4.varcov[imp,,]=vcov(mi.chemo.k4)[c("chemo:t","chemo:t.k4.1","chemo:t.k4.2"),c("chemo:t","chemo:t.k4.1","chemo:t.k4.2")]
  mi.chemo.k3.varcov[imp,,]=vcov(mi.chemo.k3)[c("chemo:t","chemo:t.k3"),c("chemo:t","chemo:t.k3")]
  mi.chemo.lin.varcov[imp,,]=vcov(mi.chemo.lin)[c("chemo:t"),c("chemo:t")]
  
}

#-------------
#logpgr

if("logpgr"%in%vars.remaining){
  formula.tve.sp=as.formula(paste(formula.base,"+logpgr:t+logpgr:t.k5.1+logpgr:t.k5.2+logpgr:t.k5.3"))
  mi.logpgr.k5=coxph(formula.tve.sp,data=temp)
  
  formula.tve.sp=as.formula(paste(formula.base,"+logpgr:t+logpgr:t.k4.1+logpgr:t.k4.2"))
  mi.logpgr.k4=coxph(formula.tve.sp,data=temp)
  
  formula.tve.sp=as.formula(paste(formula.base,"+logpgr:t+logpgr:t.k3"))
  mi.logpgr.k3=coxph(formula.tve.sp,data=temp)
  
  formula.tve.sp=as.formula(paste(formula.base,"+logpgr:t"))
  mi.logpgr.lin=coxph(formula.tve.sp,data=temp)
  
  mi.logpgr.k5.coef[imp,]=mi.logpgr.k5$coef[c("logpgr:t","logpgr:t.k5.1","logpgr:t.k5.2","logpgr:t.k5.3")]
  mi.logpgr.k4.coef[imp,]=mi.logpgr.k4$coef[c("logpgr:t","logpgr:t.k4.1","logpgr:t.k4.2")]
  mi.logpgr.k3.coef[imp,]=mi.logpgr.k3$coef[c("logpgr:t","logpgr:t.k3")]
  mi.logpgr.lin.coef[imp,]=mi.logpgr.lin$coef[c("logpgr:t")]
  
  mi.logpgr.k5.varcov[imp,,]=vcov(mi.logpgr.k5)[c("logpgr:t","logpgr:t.k5.1","logpgr:t.k5.2","logpgr:t.k5.3"),c("logpgr:t","logpgr:t.k5.1","logpgr:t.k5.2","logpgr:t.k5.3")]
  mi.logpgr.k4.varcov[imp,,]=vcov(mi.logpgr.k4)[c("logpgr:t","logpgr:t.k4.1","logpgr:t.k4.2"),c("logpgr:t","logpgr:t.k4.1","logpgr:t.k4.2")]
  mi.logpgr.k3.varcov[imp,,]=vcov(mi.logpgr.k3)[c("logpgr:t","logpgr:t.k3"),c("logpgr:t","logpgr:t.k3")]
  mi.logpgr.lin.varcov[imp,,]=vcov(mi.logpgr.lin)[c("logpgr:t"),c("logpgr:t")]
  
}
  
}

#----------------------------
#pooled results
#----------------------------

#-------------
#age

if("age"%in%vars.remaining){
  mi.age.k5.pool=MIcombine(eval(parse(text=paste0("list(",paste0("mi.age.k5.coef[",1:nimp,",]",collapse=","),")"))),
                               eval(parse(text=paste0("list(",paste0("mi.age.k5.varcov[",1:nimp,",,]",collapse=","),")"))))
  mi.age.k4.pool=MIcombine(eval(parse(text=paste0("list(",paste0("mi.age.k4.coef[",1:nimp,",]",collapse=","),")"))),
                               eval(parse(text=paste0("list(",paste0("mi.age.k4.varcov[",1:nimp,",,]",collapse=","),")"))))
  mi.age.k3.pool=MIcombine(eval(parse(text=paste0("list(",paste0("mi.age.k3.coef[",1:nimp,",]",collapse=","),")"))),
                               eval(parse(text=paste0("list(",paste0("mi.age.k3.varcov[",1:nimp,",,]",collapse=","),")"))))
  mi.age.lin.pool=MIcombine(eval(parse(text=paste0("list(",paste0("mi.age.lin.coef[",1:nimp,",]",collapse=","),")"))),
                                eval(parse(text=paste0("list(",paste0("mi.age.lin.varcov[",1:nimp,",,]",collapse=","),")"))))
  

  num.param=4
  pool.coef=mi.age.k5.pool$coefficients
  all.coef=mi.age.k5.coef
  all.var=mi.age.k5.varcov
  source("H:/MI_timedep_Cox/rotterdam_example/ruth_wald_test_mi_rotterdam_v2.R")
  wald.mi.age.k5=pval
  stat.mi.age.k5=d.x

  num.param=3
  pool.coef=mi.age.k4.pool$coefficients
  all.coef=mi.age.k4.coef
  all.var=mi.age.k4.varcov
  source("H:/MI_timedep_Cox/rotterdam_example/ruth_wald_test_mi_rotterdam_v2.R")
  wald.mi.age.k4=pval
  stat.mi.age.k4=d.x
  
  num.param=2
  pool.coef=mi.age.k3.pool$coefficients
  all.coef=mi.age.k3.coef
  all.var=mi.age.k3.varcov
  source("H:/MI_timedep_Cox/rotterdam_example/ruth_wald_test_mi_rotterdam_v2.R")
  wald.mi.age.k3=pval
  stat.mi.age.k3=d.x
  
  wald.mi.age.lin=wald.test(coef=mi.age.lin.pool$coef,Sigma=mi.age.lin.pool$var,print.coef = F)$p.value
  stat.mi.age.lin=wald.test(coef=mi.age.lin.pool$coef,Sigma=mi.age.lin.pool$var,print.coef = F)$statistic

  # waldtest.mi.age.k5=wald.test(coef=mi.age.k5.pool$coef,Sigma=mi.age.k5.pool$var,print.coef = F)
  # waldtest.mi.age.k4=wald.test(coef=mi.age.k4.pool$coef,Sigma=mi.age.k4.pool$var,print.coef = F)
  # waldtest.mi.age.k3=wald.test(coef=mi.age.k3.pool$coef,Sigma=mi.age.k3.pool$var,print.coef = F)
  # waldtest.mi.age.lin=wald.test(coef=mi.age.lin.pool$coef,Sigma=mi.age.lin.pool$var,print.coef = F)
  # 
  # wald.mi.age.k5=waldtest.mi.age.k5$p.value
  # wald.mi.age.k4=waldtest.mi.age.k4$p.value
  # wald.mi.age.k3=waldtest.mi.age.k3$p.value
  # wald.mi.age.lin=waldtest.mi.age.lin$p.value
  # 
  # stat.mi.age.k5=waldtest.mi.age.k5$statistic
  # stat.mi.age.k4=waldtest.mi.age.k4$statistic
  # stat.mi.age.k3=waldtest.mi.age.k3$statistic
  # stat.mi.age.lin=waldtest.mi.age.lin$statistic

}

#-------------
#sized1

if("sized1"%in%vars.remaining){
  mi.sized1.k5.pool=MIcombine(eval(parse(text=paste0("list(",paste0("mi.sized1.k5.coef[",1:nimp,",]",collapse=","),")"))),
                           eval(parse(text=paste0("list(",paste0("mi.sized1.k5.varcov[",1:nimp,",,]",collapse=","),")"))))
  mi.sized1.k4.pool=MIcombine(eval(parse(text=paste0("list(",paste0("mi.sized1.k4.coef[",1:nimp,",]",collapse=","),")"))),
                           eval(parse(text=paste0("list(",paste0("mi.sized1.k4.varcov[",1:nimp,",,]",collapse=","),")"))))
  mi.sized1.k3.pool=MIcombine(eval(parse(text=paste0("list(",paste0("mi.sized1.k3.coef[",1:nimp,",]",collapse=","),")"))),
                           eval(parse(text=paste0("list(",paste0("mi.sized1.k3.varcov[",1:nimp,",,]",collapse=","),")"))))
  mi.sized1.lin.pool=MIcombine(eval(parse(text=paste0("list(",paste0("mi.sized1.lin.coef[",1:nimp,",]",collapse=","),")"))),
                            eval(parse(text=paste0("list(",paste0("mi.sized1.lin.varcov[",1:nimp,",,]",collapse=","),")"))))
  
  num.param=4
  pool.coef=mi.sized1.k5.pool$coefficients
  all.coef=mi.sized1.k5.coef
  all.var=mi.sized1.k5.varcov
  source("H:/MI_timedep_Cox/rotterdam_example/ruth_wald_test_mi_rotterdam_v2.R")
  wald.mi.sized1.k5=pval
  stat.mi.sized1.k5=d.x
  
  num.param=3
  pool.coef=mi.sized1.k4.pool$coefficients
  all.coef=mi.sized1.k4.coef
  all.var=mi.sized1.k4.varcov
  source("H:/MI_timedep_Cox/rotterdam_example/ruth_wald_test_mi_rotterdam_v2.R")
  wald.mi.sized1.k4=pval
  stat.mi.sized1.k4=d.x
  
  num.param=2
  pool.coef=mi.sized1.k3.pool$coefficients
  all.coef=mi.sized1.k3.coef
  all.var=mi.sized1.k3.varcov
  source("H:/MI_timedep_Cox/rotterdam_example/ruth_wald_test_mi_rotterdam_v2.R")
  wald.mi.sized1.k3=pval
  stat.mi.sized1.k3=d.x
  
  wald.mi.sized1.lin=wald.test(coef=mi.sized1.lin.pool$coef,Sigma=mi.sized1.lin.pool$var,print.coef = F)$p.value
  stat.mi.sized1.lin=wald.test(coef=mi.sized1.lin.pool$coef,Sigma=mi.sized1.lin.pool$var,print.coef = F)$statistic
  
  # waldtest.mi.sized1.k5=wald.test(coef=mi.sized1.k5.pool$coef,Sigma=mi.sized1.k5.pool$var,print.coef = F)
  # waldtest.mi.sized1.k4=wald.test(coef=mi.sized1.k4.pool$coef,Sigma=mi.sized1.k4.pool$var,print.coef = F)
  # waldtest.mi.sized1.k3=wald.test(coef=mi.sized1.k3.pool$coef,Sigma=mi.sized1.k3.pool$var,print.coef = F)
  # waldtest.mi.sized1.lin=wald.test(coef=mi.sized1.lin.pool$coef,Sigma=mi.sized1.lin.pool$var,print.coef = F)
  # 
  # wald.mi.sized1.k5=waldtest.mi.sized1.k5$p.value
  # wald.mi.sized1.k4=waldtest.mi.sized1.k4$p.value
  # wald.mi.sized1.k3=waldtest.mi.sized1.k3$p.value
  # wald.mi.sized1.lin=waldtest.mi.sized1.lin$p.value
  # 
  # stat.mi.sized1.k5=waldtest.mi.sized1.k5$statistic
  # stat.mi.sized1.k4=waldtest.mi.sized1.k4$statistic
  # stat.mi.sized1.k3=waldtest.mi.sized1.k3$statistic
  # stat.mi.sized1.lin=waldtest.mi.sized1.lin$statistic
  
}


#-------------
#sized2

if("sized2"%in%vars.remaining){
  mi.sized2.k5.pool=MIcombine(eval(parse(text=paste0("list(",paste0("mi.sized2.k5.coef[",1:nimp,",]",collapse=","),")"))),
                           eval(parse(text=paste0("list(",paste0("mi.sized2.k5.varcov[",1:nimp,",,]",collapse=","),")"))))
  mi.sized2.k4.pool=MIcombine(eval(parse(text=paste0("list(",paste0("mi.sized2.k4.coef[",1:nimp,",]",collapse=","),")"))),
                           eval(parse(text=paste0("list(",paste0("mi.sized2.k4.varcov[",1:nimp,",,]",collapse=","),")"))))
  mi.sized2.k3.pool=MIcombine(eval(parse(text=paste0("list(",paste0("mi.sized2.k3.coef[",1:nimp,",]",collapse=","),")"))),
                           eval(parse(text=paste0("list(",paste0("mi.sized2.k3.varcov[",1:nimp,",,]",collapse=","),")"))))
  mi.sized2.lin.pool=MIcombine(eval(parse(text=paste0("list(",paste0("mi.sized2.lin.coef[",1:nimp,",]",collapse=","),")"))),
                            eval(parse(text=paste0("list(",paste0("mi.sized2.lin.varcov[",1:nimp,",,]",collapse=","),")"))))
  
  num.param=4
  pool.coef=mi.sized2.k5.pool$coefficients
  all.coef=mi.sized2.k5.coef
  all.var=mi.sized2.k5.varcov
  source("H:/MI_timedep_Cox/rotterdam_example/ruth_wald_test_mi_rotterdam_v2.R")
  wald.mi.sized2.k5=pval
  stat.mi.sized2.k5=d.x
  
  num.param=3
  pool.coef=mi.sized2.k4.pool$coefficients
  all.coef=mi.sized2.k4.coef
  all.var=mi.sized2.k4.varcov
  source("H:/MI_timedep_Cox/rotterdam_example/ruth_wald_test_mi_rotterdam_v2.R")
  wald.mi.sized2.k4=pval
  stat.mi.sized2.k4=d.x
  
  num.param=2
  pool.coef=mi.sized2.k3.pool$coefficients
  all.coef=mi.sized2.k3.coef
  all.var=mi.sized2.k3.varcov
  source("H:/MI_timedep_Cox/rotterdam_example/ruth_wald_test_mi_rotterdam_v2.R")
  wald.mi.sized2.k3=pval
  stat.mi.sized2.k3=d.x
  
  wald.mi.sized2.lin=wald.test(coef=mi.sized2.lin.pool$coef,Sigma=mi.sized2.lin.pool$var,print.coef = F)$p.value
  stat.mi.sized2.lin=wald.test(coef=mi.sized2.lin.pool$coef,Sigma=mi.sized2.lin.pool$var,print.coef = F)$statistic
  
  # waldtest.mi.sized2.k5=wald.test(coef=mi.sized2.k5.pool$coef,Sigma=mi.sized2.k5.pool$var,print.coef = F)
  # waldtest.mi.sized2.k4=wald.test(coef=mi.sized2.k4.pool$coef,Sigma=mi.sized2.k4.pool$var,print.coef = F)
  # waldtest.mi.sized2.k3=wald.test(coef=mi.sized2.k3.pool$coef,Sigma=mi.sized2.k3.pool$var,print.coef = F)
  # waldtest.mi.sized2.lin=wald.test(coef=mi.sized2.lin.pool$coef,Sigma=mi.sized2.lin.pool$var,print.coef = F)
  # 
  # wald.mi.sized2.k5=waldtest.mi.sized2.k5$p.value
  # wald.mi.sized2.k4=waldtest.mi.sized2.k4$p.value
  # wald.mi.sized2.k3=waldtest.mi.sized2.k3$p.value
  # wald.mi.sized2.lin=waldtest.mi.sized2.lin$p.value
  # 
  # stat.mi.sized2.k5=waldtest.mi.sized2.k5$statistic
  # stat.mi.sized2.k4=waldtest.mi.sized2.k4$statistic
  # stat.mi.sized2.k3=waldtest.mi.sized2.k3$statistic
  # stat.mi.sized2.lin=waldtest.mi.sized2.lin$statistic
  
}

#-------------
#grade

if("grade"%in%vars.remaining){
  mi.grade.k5.pool=MIcombine(eval(parse(text=paste0("list(",paste0("mi.grade.k5.coef[",1:nimp,",]",collapse=","),")"))),
                           eval(parse(text=paste0("list(",paste0("mi.grade.k5.varcov[",1:nimp,",,]",collapse=","),")"))))
  mi.grade.k4.pool=MIcombine(eval(parse(text=paste0("list(",paste0("mi.grade.k4.coef[",1:nimp,",]",collapse=","),")"))),
                           eval(parse(text=paste0("list(",paste0("mi.grade.k4.varcov[",1:nimp,",,]",collapse=","),")"))))
  mi.grade.k3.pool=MIcombine(eval(parse(text=paste0("list(",paste0("mi.grade.k3.coef[",1:nimp,",]",collapse=","),")"))),
                           eval(parse(text=paste0("list(",paste0("mi.grade.k3.varcov[",1:nimp,",,]",collapse=","),")"))))
  mi.grade.lin.pool=MIcombine(eval(parse(text=paste0("list(",paste0("mi.grade.lin.coef[",1:nimp,",]",collapse=","),")"))),
                            eval(parse(text=paste0("list(",paste0("mi.grade.lin.varcov[",1:nimp,",,]",collapse=","),")"))))
  
  num.param=4
  pool.coef=mi.grade.k5.pool$coefficients
  all.coef=mi.grade.k5.coef
  all.var=mi.grade.k5.varcov
  source("H:/MI_timedep_Cox/rotterdam_example/ruth_wald_test_mi_rotterdam_v2.R")
  wald.mi.grade.k5=pval
  stat.mi.grade.k5=d.x
  
  num.param=3
  pool.coef=mi.grade.k4.pool$coefficients
  all.coef=mi.grade.k4.coef
  all.var=mi.grade.k4.varcov
  source("H:/MI_timedep_Cox/rotterdam_example/ruth_wald_test_mi_rotterdam_v2.R")
  wald.mi.grade.k4=pval
  stat.mi.grade.k4=d.x
  
  num.param=2
  pool.coef=mi.grade.k3.pool$coefficients
  all.coef=mi.grade.k3.coef
  all.var=mi.grade.k3.varcov
  source("H:/MI_timedep_Cox/rotterdam_example/ruth_wald_test_mi_rotterdam_v2.R")
  wald.mi.grade.k3=pval
  stat.mi.grade.k3=d.x
  
  wald.mi.grade.lin=wald.test(coef=mi.grade.lin.pool$coef,Sigma=mi.grade.lin.pool$var,print.coef = F)$p.value
  stat.mi.grade.lin=wald.test(coef=mi.grade.lin.pool$coef,Sigma=mi.grade.lin.pool$var,print.coef = F)$statistic
  
  # waldtest.mi.grade.k5=wald.test(coef=mi.grade.k5.pool$coef,Sigma=mi.grade.k5.pool$var,print.coef = F)
  # waldtest.mi.grade.k4=wald.test(coef=mi.grade.k4.pool$coef,Sigma=mi.grade.k4.pool$var,print.coef = F)
  # waldtest.mi.grade.k3=wald.test(coef=mi.grade.k3.pool$coef,Sigma=mi.grade.k3.pool$var,print.coef = F)
  # waldtest.mi.grade.lin=wald.test(coef=mi.grade.lin.pool$coef,Sigma=mi.grade.lin.pool$var,print.coef = F)
  # 
  # wald.mi.grade.k5=waldtest.mi.grade.k5$p.value
  # wald.mi.grade.k4=waldtest.mi.grade.k4$p.value
  # wald.mi.grade.k3=waldtest.mi.grade.k3$p.value
  # wald.mi.grade.lin=waldtest.mi.grade.lin$p.value
  # 
  # stat.mi.grade.k5=waldtest.mi.grade.k5$statistic
  # stat.mi.grade.k4=waldtest.mi.grade.k4$statistic
  # stat.mi.grade.k3=waldtest.mi.grade.k3$statistic
  # stat.mi.grade.lin=waldtest.mi.grade.lin$statistic
  
}

#-------------
#enodessq

if("enodessq"%in%vars.remaining){
  mi.enodessq.k5.pool=MIcombine(eval(parse(text=paste0("list(",paste0("mi.enodessq.k5.coef[",1:nimp,",]",collapse=","),")"))),
                           eval(parse(text=paste0("list(",paste0("mi.enodessq.k5.varcov[",1:nimp,",,]",collapse=","),")"))))
  mi.enodessq.k4.pool=MIcombine(eval(parse(text=paste0("list(",paste0("mi.enodessq.k4.coef[",1:nimp,",]",collapse=","),")"))),
                           eval(parse(text=paste0("list(",paste0("mi.enodessq.k4.varcov[",1:nimp,",,]",collapse=","),")"))))
  mi.enodessq.k3.pool=MIcombine(eval(parse(text=paste0("list(",paste0("mi.enodessq.k3.coef[",1:nimp,",]",collapse=","),")"))),
                           eval(parse(text=paste0("list(",paste0("mi.enodessq.k3.varcov[",1:nimp,",,]",collapse=","),")"))))
  mi.enodessq.lin.pool=MIcombine(eval(parse(text=paste0("list(",paste0("mi.enodessq.lin.coef[",1:nimp,",]",collapse=","),")"))),
                            eval(parse(text=paste0("list(",paste0("mi.enodessq.lin.varcov[",1:nimp,",,]",collapse=","),")"))))
  
  num.param=4
  pool.coef=mi.enodessq.k5.pool$coefficients
  all.coef=mi.enodessq.k5.coef
  all.var=mi.enodessq.k5.varcov
  source("H:/MI_timedep_Cox/rotterdam_example/ruth_wald_test_mi_rotterdam_v2.R")
  wald.mi.enodessq.k5=pval
  stat.mi.enodessq.k5=d.x
  
  num.param=3
  pool.coef=mi.enodessq.k4.pool$coefficients
  all.coef=mi.enodessq.k4.coef
  all.var=mi.enodessq.k4.varcov
  source("H:/MI_timedep_Cox/rotterdam_example/ruth_wald_test_mi_rotterdam_v2.R")
  wald.mi.enodessq.k4=pval
  stat.mi.enodessq.k4=d.x
  
  num.param=2
  pool.coef=mi.enodessq.k3.pool$coefficients
  all.coef=mi.enodessq.k3.coef
  all.var=mi.enodessq.k3.varcov
  source("H:/MI_timedep_Cox/rotterdam_example/ruth_wald_test_mi_rotterdam_v2.R")
  wald.mi.enodessq.k3=pval
  stat.mi.enodessq.k3=d.x
  
  wald.mi.enodessq.lin=wald.test(coef=mi.enodessq.lin.pool$coef,Sigma=mi.enodessq.lin.pool$var,print.coef = F)$p.value
  stat.mi.enodessq.lin=wald.test(coef=mi.enodessq.lin.pool$coef,Sigma=mi.enodessq.lin.pool$var,print.coef = F)$statistic
  
  # waldtest.mi.enodessq.k5=wald.test(coef=mi.enodessq.k5.pool$coef,Sigma=mi.enodessq.k5.pool$var,print.coef = F)
  # waldtest.mi.enodessq.k4=wald.test(coef=mi.enodessq.k4.pool$coef,Sigma=mi.enodessq.k4.pool$var,print.coef = F)
  # waldtest.mi.enodessq.k3=wald.test(coef=mi.enodessq.k3.pool$coef,Sigma=mi.enodessq.k3.pool$var,print.coef = F)
  # waldtest.mi.enodessq.lin=wald.test(coef=mi.enodessq.lin.pool$coef,Sigma=mi.enodessq.lin.pool$var,print.coef = F)
  # 
  # wald.mi.enodessq.k5=waldtest.mi.enodessq.k5$p.value
  # wald.mi.enodessq.k4=waldtest.mi.enodessq.k4$p.value
  # wald.mi.enodessq.k3=waldtest.mi.enodessq.k3$p.value
  # wald.mi.enodessq.lin=waldtest.mi.enodessq.lin$p.value
  # 
  # stat.mi.enodessq.k5=waldtest.mi.enodessq.k5$statistic
  # stat.mi.enodessq.k4=waldtest.mi.enodessq.k4$statistic
  # stat.mi.enodessq.k3=waldtest.mi.enodessq.k3$statistic
  # stat.mi.enodessq.lin=waldtest.mi.enodessq.lin$statistic
  
}

#-------------
#hormon

if("hormon"%in%vars.remaining){
  mi.hormon.k5.pool=MIcombine(eval(parse(text=paste0("list(",paste0("mi.hormon.k5.coef[",1:nimp,",]",collapse=","),")"))),
                           eval(parse(text=paste0("list(",paste0("mi.hormon.k5.varcov[",1:nimp,",,]",collapse=","),")"))))
  mi.hormon.k4.pool=MIcombine(eval(parse(text=paste0("list(",paste0("mi.hormon.k4.coef[",1:nimp,",]",collapse=","),")"))),
                           eval(parse(text=paste0("list(",paste0("mi.hormon.k4.varcov[",1:nimp,",,]",collapse=","),")"))))
  mi.hormon.k3.pool=MIcombine(eval(parse(text=paste0("list(",paste0("mi.hormon.k3.coef[",1:nimp,",]",collapse=","),")"))),
                           eval(parse(text=paste0("list(",paste0("mi.hormon.k3.varcov[",1:nimp,",,]",collapse=","),")"))))
  mi.hormon.lin.pool=MIcombine(eval(parse(text=paste0("list(",paste0("mi.hormon.lin.coef[",1:nimp,",]",collapse=","),")"))),
                            eval(parse(text=paste0("list(",paste0("mi.hormon.lin.varcov[",1:nimp,",,]",collapse=","),")"))))
  
  num.param=4
  pool.coef=mi.hormon.k5.pool$coefficients
  all.coef=mi.hormon.k5.coef
  all.var=mi.hormon.k5.varcov
  source("H:/MI_timedep_Cox/rotterdam_example/ruth_wald_test_mi_rotterdam_v2.R")
  wald.mi.hormon.k5=pval
  stat.mi.hormon.k5=d.x
  
  num.param=3
  pool.coef=mi.hormon.k4.pool$coefficients
  all.coef=mi.hormon.k4.coef
  all.var=mi.hormon.k4.varcov
  source("H:/MI_timedep_Cox/rotterdam_example/ruth_wald_test_mi_rotterdam_v2.R")
  wald.mi.hormon.k4=pval
  stat.mi.hormon.k4=d.x
  
  num.param=2
  pool.coef=mi.hormon.k3.pool$coefficients
  all.coef=mi.hormon.k3.coef
  all.var=mi.hormon.k3.varcov
  source("H:/MI_timedep_Cox/rotterdam_example/ruth_wald_test_mi_rotterdam_v2.R")
  wald.mi.hormon.k3=pval
  stat.mi.hormon.k3=d.x
  
  wald.mi.hormon.lin=wald.test(coef=mi.hormon.lin.pool$coef,Sigma=mi.hormon.lin.pool$var,print.coef = F)$p.value
  stat.mi.hormon.lin=wald.test(coef=mi.hormon.lin.pool$coef,Sigma=mi.hormon.lin.pool$var,print.coef = F)$statistic
  
  # waldtest.mi.hormon.k5=wald.test(coef=mi.hormon.k5.pool$coef,Sigma=mi.hormon.k5.pool$var,print.coef = F)
  # waldtest.mi.hormon.k4=wald.test(coef=mi.hormon.k4.pool$coef,Sigma=mi.hormon.k4.pool$var,print.coef = F)
  # waldtest.mi.hormon.k3=wald.test(coef=mi.hormon.k3.pool$coef,Sigma=mi.hormon.k3.pool$var,print.coef = F)
  # waldtest.mi.hormon.lin=wald.test(coef=mi.hormon.lin.pool$coef,Sigma=mi.hormon.lin.pool$var,print.coef = F)
  # 
  # wald.mi.hormon.k5=waldtest.mi.hormon.k5$p.value
  # wald.mi.hormon.k4=waldtest.mi.hormon.k4$p.value
  # wald.mi.hormon.k3=waldtest.mi.hormon.k3$p.value
  # wald.mi.hormon.lin=waldtest.mi.hormon.lin$p.value
  # 
  # stat.mi.hormon.k5=waldtest.mi.hormon.k5$statistic
  # stat.mi.hormon.k4=waldtest.mi.hormon.k4$statistic
  # stat.mi.hormon.k3=waldtest.mi.hormon.k3$statistic
  # stat.mi.hormon.lin=waldtest.mi.hormon.lin$statistic
  
}

#-------------
#chemo

if("chemo"%in%vars.remaining){
  mi.chemo.k5.pool=MIcombine(eval(parse(text=paste0("list(",paste0("mi.chemo.k5.coef[",1:nimp,",]",collapse=","),")"))),
                           eval(parse(text=paste0("list(",paste0("mi.chemo.k5.varcov[",1:nimp,",,]",collapse=","),")"))))
  mi.chemo.k4.pool=MIcombine(eval(parse(text=paste0("list(",paste0("mi.chemo.k4.coef[",1:nimp,",]",collapse=","),")"))),
                           eval(parse(text=paste0("list(",paste0("mi.chemo.k4.varcov[",1:nimp,",,]",collapse=","),")"))))
  mi.chemo.k3.pool=MIcombine(eval(parse(text=paste0("list(",paste0("mi.chemo.k3.coef[",1:nimp,",]",collapse=","),")"))),
                           eval(parse(text=paste0("list(",paste0("mi.chemo.k3.varcov[",1:nimp,",,]",collapse=","),")"))))
  mi.chemo.lin.pool=MIcombine(eval(parse(text=paste0("list(",paste0("mi.chemo.lin.coef[",1:nimp,",]",collapse=","),")"))),
                            eval(parse(text=paste0("list(",paste0("mi.chemo.lin.varcov[",1:nimp,",,]",collapse=","),")"))))
  
  num.param=4
  pool.coef=mi.chemo.k5.pool$coefficients
  all.coef=mi.chemo.k5.coef
  all.var=mi.chemo.k5.varcov
  source("H:/MI_timedep_Cox/rotterdam_example/ruth_wald_test_mi_rotterdam_v2.R")
  wald.mi.chemo.k5=pval
  stat.mi.chemo.k5=d.x
  
  num.param=3
  pool.coef=mi.chemo.k4.pool$coefficients
  all.coef=mi.chemo.k4.coef
  all.var=mi.chemo.k4.varcov
  source("H:/MI_timedep_Cox/rotterdam_example/ruth_wald_test_mi_rotterdam_v2.R")
  wald.mi.chemo.k4=pval
  stat.mi.chemo.k4=d.x
  
  num.param=2
  pool.coef=mi.chemo.k3.pool$coefficients
  all.coef=mi.chemo.k3.coef
  all.var=mi.chemo.k3.varcov
  source("H:/MI_timedep_Cox/rotterdam_example/ruth_wald_test_mi_rotterdam_v2.R")
  wald.mi.chemo.k3=pval
  stat.mi.chemo.k3=d.x
  
  wald.mi.chemo.lin=wald.test(coef=mi.chemo.lin.pool$coef,Sigma=mi.chemo.lin.pool$var,print.coef = F)$p.value
  stat.mi.chemo.lin=wald.test(coef=mi.chemo.lin.pool$coef,Sigma=mi.chemo.lin.pool$var,print.coef = F)$statistic
  
  # waldtest.mi.chemo.k5=wald.test(coef=mi.chemo.k5.pool$coef,Sigma=mi.chemo.k5.pool$var,print.coef = F)
  # waldtest.mi.chemo.k4=wald.test(coef=mi.chemo.k4.pool$coef,Sigma=mi.chemo.k4.pool$var,print.coef = F)
  # waldtest.mi.chemo.k3=wald.test(coef=mi.chemo.k3.pool$coef,Sigma=mi.chemo.k3.pool$var,print.coef = F)
  # waldtest.mi.chemo.lin=wald.test(coef=mi.chemo.lin.pool$coef,Sigma=mi.chemo.lin.pool$var,print.coef = F)
  # 
  # wald.mi.chemo.k5=waldtest.mi.chemo.k5$p.value
  # wald.mi.chemo.k4=waldtest.mi.chemo.k4$p.value
  # wald.mi.chemo.k3=waldtest.mi.chemo.k3$p.value
  # wald.mi.chemo.lin=waldtest.mi.chemo.lin$p.value
  # 
  # stat.mi.chemo.k5=waldtest.mi.chemo.k5$statistic
  # stat.mi.chemo.k4=waldtest.mi.chemo.k4$statistic
  # stat.mi.chemo.k3=waldtest.mi.chemo.k3$statistic
  # stat.mi.chemo.lin=waldtest.mi.chemo.lin$statistic
  
}

#-------------
#logpgr

if("logpgr"%in%vars.remaining){
  mi.logpgr.k5.pool=MIcombine(eval(parse(text=paste0("list(",paste0("mi.logpgr.k5.coef[",1:nimp,",]",collapse=","),")"))),
                           eval(parse(text=paste0("list(",paste0("mi.logpgr.k5.varcov[",1:nimp,",,]",collapse=","),")"))))
  mi.logpgr.k4.pool=MIcombine(eval(parse(text=paste0("list(",paste0("mi.logpgr.k4.coef[",1:nimp,",]",collapse=","),")"))),
                           eval(parse(text=paste0("list(",paste0("mi.logpgr.k4.varcov[",1:nimp,",,]",collapse=","),")"))))
  mi.logpgr.k3.pool=MIcombine(eval(parse(text=paste0("list(",paste0("mi.logpgr.k3.coef[",1:nimp,",]",collapse=","),")"))),
                           eval(parse(text=paste0("list(",paste0("mi.logpgr.k3.varcov[",1:nimp,",,]",collapse=","),")"))))
  mi.logpgr.lin.pool=MIcombine(eval(parse(text=paste0("list(",paste0("mi.logpgr.lin.coef[",1:nimp,",]",collapse=","),")"))),
                            eval(parse(text=paste0("list(",paste0("mi.logpgr.lin.varcov[",1:nimp,",,]",collapse=","),")"))))
  
  num.param=4
  pool.coef=mi.logpgr.k5.pool$coefficients
  all.coef=mi.logpgr.k5.coef
  all.var=mi.logpgr.k5.varcov
  source("H:/MI_timedep_Cox/rotterdam_example/ruth_wald_test_mi_rotterdam_v2.R")
  wald.mi.logpgr.k5=pval
  stat.mi.logpgr.k5=d.x
  
  num.param=3
  pool.coef=mi.logpgr.k4.pool$coefficients
  all.coef=mi.logpgr.k4.coef
  all.var=mi.logpgr.k4.varcov
  source("H:/MI_timedep_Cox/rotterdam_example/ruth_wald_test_mi_rotterdam_v2.R")
  wald.mi.logpgr.k4=pval
  stat.mi.logpgr.k4=d.x
  
  num.param=2
  pool.coef=mi.logpgr.k3.pool$coefficients
  all.coef=mi.logpgr.k3.coef
  all.var=mi.logpgr.k3.varcov
  source("H:/MI_timedep_Cox/rotterdam_example/ruth_wald_test_mi_rotterdam_v2.R")
  wald.mi.logpgr.k3=pval
  stat.mi.logpgr.k3=d.x
  
  wald.mi.logpgr.lin=wald.test(coef=mi.logpgr.lin.pool$coef,Sigma=mi.logpgr.lin.pool$var,print.coef = F)$p.value
  stat.mi.logpgr.lin=wald.test(coef=mi.logpgr.lin.pool$coef,Sigma=mi.logpgr.lin.pool$var,print.coef = F)$statistic
  
  # waldtest.mi.logpgr.k5=wald.test(coef=mi.logpgr.k5.pool$coef,Sigma=mi.logpgr.k5.pool$var,print.coef = F)
  # waldtest.mi.logpgr.k4=wald.test(coef=mi.logpgr.k4.pool$coef,Sigma=mi.logpgr.k4.pool$var,print.coef = F)
  # waldtest.mi.logpgr.k3=wald.test(coef=mi.logpgr.k3.pool$coef,Sigma=mi.logpgr.k3.pool$var,print.coef = F)
  # waldtest.mi.logpgr.lin=wald.test(coef=mi.logpgr.lin.pool$coef,Sigma=mi.logpgr.lin.pool$var,print.coef = F)
  # 
  # wald.mi.logpgr.k5=waldtest.mi.logpgr.k5$p.value
  # wald.mi.logpgr.k4=waldtest.mi.logpgr.k4$p.value
  # wald.mi.logpgr.k3=waldtest.mi.logpgr.k3$p.value
  # wald.mi.logpgr.lin=waldtest.mi.logpgr.lin$p.value
  # 
  # stat.mi.logpgr.k5=waldtest.mi.logpgr.k5$statistic
  # stat.mi.logpgr.k4=waldtest.mi.logpgr.k4$statistic
  # stat.mi.logpgr.k3=waldtest.mi.logpgr.k3$statistic
  # stat.mi.logpgr.lin=waldtest.mi.logpgr.lin$statistic
  
}

#-------------
#table of results

wald.mi.age.whichmin=which.min(c(wald.mi.age.k5,wald.mi.age.k4,wald.mi.age.k3,wald.mi.age.lin))
wald.mi.sized1.whichmin=which.min(c(wald.mi.sized1.k5,wald.mi.sized1.k4,wald.mi.sized1.k3,wald.mi.sized1.lin))
wald.mi.sized2.whichmin=which.min(c(wald.mi.sized2.k5,wald.mi.sized2.k4,wald.mi.sized2.k3,wald.mi.sized2.lin))
wald.mi.grade.whichmin=which.min(c(wald.mi.grade.k5,wald.mi.grade.k4,wald.mi.grade.k3,wald.mi.grade.lin))
wald.mi.enodessq.whichmin=which.min(c(wald.mi.enodessq.k5,wald.mi.enodessq.k4,wald.mi.enodessq.k3,wald.mi.enodessq.lin))
wald.mi.hormon.whichmin=which.min(c(wald.mi.hormon.k5,wald.mi.hormon.k4,wald.mi.hormon.k3,wald.mi.hormon.lin))
wald.mi.chemo.whichmin=which.min(c(wald.mi.chemo.k5,wald.mi.chemo.k4,wald.mi.chemo.k3,wald.mi.chemo.lin))
wald.mi.logpgr.whichmin=which.min(c(wald.mi.logpgr.k5,wald.mi.logpgr.k4,wald.mi.logpgr.k3,wald.mi.logpgr.lin))

wald.mi.age.min=c(wald.mi.age.k5,wald.mi.age.k4,wald.mi.age.k3,wald.mi.age.lin)[wald.mi.age.whichmin]
wald.mi.sized1.min=c(wald.mi.sized1.k5,wald.mi.sized1.k4,wald.mi.sized1.k3,wald.mi.sized1.lin)[wald.mi.sized1.whichmin]
wald.mi.sized2.min=c(wald.mi.sized2.k5,wald.mi.sized2.k4,wald.mi.sized2.k3,wald.mi.sized2.lin)[wald.mi.sized2.whichmin]
wald.mi.grade.min=c(wald.mi.grade.k5,wald.mi.grade.k4,wald.mi.grade.k3,wald.mi.grade.lin)[wald.mi.grade.whichmin]
wald.mi.enodessq.min=c(wald.mi.enodessq.k5,wald.mi.enodessq.k4,wald.mi.enodessq.k3,wald.mi.enodessq.lin)[wald.mi.enodessq.whichmin]
wald.mi.hormon.min=c(wald.mi.hormon.k5,wald.mi.hormon.k4,wald.mi.hormon.k3,wald.mi.hormon.lin)[wald.mi.hormon.whichmin]
wald.mi.chemo.min=c(wald.mi.chemo.k5,wald.mi.chemo.k4,wald.mi.chemo.k3,wald.mi.chemo.lin)[wald.mi.chemo.whichmin]
wald.mi.logpgr.min=c(wald.mi.logpgr.k5,wald.mi.logpgr.k4,wald.mi.logpgr.k3,wald.mi.logpgr.lin)[wald.mi.logpgr.whichmin]

if(wald.mi.age.whichmin==1){wald.mi.age.whichmin="k5"} else 
  if(wald.mi.age.whichmin==2){wald.mi.age.whichmin="k4"} else
    if(wald.mi.age.whichmin==3){wald.mi.age.whichmin="k3"} else
      if(wald.mi.age.whichmin==4){wald.mi.age.whichmin="lin"} 

if(wald.mi.sized1.whichmin==1){wald.mi.sized1.whichmin="k5"} else 
  if(wald.mi.sized1.whichmin==2){wald.mi.sized1.whichmin="k4"} else
    if(wald.mi.sized1.whichmin==3){wald.mi.sized1.whichmin="k3"} else
      if(wald.mi.sized1.whichmin==4){wald.mi.sized1.whichmin="lin"} 

if(wald.mi.sized2.whichmin==1){wald.mi.sized2.whichmin="k5"} else 
  if(wald.mi.sized2.whichmin==2){wald.mi.sized2.whichmin="k4"} else
    if(wald.mi.sized2.whichmin==3){wald.mi.sized2.whichmin="k3"} else
      if(wald.mi.sized2.whichmin==4){wald.mi.sized2.whichmin="lin"} 

if(wald.mi.grade.whichmin==1){wald.mi.grade.whichmin="k5"} else 
  if(wald.mi.grade.whichmin==2){wald.mi.grade.whichmin="k4"} else
    if(wald.mi.grade.whichmin==3){wald.mi.grade.whichmin="k3"} else
      if(wald.mi.grade.whichmin==4){wald.mi.grade.whichmin="lin"} 

if(wald.mi.enodessq.whichmin==1){wald.mi.enodessq.whichmin="k5"} else 
  if(wald.mi.enodessq.whichmin==2){wald.mi.enodessq.whichmin="k4"} else
    if(wald.mi.enodessq.whichmin==3){wald.mi.enodessq.whichmin="k3"} else
      if(wald.mi.enodessq.whichmin==4){wald.mi.enodessq.whichmin="lin"} 

if(wald.mi.hormon.whichmin==1){wald.mi.hormon.whichmin="k5"} else 
  if(wald.mi.hormon.whichmin==2){wald.mi.hormon.whichmin="k4"} else
    if(wald.mi.hormon.whichmin==3){wald.mi.hormon.whichmin="k3"} else
      if(wald.mi.hormon.whichmin==4){wald.mi.hormon.whichmin="lin"} 

if(wald.mi.chemo.whichmin==1){wald.mi.chemo.whichmin="k5"} else 
  if(wald.mi.chemo.whichmin==2){wald.mi.chemo.whichmin="k4"} else
    if(wald.mi.chemo.whichmin==3){wald.mi.chemo.whichmin="k3"} else
      if(wald.mi.chemo.whichmin==4){wald.mi.chemo.whichmin="lin"} 

if(wald.mi.logpgr.whichmin==1){wald.mi.logpgr.whichmin="k5"} else 
  if(wald.mi.logpgr.whichmin==2){wald.mi.logpgr.whichmin="k4"} else
    if(wald.mi.logpgr.whichmin==3){wald.mi.logpgr.whichmin="k3"} else
      if(wald.mi.logpgr.whichmin==4){wald.mi.logpgr.whichmin="lin"} 


mi.res=rbind(c(wald.mi.age.k5,wald.mi.age.k4,wald.mi.age.k3,wald.mi.age.lin,wald.mi.age.whichmin,wald.mi.age.min),
  c(wald.mi.sized1.k5,wald.mi.sized1.k4,wald.mi.sized1.k3,wald.mi.sized1.lin,wald.mi.sized1.whichmin,wald.mi.sized1.min),
  c(wald.mi.sized2.k5,wald.mi.sized2.k4,wald.mi.sized2.k3,wald.mi.sized2.lin,wald.mi.sized2.whichmin,wald.mi.sized2.min),
  c(wald.mi.grade.k5,wald.mi.grade.k4,wald.mi.grade.k3,wald.mi.grade.lin,wald.mi.grade.whichmin,wald.mi.grade.min),
  c(wald.mi.enodessq.k5,wald.mi.enodessq.k4,wald.mi.enodessq.k3,wald.mi.enodessq.lin,wald.mi.enodessq.whichmin,wald.mi.enodessq.min),
  c(wald.mi.hormon.k5,wald.mi.hormon.k4,wald.mi.hormon.k3,wald.mi.hormon.lin,wald.mi.hormon.whichmin,wald.mi.hormon.min),
  c(wald.mi.chemo.k5,wald.mi.chemo.k4,wald.mi.chemo.k3,wald.mi.chemo.lin,wald.mi.chemo.whichmin,wald.mi.chemo.min),
  c(wald.mi.logpgr.k5,wald.mi.logpgr.k4,wald.mi.logpgr.k3,wald.mi.logpgr.lin,wald.mi.logpgr.whichmin,wald.mi.logpgr.min))

colnames(mi.res)=c("k5","k4","k3","lin","which.min","MIN")
rownames(mi.res)=c("age","sized1","sized2","grade","enodessq","hormon","chemo","logpgr")

mi.stat.res=rbind(c(stat.mi.age.k5,stat.mi.age.k4,stat.mi.age.k3,stat.mi.age.lin),
                         c(stat.mi.sized1.k5,stat.mi.sized1.k4,stat.mi.sized1.k3,stat.mi.sized1.lin),
                         c(stat.mi.sized2.k5,stat.mi.sized2.k4,stat.mi.sized2.k3,stat.mi.sized2.lin),
                         c(stat.mi.grade.k5,stat.mi.grade.k4,stat.mi.grade.k3,stat.mi.grade.lin),
                         c(stat.mi.enodessq.k5,stat.mi.enodessq.k4,stat.mi.enodessq.k3,stat.mi.enodessq.lin),
                         c(stat.mi.hormon.k5,stat.mi.hormon.k4,stat.mi.hormon.k3,stat.mi.hormon.lin),
                         c(stat.mi.chemo.k5,stat.mi.chemo.k4,stat.mi.chemo.k3,stat.mi.chemo.lin),
                         c(stat.mi.logpgr.k5,stat.mi.logpgr.k4,stat.mi.logpgr.k3,stat.mi.logpgr.lin))
colnames(mi.stat.res)=c("k5","k4","k3","lin")
rownames(mi.stat.res)=c("age","sized1","sized2","grade","enodessq","hormon","chemo","logpgr")

mi.res=mi.res[vars.remaining,]
mi.stat.res=mi.stat.res[vars.remaining,]

# mi.res
# 
# mi.stat.res

