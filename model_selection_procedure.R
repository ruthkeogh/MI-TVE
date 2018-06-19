
#-------------
#age

if("age"%in%vars.remaining){
  formula.tve.sp=as.formula(paste(formula.base,"+age:t+age:t.k5.1+age:t.k5.2+age:t.k5.3"))
  compdata.age.k5=coxph(formula.tve.sp,data=mydata.split)
  wald.age.k5=wald.test(coef=compdata.age.k5$coef[c("age:t","age:t.k5.1","age:t.k5.2","age:t.k5.3")],
            Sigma=vcov(compdata.age.k5)[c("age:t","age:t.k5.1","age:t.k5.2","age:t.k5.3"),c("age:t","age:t.k5.1","age:t.k5.2","age:t.k5.3")],print.coef = F)
  wald.compdata.age.k5=wald.age.k5$p.value
  stat.compdata.age.k5=wald.age.k5$statistic
  
  formula.tve.sp=as.formula(paste(formula.base,"+age:t+age:t.k4.1+age:t.k4.2"))
  compdata.age.k4=coxph(formula.tve.sp,data=mydata.split)
  wald.age.k4=wald.test(coef=compdata.age.k4$coef[c("age:t","age:t.k4.1","age:t.k4.2")],
                        Sigma=vcov(compdata.age.k4)[c("age:t","age:t.k4.1","age:t.k4.2"),c("age:t","age:t.k4.1","age:t.k4.2")],print.coef = F)
  wald.compdata.age.k4=wald.age.k4$p.value
  stat.compdata.age.k4=wald.age.k4$statistic
  
  formula.tve.sp=as.formula(paste(formula.base,"+age:t+age:t.k3"))
  compdata.age.k3=coxph(formula.tve.sp,data=mydata.split)
  wald.age.k3=wald.test(coef=compdata.age.k3$coef[c("age:t","age:t.k3")],
                        Sigma=vcov(compdata.age.k3)[c("age:t","age:t.k3"),c("age:t","age:t.k3")],print.coef = F)
  wald.compdata.age.k3=wald.age.k3$p.value
  stat.compdata.age.k3=wald.age.k3$statistic
  
  formula.tve.sp=as.formula(paste(formula.base,"+age:t"))
  compdata.age.lin=coxph(formula.tve.sp,data=mydata.split)
  wald.age.lin=wald.test(coef=compdata.age.lin$coef[c("age:t")],
                        Sigma=vcov(compdata.age.lin)[c("age:t"),c("age:t")],print.coef = F)
  wald.compdata.age.lin=wald.age.lin$p.value
  stat.compdata.age.lin=wald.age.lin$statistic
}

#-------------
#sized1

if("sized1"%in%vars.remaining){
  formula.tve.sp=as.formula(paste(formula.base,"+sized1:t+sized1:t.k5.1+sized1:t.k5.2+sized1:t.k5.3"))
  compdata.sized1.k5=coxph(formula.tve.sp,data=mydata.split)
  wald.sized1.k5=wald.test(coef=compdata.sized1.k5$coef[c("sized1:t","sized1:t.k5.1","sized1:t.k5.2","sized1:t.k5.3")],
                        Sigma=vcov(compdata.sized1.k5)[c("sized1:t","sized1:t.k5.1","sized1:t.k5.2","sized1:t.k5.3"),c("sized1:t","sized1:t.k5.1","sized1:t.k5.2","sized1:t.k5.3")],print.coef = F)
  wald.compdata.sized1.k5=wald.sized1.k5$p.value
  stat.compdata.sized1.k5=wald.sized1.k5$statistic
  
  formula.tve.sp=as.formula(paste(formula.base,"+sized1:t+sized1:t.k4.1+sized1:t.k4.2"))
  compdata.sized1.k4=coxph(formula.tve.sp,data=mydata.split)
  wald.sized1.k4=wald.test(coef=compdata.sized1.k4$coef[c("sized1:t","sized1:t.k4.1","sized1:t.k4.2")],
                        Sigma=vcov(compdata.sized1.k4)[c("sized1:t","sized1:t.k4.1","sized1:t.k4.2"),c("sized1:t","sized1:t.k4.1","sized1:t.k4.2")],print.coef = F)
  wald.compdata.sized1.k4=wald.sized1.k4$p.value
  stat.compdata.sized1.k4=wald.sized1.k4$statistic
  
  formula.tve.sp=as.formula(paste(formula.base,"+sized1:t+sized1:t.k3"))
  compdata.sized1.k3=coxph(formula.tve.sp,data=mydata.split)
  wald.sized1.k3=wald.test(coef=compdata.sized1.k3$coef[c("sized1:t","sized1:t.k3")],
                        Sigma=vcov(compdata.sized1.k3)[c("sized1:t","sized1:t.k3"),c("sized1:t","sized1:t.k3")],print.coef = F)
  wald.compdata.sized1.k3=wald.sized1.k3$p.value
  stat.compdata.sized1.k3=wald.sized1.k3$statistic
  
  formula.tve.sp=as.formula(paste(formula.base,"+sized1:t"))
  compdata.sized1.lin=coxph(formula.tve.sp,data=mydata.split)
  wald.sized1.lin=wald.test(coef=compdata.sized1.lin$coef[c("sized1:t")],
                         Sigma=vcov(compdata.sized1.lin)[c("sized1:t"),c("sized1:t")],print.coef = F)
  wald.compdata.sized1.lin=wald.sized1.lin$p.value
  stat.compdata.sized1.lin=wald.sized1.lin$statistic
}

#-------------
#sized2

if("sized2"%in%vars.remaining){
  formula.tve.sp=as.formula(paste(formula.base,"+sized2:t+sized2:t.k5.1+sized2:t.k5.2+sized2:t.k5.3"))
  compdata.sized2.k5=coxph(formula.tve.sp,data=mydata.split)
  wald.sized2.k5=wald.test(coef=compdata.sized2.k5$coef[c("sized2:t","sized2:t.k5.1","sized2:t.k5.2","sized2:t.k5.3")],
                        Sigma=vcov(compdata.sized2.k5)[c("sized2:t","sized2:t.k5.1","sized2:t.k5.2","sized2:t.k5.3"),c("sized2:t","sized2:t.k5.1","sized2:t.k5.2","sized2:t.k5.3")],print.coef = F)
  wald.compdata.sized2.k5=wald.sized2.k5$p.value
  stat.compdata.sized2.k5=wald.sized2.k5$statistic
  
  formula.tve.sp=as.formula(paste(formula.base,"+sized2:t+sized2:t.k4.1+sized2:t.k4.2"))
  compdata.sized2.k4=coxph(formula.tve.sp,data=mydata.split)
  wald.sized2.k4=wald.test(coef=compdata.sized2.k4$coef[c("sized2:t","sized2:t.k4.1","sized2:t.k4.2")],
                        Sigma=vcov(compdata.sized2.k4)[c("sized2:t","sized2:t.k4.1","sized2:t.k4.2"),c("sized2:t","sized2:t.k4.1","sized2:t.k4.2")],print.coef = F)
  wald.compdata.sized2.k4=wald.sized2.k4$p.value
  stat.compdata.sized2.k4=wald.sized2.k4$statistic
  
  formula.tve.sp=as.formula(paste(formula.base,"+sized2:t+sized2:t.k3"))
  compdata.sized2.k3=coxph(formula.tve.sp,data=mydata.split)
  wald.sized2.k3=wald.test(coef=compdata.sized2.k3$coef[c("sized2:t","sized2:t.k3")],
                        Sigma=vcov(compdata.sized2.k3)[c("sized2:t","sized2:t.k3"),c("sized2:t","sized2:t.k3")],print.coef = F)
  wald.compdata.sized2.k3=wald.sized2.k3$p.value
  stat.compdata.sized2.k3=wald.sized2.k3$statistic
  
  formula.tve.sp=as.formula(paste(formula.base,"+sized2:t"))
  compdata.sized2.lin=coxph(formula.tve.sp,data=mydata.split)
  wald.sized2.lin=wald.test(coef=compdata.sized2.lin$coef[c("sized2:t")],
                         Sigma=vcov(compdata.sized2.lin)[c("sized2:t"),c("sized2:t")],print.coef = F)
  wald.compdata.sized2.lin=wald.sized2.lin$p.value
  stat.compdata.sized2.lin=wald.sized2.lin$statistic
}

#-------------
#grade

if("grade"%in%vars.remaining){
  formula.tve.sp=as.formula(paste(formula.base,"+grade:t+grade:t.k5.1+grade:t.k5.2+grade:t.k5.3"))
  compdata.grade.k5=coxph(formula.tve.sp,data=mydata.split)
  wald.grade.k5=wald.test(coef=compdata.grade.k5$coef[c("grade:t","grade:t.k5.1","grade:t.k5.2","grade:t.k5.3")],
                        Sigma=vcov(compdata.grade.k5)[c("grade:t","grade:t.k5.1","grade:t.k5.2","grade:t.k5.3"),c("grade:t","grade:t.k5.1","grade:t.k5.2","grade:t.k5.3")],print.coef = F)
  wald.compdata.grade.k5=wald.grade.k5$p.value
  stat.compdata.grade.k5=wald.grade.k5$statistic
  
  formula.tve.sp=as.formula(paste(formula.base,"+grade:t+grade:t.k4.1+grade:t.k4.2"))
  compdata.grade.k4=coxph(formula.tve.sp,data=mydata.split)
  wald.grade.k4=wald.test(coef=compdata.grade.k4$coef[c("grade:t","grade:t.k4.1","grade:t.k4.2")],
                        Sigma=vcov(compdata.grade.k4)[c("grade:t","grade:t.k4.1","grade:t.k4.2"),c("grade:t","grade:t.k4.1","grade:t.k4.2")],print.coef = F)
  wald.compdata.grade.k4=wald.grade.k4$p.value
  stat.compdata.grade.k4=wald.grade.k4$statistic
  
  formula.tve.sp=as.formula(paste(formula.base,"+grade:t+grade:t.k3"))
  compdata.grade.k3=coxph(formula.tve.sp,data=mydata.split)
  wald.grade.k3=wald.test(coef=compdata.grade.k3$coef[c("grade:t","grade:t.k3")],
                        Sigma=vcov(compdata.grade.k3)[c("grade:t","grade:t.k3"),c("grade:t","grade:t.k3")],print.coef = F)
  wald.compdata.grade.k3=wald.grade.k3$p.value
  stat.compdata.grade.k3=wald.grade.k3$statistic
  
  formula.tve.sp=as.formula(paste(formula.base,"+grade:t"))
  compdata.grade.lin=coxph(formula.tve.sp,data=mydata.split)
  wald.grade.lin=wald.test(coef=compdata.grade.lin$coef[c("grade:t")],
                         Sigma=vcov(compdata.grade.lin)[c("grade:t"),c("grade:t")],print.coef = F)
  wald.compdata.grade.lin=wald.grade.lin$p.value
  stat.compdata.grade.lin=wald.grade.lin$statistic
}

#-------------
#enodessq

if("enodessq"%in%vars.remaining){
  formula.tve.sp=as.formula(paste(formula.base,"+enodessq:t+enodessq:t.k5.1+enodessq:t.k5.2+enodessq:t.k5.3"))
  compdata.enodessq.k5=coxph(formula.tve.sp,data=mydata.split)
  wald.enodessq.k5=wald.test(coef=compdata.enodessq.k5$coef[c("enodessq:t","enodessq:t.k5.1","enodessq:t.k5.2","enodessq:t.k5.3")],
                        Sigma=vcov(compdata.enodessq.k5)[c("enodessq:t","enodessq:t.k5.1","enodessq:t.k5.2","enodessq:t.k5.3"),c("enodessq:t","enodessq:t.k5.1","enodessq:t.k5.2","enodessq:t.k5.3")],print.coef = F)
  wald.compdata.enodessq.k5=wald.enodessq.k5$p.value
  stat.compdata.enodessq.k5=wald.enodessq.k5$statistic
  
  formula.tve.sp=as.formula(paste(formula.base,"+enodessq:t+enodessq:t.k4.1+enodessq:t.k4.2"))
  compdata.enodessq.k4=coxph(formula.tve.sp,data=mydata.split)
  wald.enodessq.k4=wald.test(coef=compdata.enodessq.k4$coef[c("enodessq:t","enodessq:t.k4.1","enodessq:t.k4.2")],
                        Sigma=vcov(compdata.enodessq.k4)[c("enodessq:t","enodessq:t.k4.1","enodessq:t.k4.2"),c("enodessq:t","enodessq:t.k4.1","enodessq:t.k4.2")],print.coef = F)
  wald.compdata.enodessq.k4=wald.enodessq.k4$p.value
  stat.compdata.enodessq.k4=wald.enodessq.k4$statistic
  
  formula.tve.sp=as.formula(paste(formula.base,"+enodessq:t+enodessq:t.k3"))
  compdata.enodessq.k3=coxph(formula.tve.sp,data=mydata.split)
  wald.enodessq.k3=wald.test(coef=compdata.enodessq.k3$coef[c("enodessq:t","enodessq:t.k3")],
                        Sigma=vcov(compdata.enodessq.k3)[c("enodessq:t","enodessq:t.k3"),c("enodessq:t","enodessq:t.k3")],print.coef = F)
  wald.compdata.enodessq.k3=wald.enodessq.k3$p.value
  stat.compdata.enodessq.k3=wald.enodessq.k3$statistic
  
  formula.tve.sp=as.formula(paste(formula.base,"+enodessq:t"))
  compdata.enodessq.lin=coxph(formula.tve.sp,data=mydata.split)
  wald.enodessq.lin=wald.test(coef=compdata.enodessq.lin$coef[c("enodessq:t")],
                         Sigma=vcov(compdata.enodessq.lin)[c("enodessq:t"),c("enodessq:t")],print.coef = F)
  wald.compdata.enodessq.lin=wald.enodessq.lin$p.value
  stat.compdata.enodessq.lin=wald.enodessq.lin$statistic
}

#-------------
#hormon

if("hormon"%in%vars.remaining){
  formula.tve.sp=as.formula(paste(formula.base,"+hormon:t+hormon:t.k5.1+hormon:t.k5.2+hormon:t.k5.3"))
  compdata.hormon.k5=coxph(formula.tve.sp,data=mydata.split)
  wald.hormon.k5=wald.test(coef=compdata.hormon.k5$coef[c("hormon:t","hormon:t.k5.1","hormon:t.k5.2","hormon:t.k5.3")],
                        Sigma=vcov(compdata.hormon.k5)[c("hormon:t","hormon:t.k5.1","hormon:t.k5.2","hormon:t.k5.3"),c("hormon:t","hormon:t.k5.1","hormon:t.k5.2","hormon:t.k5.3")],print.coef = F)
  wald.compdata.hormon.k5=wald.hormon.k5$p.value
  stat.compdata.hormon.k5=wald.hormon.k5$statistic
  
  formula.tve.sp=as.formula(paste(formula.base,"+hormon:t+hormon:t.k4.1+hormon:t.k4.2"))
  compdata.hormon.k4=coxph(formula.tve.sp,data=mydata.split)
  wald.hormon.k4=wald.test(coef=compdata.hormon.k4$coef[c("hormon:t","hormon:t.k4.1","hormon:t.k4.2")],
                        Sigma=vcov(compdata.hormon.k4)[c("hormon:t","hormon:t.k4.1","hormon:t.k4.2"),c("hormon:t","hormon:t.k4.1","hormon:t.k4.2")],print.coef = F)
  wald.compdata.hormon.k4=wald.hormon.k4$p.value
  stat.compdata.hormon.k4=wald.hormon.k4$statistic
  
  formula.tve.sp=as.formula(paste(formula.base,"+hormon:t+hormon:t.k3"))
  compdata.hormon.k3=coxph(formula.tve.sp,data=mydata.split)
  wald.hormon.k3=wald.test(coef=compdata.hormon.k3$coef[c("hormon:t","hormon:t.k3")],
                        Sigma=vcov(compdata.hormon.k3)[c("hormon:t","hormon:t.k3"),c("hormon:t","hormon:t.k3")],print.coef = F)
  wald.compdata.hormon.k3=wald.hormon.k3$p.value
  stat.compdata.hormon.k3=wald.hormon.k3$statistic
  
  formula.tve.sp=as.formula(paste(formula.base,"+hormon:t"))
  compdata.hormon.lin=coxph(formula.tve.sp,data=mydata.split)
  wald.hormon.lin=wald.test(coef=compdata.hormon.lin$coef[c("hormon:t")],
                         Sigma=vcov(compdata.hormon.lin)[c("hormon:t"),c("hormon:t")],print.coef = F)
  wald.compdata.hormon.lin=wald.hormon.lin$p.value
  stat.compdata.hormon.lin=wald.hormon.lin$statistic
}

#-------------
#chemo

if("chemo"%in%vars.remaining){
  formula.tve.sp=as.formula(paste(formula.base,"+chemo:t+chemo:t.k5.1+chemo:t.k5.2+chemo:t.k5.3"))
  compdata.chemo.k5=coxph(formula.tve.sp,data=mydata.split)
  wald.chemo.k5=wald.test(coef=compdata.chemo.k5$coef[c("chemo:t","chemo:t.k5.1","chemo:t.k5.2","chemo:t.k5.3")],
                        Sigma=vcov(compdata.chemo.k5)[c("chemo:t","chemo:t.k5.1","chemo:t.k5.2","chemo:t.k5.3"),c("chemo:t","chemo:t.k5.1","chemo:t.k5.2","chemo:t.k5.3")],print.coef = F)
  wald.compdata.chemo.k5=wald.chemo.k5$p.value
  stat.compdata.chemo.k5=wald.chemo.k5$statistic
  
  formula.tve.sp=as.formula(paste(formula.base,"+chemo:t+chemo:t.k4.1+chemo:t.k4.2"))
  compdata.chemo.k4=coxph(formula.tve.sp,data=mydata.split)
  wald.chemo.k4=wald.test(coef=compdata.chemo.k4$coef[c("chemo:t","chemo:t.k4.1","chemo:t.k4.2")],
                        Sigma=vcov(compdata.chemo.k4)[c("chemo:t","chemo:t.k4.1","chemo:t.k4.2"),c("chemo:t","chemo:t.k4.1","chemo:t.k4.2")],print.coef = F)
  wald.compdata.chemo.k4=wald.chemo.k4$p.value
  stat.compdata.chemo.k4=wald.chemo.k4$statistic
  
  formula.tve.sp=as.formula(paste(formula.base,"+chemo:t+chemo:t.k3"))
  compdata.chemo.k3=coxph(formula.tve.sp,data=mydata.split)
  wald.chemo.k3=wald.test(coef=compdata.chemo.k3$coef[c("chemo:t","chemo:t.k3")],
                        Sigma=vcov(compdata.chemo.k3)[c("chemo:t","chemo:t.k3"),c("chemo:t","chemo:t.k3")],print.coef = F)
  wald.compdata.chemo.k3=wald.chemo.k3$p.value
  stat.compdata.chemo.k3=wald.chemo.k3$statistic
  
  formula.tve.sp=as.formula(paste(formula.base,"+chemo:t"))
  compdata.chemo.lin=coxph(formula.tve.sp,data=mydata.split)
  wald.chemo.lin=wald.test(coef=compdata.chemo.lin$coef[c("chemo:t")],
                         Sigma=vcov(compdata.chemo.lin)[c("chemo:t"),c("chemo:t")],print.coef = F)
  wald.compdata.chemo.lin=wald.chemo.lin$p.value
  stat.compdata.chemo.lin=wald.chemo.lin$statistic
}

#-------------
#logpgr

if("logpgr"%in%vars.remaining){
  formula.tve.sp=as.formula(paste(formula.base,"+logpgr:t+logpgr:t.k5.1+logpgr:t.k5.2+logpgr:t.k5.3"))
  compdata.logpgr.k5=coxph(formula.tve.sp,data=mydata.split)
  wald.logpgr.k5=wald.test(coef=compdata.logpgr.k5$coef[c("logpgr:t","logpgr:t.k5.1","logpgr:t.k5.2","logpgr:t.k5.3")],
                        Sigma=vcov(compdata.logpgr.k5)[c("logpgr:t","logpgr:t.k5.1","logpgr:t.k5.2","logpgr:t.k5.3"),c("logpgr:t","logpgr:t.k5.1","logpgr:t.k5.2","logpgr:t.k5.3")],print.coef = F)
  wald.compdata.logpgr.k5=wald.logpgr.k5$p.value
  stat.compdata.logpgr.k5=wald.logpgr.k5$statistic
  
  formula.tve.sp=as.formula(paste(formula.base,"+logpgr:t+logpgr:t.k4.1+logpgr:t.k4.2"))
  compdata.logpgr.k4=coxph(formula.tve.sp,data=mydata.split)
  wald.logpgr.k4=wald.test(coef=compdata.logpgr.k4$coef[c("logpgr:t","logpgr:t.k4.1","logpgr:t.k4.2")],
                        Sigma=vcov(compdata.logpgr.k4)[c("logpgr:t","logpgr:t.k4.1","logpgr:t.k4.2"),c("logpgr:t","logpgr:t.k4.1","logpgr:t.k4.2")],print.coef = F)
  wald.compdata.logpgr.k4=wald.logpgr.k4$p.value
  stat.compdata.logpgr.k4=wald.logpgr.k4$statistic
  
  formula.tve.sp=as.formula(paste(formula.base,"+logpgr:t+logpgr:t.k3"))
  compdata.logpgr.k3=coxph(formula.tve.sp,data=mydata.split)
  wald.logpgr.k3=wald.test(coef=compdata.logpgr.k3$coef[c("logpgr:t","logpgr:t.k3")],
                        Sigma=vcov(compdata.logpgr.k3)[c("logpgr:t","logpgr:t.k3"),c("logpgr:t","logpgr:t.k3")],print.coef = F)
  wald.compdata.logpgr.k3=wald.logpgr.k3$p.value
  stat.compdata.logpgr.k3=wald.logpgr.k3$statistic
  
  formula.tve.sp=as.formula(paste(formula.base,"+logpgr:t"))
  compdata.logpgr.lin=coxph(formula.tve.sp,data=mydata.split)
  wald.logpgr.lin=wald.test(coef=compdata.logpgr.lin$coef[c("logpgr:t")],
                         Sigma=vcov(compdata.logpgr.lin)[c("logpgr:t"),c("logpgr:t")],print.coef = F)
  wald.compdata.logpgr.lin=wald.logpgr.lin$p.value
  stat.compdata.logpgr.lin=wald.logpgr.lin$statistic
}

#-------------
#table of results

wald.compdata.age.whichmin=which.min(c(wald.compdata.age.k5,wald.compdata.age.k4,wald.compdata.age.k3,wald.compdata.age.lin))
wald.compdata.sized1.whichmin=which.min(c(wald.compdata.sized1.k5,wald.compdata.sized1.k4,wald.compdata.sized1.k3,wald.compdata.sized1.lin))
wald.compdata.sized2.whichmin=which.min(c(wald.compdata.sized2.k5,wald.compdata.sized2.k4,wald.compdata.sized2.k3,wald.compdata.sized2.lin))
wald.compdata.grade.whichmin=which.min(c(wald.compdata.grade.k5,wald.compdata.grade.k4,wald.compdata.grade.k3,wald.compdata.grade.lin))
wald.compdata.enodessq.whichmin=which.min(c(wald.compdata.enodessq.k5,wald.compdata.enodessq.k4,wald.compdata.enodessq.k3,wald.compdata.enodessq.lin))
wald.compdata.hormon.whichmin=which.min(c(wald.compdata.hormon.k5,wald.compdata.hormon.k4,wald.compdata.hormon.k3,wald.compdata.hormon.lin))
wald.compdata.chemo.whichmin=which.min(c(wald.compdata.chemo.k5,wald.compdata.chemo.k4,wald.compdata.chemo.k3,wald.compdata.chemo.lin))
wald.compdata.logpgr.whichmin=which.min(c(wald.compdata.logpgr.k5,wald.compdata.logpgr.k4,wald.compdata.logpgr.k3,wald.compdata.logpgr.lin))

wald.compdata.age.min=c(wald.compdata.age.k5,wald.compdata.age.k4,wald.compdata.age.k3,wald.compdata.age.lin)[wald.compdata.age.whichmin]
wald.compdata.sized1.min=c(wald.compdata.sized1.k5,wald.compdata.sized1.k4,wald.compdata.sized1.k3,wald.compdata.sized1.lin)[wald.compdata.sized1.whichmin]
wald.compdata.sized2.min=c(wald.compdata.sized2.k5,wald.compdata.sized2.k4,wald.compdata.sized2.k3,wald.compdata.sized2.lin)[wald.compdata.sized2.whichmin]
wald.compdata.grade.min=c(wald.compdata.grade.k5,wald.compdata.grade.k4,wald.compdata.grade.k3,wald.compdata.grade.lin)[wald.compdata.grade.whichmin]
wald.compdata.enodessq.min=c(wald.compdata.enodessq.k5,wald.compdata.enodessq.k4,wald.compdata.enodessq.k3,wald.compdata.enodessq.lin)[wald.compdata.enodessq.whichmin]
wald.compdata.hormon.min=c(wald.compdata.hormon.k5,wald.compdata.hormon.k4,wald.compdata.hormon.k3,wald.compdata.hormon.lin)[wald.compdata.hormon.whichmin]
wald.compdata.chemo.min=c(wald.compdata.chemo.k5,wald.compdata.chemo.k4,wald.compdata.chemo.k3,wald.compdata.chemo.lin)[wald.compdata.chemo.whichmin]
wald.compdata.logpgr.min=c(wald.compdata.logpgr.k5,wald.compdata.logpgr.k4,wald.compdata.logpgr.k3,wald.compdata.logpgr.lin)[wald.compdata.logpgr.whichmin]

if(wald.compdata.age.whichmin==1){wald.compdata.age.whichmin="k5"} else 
  if(wald.compdata.age.whichmin==2){wald.compdata.age.whichmin="k4"} else
    if(wald.compdata.age.whichmin==3){wald.compdata.age.whichmin="k3"} else
      if(wald.compdata.age.whichmin==4){wald.compdata.age.whichmin="lin"} 

if(wald.compdata.sized1.whichmin==1){wald.compdata.sized1.whichmin="k5"} else 
  if(wald.compdata.sized1.whichmin==2){wald.compdata.sized1.whichmin="k4"} else
    if(wald.compdata.sized1.whichmin==3){wald.compdata.sized1.whichmin="k3"} else
      if(wald.compdata.sized1.whichmin==4){wald.compdata.sized1.whichmin="lin"} 

if(wald.compdata.sized2.whichmin==1){wald.compdata.sized2.whichmin="k5"} else 
  if(wald.compdata.sized2.whichmin==2){wald.compdata.sized2.whichmin="k4"} else
    if(wald.compdata.sized2.whichmin==3){wald.compdata.sized2.whichmin="k3"} else
      if(wald.compdata.sized2.whichmin==4){wald.compdata.sized2.whichmin="lin"} 

if(wald.compdata.grade.whichmin==1){wald.compdata.grade.whichmin="k5"} else 
  if(wald.compdata.grade.whichmin==2){wald.compdata.grade.whichmin="k4"} else
    if(wald.compdata.grade.whichmin==3){wald.compdata.grade.whichmin="k3"} else
      if(wald.compdata.grade.whichmin==4){wald.compdata.grade.whichmin="lin"} 

if(wald.compdata.enodessq.whichmin==1){wald.compdata.enodessq.whichmin="k5"} else 
  if(wald.compdata.enodessq.whichmin==2){wald.compdata.enodessq.whichmin="k4"} else
    if(wald.compdata.enodessq.whichmin==3){wald.compdata.enodessq.whichmin="k3"} else
      if(wald.compdata.enodessq.whichmin==4){wald.compdata.enodessq.whichmin="lin"} 

if(wald.compdata.hormon.whichmin==1){wald.compdata.hormon.whichmin="k5"} else 
  if(wald.compdata.hormon.whichmin==2){wald.compdata.hormon.whichmin="k4"} else
    if(wald.compdata.hormon.whichmin==3){wald.compdata.hormon.whichmin="k3"} else
      if(wald.compdata.hormon.whichmin==4){wald.compdata.hormon.whichmin="lin"} 

if(wald.compdata.chemo.whichmin==1){wald.compdata.chemo.whichmin="k5"} else 
  if(wald.compdata.chemo.whichmin==2){wald.compdata.chemo.whichmin="k4"} else
    if(wald.compdata.chemo.whichmin==3){wald.compdata.chemo.whichmin="k3"} else
      if(wald.compdata.chemo.whichmin==4){wald.compdata.chemo.whichmin="lin"} 

if(wald.compdata.logpgr.whichmin==1){wald.compdata.logpgr.whichmin="k5"} else 
  if(wald.compdata.logpgr.whichmin==2){wald.compdata.logpgr.whichmin="k4"} else
    if(wald.compdata.logpgr.whichmin==3){wald.compdata.logpgr.whichmin="k3"} else
      if(wald.compdata.logpgr.whichmin==4){wald.compdata.logpgr.whichmin="lin"} 


compdata.res=rbind(c(wald.compdata.age.k5,wald.compdata.age.k4,wald.compdata.age.k3,wald.compdata.age.lin,wald.compdata.age.whichmin,wald.compdata.age.min),
  c(wald.compdata.sized1.k5,wald.compdata.sized1.k4,wald.compdata.sized1.k3,wald.compdata.sized1.lin,wald.compdata.sized1.whichmin,wald.compdata.sized1.min),
  c(wald.compdata.sized2.k5,wald.compdata.sized2.k4,wald.compdata.sized2.k3,wald.compdata.sized2.lin,wald.compdata.sized2.whichmin,wald.compdata.sized2.min),
  c(wald.compdata.grade.k5,wald.compdata.grade.k4,wald.compdata.grade.k3,wald.compdata.grade.lin,wald.compdata.grade.whichmin,wald.compdata.grade.min),
  c(wald.compdata.enodessq.k5,wald.compdata.enodessq.k4,wald.compdata.enodessq.k3,wald.compdata.enodessq.lin,wald.compdata.enodessq.whichmin,wald.compdata.enodessq.min),
  c(wald.compdata.hormon.k5,wald.compdata.hormon.k4,wald.compdata.hormon.k3,wald.compdata.hormon.lin,wald.compdata.hormon.whichmin,wald.compdata.hormon.min),
  c(wald.compdata.chemo.k5,wald.compdata.chemo.k4,wald.compdata.chemo.k3,wald.compdata.chemo.lin,wald.compdata.chemo.whichmin,wald.compdata.chemo.min),
  c(wald.compdata.logpgr.k5,wald.compdata.logpgr.k4,wald.compdata.logpgr.k3,wald.compdata.logpgr.lin,wald.compdata.logpgr.whichmin,wald.compdata.logpgr.min))

colnames(compdata.res)=c("k5","k4","k3","lin","which.min","MIN")
rownames(compdata.res)=c("age","sized1","sized2","grade","enodessq","hormon","chemo","logpgr")

compdata.stat.res=rbind(c(stat.compdata.age.k5,stat.compdata.age.k4,stat.compdata.age.k3,stat.compdata.age.lin),
                         c(stat.compdata.sized1.k5,stat.compdata.sized1.k4,stat.compdata.sized1.k3,stat.compdata.sized1.lin),
                         c(stat.compdata.sized2.k5,stat.compdata.sized2.k4,stat.compdata.sized2.k3,stat.compdata.sized2.lin),
                         c(stat.compdata.grade.k5,stat.compdata.grade.k4,stat.compdata.grade.k3,stat.compdata.grade.lin),
                         c(stat.compdata.enodessq.k5,stat.compdata.enodessq.k4,stat.compdata.enodessq.k3,stat.compdata.enodessq.lin),
                         c(stat.compdata.hormon.k5,stat.compdata.hormon.k4,stat.compdata.hormon.k3,stat.compdata.hormon.lin),
                         c(stat.compdata.chemo.k5,stat.compdata.chemo.k4,stat.compdata.chemo.k3,stat.compdata.chemo.lin),
                         c(stat.compdata.logpgr.k5,stat.compdata.logpgr.k4,stat.compdata.logpgr.k3,stat.compdata.logpgr.lin))
colnames(compdata.stat.res)=c("k5","k4","k3","lin")
rownames(compdata.stat.res)=c("age","sized1","sized2","grade","enodessq","hormon","chemo","logpgr")

compdata.res=compdata.res[vars.remaining,]
compdata.stat.res=compdata.stat.res[vars.remaining,]

# compdata.res
# 
# compdata.stat.res

