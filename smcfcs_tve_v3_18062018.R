#' @param originaldata The original data. The data set should NOT be split on event times.
#' @param id Unique identifier for individuals in the data. There is assumed to be only 1 record per individual. 
#' @param smformula.timefixed Formula for the time-fixed part of the survival model, e.g. smformula.timefixed="Surv(t, d) ~ x1+x2"
#' @param method As in standard smcfcs
#' @param tvar Vector of expanatory variables that are to have a TVE in the imputation. The default is for all variables in smformula.timefixed to be included in tvar. 
#' @param tvarfunc Function of time that each variable in tvar will be multiplied by.
#' @param spline Indicator of whether the TVE function should be a restricted cubic spline. Default is FALSE. tvarfunc and spline cannot both be specified.
#' @param knots If spline=T, either knots or nknots must be specified (not both). Knots indicates the position of the internal knots for a restricted cubic spline.
#' @param nknots If spline=T, either knots or nknots must be specified (not both). nknots indicates the number internal knots for a restricted cubic spline. If nknots is selected, the knots will be placed at ...
#' @param predictorMatrix As in standard smcfcs
#' @param m As in standard smcfcs
#' @param numit As in standard smcfcs
#' @param rjlimit As in standard smcfcs
#' @param noisy As in standard smcfcs
#' @param errorProneMatrix As in standard smcfcs
#' @param quick If quick=T the cox models fitted as part of the imputation process will be split on a more coarse grid of time points, instead of on the event times. This will save some in settings with a large number of event.
#' @param quick.cuts A vector of cut points to be used in the data splitting, e.g. quick.cuts=seq(0,10,0.5) will split the survival dat on 21 times.


smcfcs.tve<- function(originaldata,id,smformula.timefixed,method,tvar=NULL,tvarfunc=NULL,spline=NULL,knots=NULL,nknots=NULL,predictorMatrix=NULL,m=5,numit=10,rjlimit=1000,noisy=FALSE,quick=F,quick.cuts=NULL,
                      errorProneMatrix=NULL) {
  
  stopifnot(is.data.frame(originaldata))
  if (ncol(originaldata)!=length(method)) stop("Method argument must have the same length as the number of columns in the data frame.")
  
  n <- dim(originaldata)[1]
  smformula=smformula.timefixed #RUTH
  
  #------------------
  #find column numbers of partially observed, fully observed variables, outcome, and explanatory variables collectively
  
  timeCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[2]])]
  dCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[3]])]
  outcomeCol <- c(timeCol, dCol)
  
  originaldata<-originaldata[order(originaldata[,timeCol]),] #RUTH: order data by event/censoring times - this is required for calculating the estimated survival probs
  
  d <- originaldata[,dCol]
  #explanCol<-(1:dim(originaldata)[2])[colnames(originaldata) %in% unlist(strsplit(toString(as.formula(smformula)[[3]]),", "))] #RUTH. Change from v1
  temp.explancol=sapply(c(unlist(strsplit(capture.output(print(as.formula(smformula)[[3]])), "[+]"))),FUN=function(x)stringr::str_trim(x))#RUTH. Change from v1
  explanCol<-(1:dim(originaldata)[2])[colnames(originaldata) %in% temp.explancol]#RUTH. Change from v1
  names.explanCol=names(originaldata[,explanCol])

  #------------------
  #RUTH: specification of the TVE function
  
  if(is.null(tvarfunc) & is.null(spline))stop("You must specify either tvarfunc or spline. If TVEs are not required, use standard smcfcs.")
  if(!is.null(tvarfunc) & !is.null(spline))stop("You must specify either tvarfunc or spline, not both.")
  if(!is.null(tvarfunc) & !is.null(spline))stop("You must specify either tvarfunc or spline, not both.")
  if(is.null(spline) & !is.null(knots))stop("You must specify spline =T if knots is used.")
  if(is.null(spline) & !is.null(nknots))stop("You must specify spline =T if nknots is used.")
  if(!is.null(knots) & !is.null(nknots))stop("You cannot specify both knots and nknots.")
  
  if((!is.null(tvarfunc) & is.null(tvar))|(!is.null(spline) & is.null(tvar))){
    tvar= names.explanCol
  }
  print(paste("TVE for variable ",tvar))
  
  if(!is.null(tvarfunc)) print(paste("The TVE function is as specified in tvarfunc."))
  if(!is.null(spline) & !is.null(nknots)) print(paste("The TVE function is a restricted cubic spline with ",nknots," knots."))
  if(!is.null(spline) & !is.null(knots)) print(paste("The TVE function is a restricted cubic spline with knots ",knots))
  
  #------------------
  #name the time and event indicator columns "t" and "d"
  
  colnames(originaldata)[timeCol]="t"
  colnames(originaldata)[dCol]="d"
  
  #------------------
  #cutspoints to be used later in survSplit
  cut.points.all=unique(originaldata[,timeCol][originaldata[,dCol]==1])
  if(quick==T){
    cut.points=quick.cuts
  }else{
    cut.points=cut.points.all
  }
  
  if(!is.null(tvarfunc)){
    cut.points.tvarfunc=sapply(cut.points.all,function(t)tvarfunc(t))
  }else if(!is.null(spline) & !is.null(nknots)) {
    cut.points.sp=rcspline.eval(cut.points.all,nk=nknots,inclx=T)
    knots.sp=attributes(cut.points.sp)$knots
  }else if(!is.null(spline) & !is.null(knots)) {
    nknots=length(knots)
    knots.sp=knots
    cut.points.sp=rcspline.eval(cut.points.all,knots=knots,inclx=T)
  }
  
  #------------------
  #transformed times: these are used in later calculations
  
  if(!is.null(tvarfunc)){
    timetransform.tvarfunc=sapply(originaldata[,timeCol],function(t)tvarfunc(t)) 
  }else if(!is.null(spline)) {
    timetransform.sp=rcspline.eval(originaldata[,timeCol],knots=knots.sp,inclx=T)
  }
  
  
  originaldata.split<- survSplit(Surv(t,d)~.,data = originaldata, cut = cut.points, end = "t",start = "t0", event = "d") 
  
  #originaldata.split.reduced<-originaldata.split[,names(originaldata.split)[!names(originaldata.split)%in%names.explanCol]]#this is originaldata.split with the explanatory variables with missingness removed
  #this is a bit tortured but used to speed up one of the sloopes below (see creation of imp.split)
  
  #imp.split.agg=aggregate(rep(1,dim(originaldata.split)[1])~originaldata.split$id,FUN=sum)[rank(originaldata$id),][,2]#used in creation of imp.split
  imp.split.agg=aggregate(rep(1,dim(originaldata.split)[1])~originaldata.split[,id],FUN=sum)[rank(originaldata[,id]),][,2]#RUTH: CHANGE MADE FROM V2 SO THAT THE ID VARIABLE DOESN'T HAVE TO BE NAMED "ID"
  
  if(!is.null(tvarfunc)){
    split.times.tvarfunc=sapply(originaldata.split$t,function(t)tvarfunc(t))
  }else if(!is.null(spline)) {
    split.times.t=originaldata.split$t
    split.times.t.sp=rcspline.eval(originaldata.split$t,knots=knots.sp,inclx=T)
  }
  
  #-----------------
  #create matrix of response indicators
  r <- 1*(is.na(originaldata)==0)
  
  smcovnames <- attr(terms(as.formula(smformula)), "term.labels")
  smcovcols <- (1:ncol(originaldata))[colnames(originaldata) %in% smcovnames]
  
  #-----------------
  #partial vars are those variables for which an imputation method has been specified among the available regression types
  partialVars <- which((method=="norm") | (method=="logreg") | (method=="poisson") | (method=="podds") | (method=="mlogit") | (method=="latnorm"))
  
  if (length(outcomeCol)==1) {
    if (method[outcomeCol]!="") stop("The element of the method argument corresponding to the outcome variable should be empty.")
  } else {
    if (sum(method[outcomeCol]!=c("",""))>0) stop("The elements of the method argument corresponding to the outcome variables should be empty.")
  }
  
  #-----------------
  #fully observed vars are those that are fully observed and are covariates in the substantive model
  fullObsVars <- which((colSums(r)==n) & (colnames(originaldata) %in% smcovnames))
  
  #-----------------
  #passive variables
  passiveVars <- which((method!="") & (method!="norm") & (method!="logreg") & (method!="poisson") & (method!="podds") & (method!="mlogit") & (method!="latnorm"))
  
  print(paste("Outcome variable(s):", paste(colnames(originaldata)[outcomeCol],collapse=',')))
  print(paste("Passive variables:", paste(colnames(originaldata)[passiveVars],collapse=',')))
  print(paste("Partially obs. variables:", paste(colnames(originaldata)[partialVars],collapse=',')))
  print(paste("Fully obs. substantive model variables:", paste(colnames(originaldata)[fullObsVars],collapse=',')))
  
  imputations <- list()
  for (imp in 1:m) {
    imputations[[imp]] <- originaldata
  }
  
  for (imp in 1:m) {
    
    print(paste("Imputation ",imp))
    
    #-----------------
    #initial imputation of each partially observed variable based on observed values
    for (var in 1:length(partialVars)) {
      targetCol <- partialVars[var]
      if (method[targetCol]=="latnorm") {
        #initialize latent predictors with single error prone measurement
        errorProneCols <- which(errorProneMatrix[targetCol,]==1)
        imputations[[imp]][,targetCol] <- apply(imputations[[imp]][,errorProneCols], 1, firstnonna)
      } else {
        imputations[[imp]][r[,targetCol]==0,targetCol] <- sample(imputations[[imp]][r[,targetCol]==1,targetCol], size=sum(r[,targetCol]==0), replace=TRUE)
      }
    }
    
    #-----------------
    #initial imputations of missing outcomes, if present (using improper imputation)
    
    for (cyclenum in 1:numit) {
      
      #if (noisy==TRUE) {
      print(paste("Iteration ",cyclenum))
      #}
      #update passive variable(s)
      imputations[[imp]] <- updatePassiveVars(imputations[[imp]], method, passiveVars)
      
      for (var in 1:length(partialVars)) {
        targetCol <- partialVars[var]
        if (is.null(predictorMatrix)) {
          predictorCols <- c(partialVars[! partialVars %in% targetCol], fullObsVars)
        }else {
          predictorCols <- which(predictorMatrix[targetCol,]==1)
          #ensure that user has not included outcome variable(s) here
          predictorCols <- predictorCols[! predictorCols %in% outcomeCol]
        }
        if ((imp==1) & (cyclenum==1)) {
          if (method[targetCol]=="latnorm") {
            print(paste("Imputing: ",colnames(imputations[[imp]])[targetCol]," using ",paste(colnames(imputations[[imp]])[c(predictorCols,which(errorProneMatrix[targetCol,]==1))],collapse=',')," plus outcome",collapse=','))
          } else {
            print(paste("Imputing: ",colnames(imputations[[imp]])[targetCol]," using ",paste(colnames(imputations[[imp]])[predictorCols],collapse=',')," plus outcome",collapse=','))
          }
        }
        if (length(predictorCols)>0) {
          xmodformula <- as.formula(paste(colnames(imputations[[imp]])[targetCol], "~", paste(colnames(imputations[[imp]])[predictorCols], collapse="+"),sep=""))
        }else {
          xmodformula <- as.formula(paste(colnames(imputations[[imp]])[targetCol], "~1",sep=""))
        }
        if ((method[targetCol]=="norm") | (method[targetCol]=="latnorm")) {
          #estimate parameters of covariate model
          xmod <- lm(xmodformula, data=imputations[[imp]])
          #take draw from posterior of covariate model parameters
          beta <- xmod$coef
          sigmasq <- summary(xmod)$sigma^2
          newsigmasq <- (sigmasq*xmod$df) / rchisq(1,xmod$df)
          covariance <- (newsigmasq/sigmasq)*vcov(xmod)
          newbeta = beta + MASS::mvrnorm(1, mu=rep(0,ncol(covariance)), Sigma=covariance)
          #calculate fitted values
          xfitted <- model.matrix(xmod) %*% newbeta
        }else if (method[targetCol]=="logreg") {
          xmod <- glm(xmodformula, family="binomial",data=imputations[[imp]])
          newbeta = modPostDraw(xmod)
          xfitted <- expit(model.matrix(xmod) %*% newbeta)
        } else if (method[targetCol]=="poisson") {
          xmod <- glm(xmodformula, family="poisson", data=imputations[[imp]])
          newbeta = modPostDraw(xmod)
          xfitted <- exp(model.matrix(xmod) %*% newbeta)
        } else if (method[targetCol]=="podds") {
          xmod <- VGAM::vglm(xmodformula, VGAM::propodds, data=imputations[[imp]])
          newbeta <- coefficients(xmod) + MASS::mvrnorm(1, mu=rep(0,ncol(vcov(xmod))), Sigma=vcov(xmod))
          linpreds <- matrix(model.matrix(xmod) %*% newbeta, byrow=TRUE, ncol=(nlevels(imputations[[imp]][,targetCol])-1))
          cumprobs <- cbind(1/(1+exp(linpreds)), rep(1,nrow(linpreds)))
          xfitted <- cbind(cumprobs[,1] ,cumprobs[,2:ncol(cumprobs)] - cumprobs[,1:(ncol(cumprobs)-1)])
        } else if (method[targetCol]=="mlogit") {
          xmod <- VGAM::vglm(xmodformula, VGAM::multinomial(refLevel=1), data=imputations[[imp]])
          newbeta <- coef(xmod) + MASS::mvrnorm(1, mu=rep(0,ncol(vcov(xmod))), Sigma=vcov(xmod))
          linpreds <- matrix(model.matrix(xmod) %*% newbeta, byrow=TRUE, ncol=(nlevels(imputations[[imp]][,targetCol])-1))
          denom <- 1+rowSums(exp(linpreds))
          xfitted <-cbind(1/denom, exp(linpreds) / denom)
        }
        if (noisy==TRUE) {
          print(summary(xmod))
        }
        
        #if latent normal, estimate error variance and calculate required conditional expectation and variance
        if (method[targetCol]=="latnorm") {
          errorProneCols <- which(errorProneMatrix[targetCol,]==1)
          wmean <- rowMeans(imputations[[imp]][,errorProneCols], na.rm=TRUE)
          n_i <- apply(imputations[[imp]][,errorProneCols], 1, sumna)
          sum_ni <- sum(n_i)
          #estimate error variance
          if (cyclenum==1) {
            xmat <- matrix(wmean, nrow=nrow(imputations[[imp]]), ncol=length(errorProneCols))
            uVec <- c(as.matrix(imputations[[imp]][,errorProneCols] - xmat))
            sigmausq <- sum(uVec^2, na.rm=TRUE) / (sum_ni - n)
          }
          else {
            xmat <- matrix(imputations[[imp]][,targetCol], nrow=nrow(imputations[[imp]]), ncol=length(errorProneCols))
            uVec <- c(as.matrix(imputations[[imp]][,errorProneCols] - xmat))
            sigmausq <- sum(uVec^2, na.rm=TRUE) / sum_ni
          }
          #take draw from posterior of error variance
          sigmausq <- sigmausq*sum_ni/rchisq(1,sum_ni)
          #calculate conditional mean and variance
          lambda <- newsigmasq/(newsigmasq+sigmausq/n_i)
          xfitted <- xfitted + lambda * (wmean - xfitted)
          newsigmasq <- newsigmasq*(1-lambda)
        }
        
        #------------------------
        #estimate parameters of substantive model
        #RUTH:SPLIT THE DATA ON ALL EVENT AND CENSORING TIMES
        #There are different ways of doing this (see options commented out below. The one chosen here is fastest (at least for not too many explanatory variables, as in my simulations.
        #imp.split<- survSplit(Surv(t,d)~.,data = imputations[[imp]], cut = cut.points, end = "t",start = "t0", event = "d") 
        #imp.split=merge(originaldata.split.reduced,imputations[[imp]][,c(id,names.explanCol)],by="id",sort=F)
        imp.split=originaldata.split
        for(n.explan in 1:length(names.explanCol)){
          imp.split[,names.explanCol[n.explan]]=rep(imputations[[imp]][,names.explanCol[n.explan]],imp.split.agg)
        }
        
        #RUTH: update the smformula to include time-varying terms and generate covariate-by-time interactions
        
        if(!is.null(tvarfunc)){
          names.explanCol.new=paste0(names.explanCol,".t")
          extra.terms=paste0(names.explanCol.new,collapse = "+")
          #smformula2<-paste0(capture.output(as.formula(smformula)),"+",extra.terms)
          
          smformula2<-as.formula(paste0(smformula,"+",extra.terms))
          for(n.explan in 1:length(names.explanCol)){
            eval(parse(text=paste0("imp.split$",names.explanCol[n.explan],".t<-imp.split[,\"",
                                   names.explanCol[n.explan],"\"]*tvarfunc(imp.split$t)")))}
        }else if(!is.null(spline)){
          names.explanCol.new=paste0(names.explanCol,".t")
          names.explanCol.new.mat=outer(names.explanCol.new,1:(nknots-1),paste,sep="")
          names.explanCol.new=as.vector(t(names.explanCol.new.mat))
          extra.terms=paste0(names.explanCol.new,collapse = "+")
          #smformula2<-paste(capture.output(as.formula(smformula)),"+",extra.terms)
          smformula2<-as.formula(paste0(smformula,"+",extra.terms))
          
          temp=lapply(1:length(names.explanCol), function(k)imp.split[,names.explanCol[k]]*as.matrix(rcspline.eval(imp.split$t,knots=knots.sp,inclx=T)))
          imp.split[,names.explanCol.new]=matrix(unlist(temp),ncol=2*(nknots-1))
        }
        
        #we now alter smformula2 to add t0 into the Surv formula, i.e. Surv(t, d) becomes Surv(t0,t, d)
        smformula3<-smformula2 #this is used later
        if(!is.null(tvarfunc)){
          orig.terms=paste0(names.explanCol,collapse = "+")
          extra.terms=paste0(names.explanCol,".t",collapse = "+")
          smformula2<-paste0("Surv (t0, t , d )~",orig.terms,"+",extra.terms)
        }else if(!is.null(spline)){
          orig.terms=paste0(names.explanCol,collapse = "+")
          extra.terms=paste0(names.explanCol.new,collapse = "+")
          smformula2<-paste0("Surv (t0, t , d )~",orig.terms,"+",extra.terms)
        }
        #---
        #fit the substantive model with time-varying effects
        ymod <- survival::coxph(as.formula(smformula2), imp.split)
        outcomeModBeta <- modPostDraw(ymod)
        ymod$coefficients <- outcomeModBeta
        
        #---
        #RUTH: obtain baseline cumulative hazard
        #this doesn't depend on the imputed values, so only needs to be calculated once
        
        #method 1
        # eval(parse(text=paste0("my.newdata=data.frame(t=cut.points,t0=c(0,cut.points[-length(cut.points)]),d=1,id=1,",
        #                       paste0(names.explanCol,"=0",collapse=","),",",paste0(names.explanCol.new,"=0",collapse=","),")")))
        # basecumhaz=unique(survfit(ymod,newdata=my.newdata,id=id)$cumhaz) 
        #method 2
        basecumhaz=unique(basehaz(ymod,centered=F)$hazard)
        
        if(originaldata$d[1]==0){basehaz=diff(basecumhaz,lag=1)} else{ #this is the increments in the cumulative baseline hazard
          basehaz=c(basecumhaz[1],diff(basecumhaz,lag=1))}
        #RUTH: now calculate exp(beta0+beta1*t+splines funcs) for each event time. Separately for each explanatory variable with a time-varying effect These only need to be calculated once. 
        
        if(!is.null(tvarfunc)){
          tfunc.matrix=as.matrix(cbind(rep(1,length(cut.points.all)),cut.points.tvarfunc))
          for(n.explan in 1:length(names.explanCol)){
            outcomeModBeta.temp=as.matrix(outcomeModBeta[c(names.explanCol[n.explan],names.explanCol.new[n.explan])])
            eval(parse(text=paste0("exp.t.",names.explanCol[n.explan],"=c(exp(t(outcomeModBeta.temp)%*%t(tfunc.matrix)))")))
            eval(parse(text=paste0("outer.exp.t.",names.explanCol[n.explan],
                                   "=outer(exp.t.",names.explanCol[n.explan],",imputations[[imp]][,",explanCol[n.explan],"],FUN=\"^\")")))
          }
        }else if(!is.null(spline)){
          tfunc.matrix=as.matrix(cbind(rep(1,length(cut.points.all)),cut.points.sp))
          for(n.explan in 1:length(names.explanCol)){
            outcomeModBeta.temp=as.matrix(outcomeModBeta[c(names.explanCol[n.explan],names.explanCol.new.mat[n.explan,])])
            eval(parse(text=paste0("exp.t.",names.explanCol[n.explan],"=c(exp(t(outcomeModBeta.temp)%*%t(tfunc.matrix)))")))
            eval(parse(text=paste0("outer.exp.t.",names.explanCol[n.explan],
                                   "=outer(exp.t.",names.explanCol[n.explan],",imputations[[imp]][,",explanCol[n.explan],"],FUN=\"^\")")))
          }
        }
        
        #RUTH: for each individual calculate the sum of basehaz*exp.t^x for all times up to the event/censoring time. This is their cumulative hazard.
        eval(parse(text=paste0("cumhaz.elements=basehaz*",paste0("outer.exp.t.",names.explanCol,collapse="*"))))
        
        #cumhaz.elements are the terms at each event time, which are to be summed up to the event or censoring time of the individual
        list.times=originaldata[,timeCol]
        cumhaz=sapply(1:n,function(x){sum(cumhaz.elements[which(cut.points<=list.times[x]),x])})
        survest=exp(-cumhaz)#this is the estimated survival function
        #now assign the cumulative baseline hazard to each individual
        basecumhaz.all=rep(NA,n)
        if(originaldata$d[1]==0){
          basecumhaz.all[which(imputations[[imp]][,timeCol]%in%cut.points.all & duplicated(imputations[[imp]][,timeCol])==F)]=  basecumhaz[-1]
        }else {
          basecumhaz.all[which(imputations[[imp]][,timeCol]%in%cut.points.all & duplicated(imputations[[imp]][,timeCol])==F)]=  basecumhaz
        } #the 'duplicated bit here should deal with ties
        basecumhaz.all=zoo::na.locf(basecumhaz.all,na.rm=F)
        basecumhaz.all=ifelse(is.na(basecumhaz.all)==1,0,basecumhaz.all)
        
        basehaz.all=rep(NA,n)
        basehaz.all[which(imputations[[imp]][,timeCol]%in%cut.points.all & duplicated(imputations[[imp]][,timeCol])==F)]=  basehaz #the 'duplicated bit here should deal with ties
        basehaz.all=zoo::na.locf(basehaz.all,na.rm=F)
        basehaz.all=ifelse(is.na(basehaz.all)==1,0,basehaz.all)
        #end of baseline hazard calcs
        #---
        
        if (noisy==TRUE) {
          print(summary(ymod))
        }
        
        
        
        if ((imp==1) & (cyclenum==1) & (var==1)) {
          
          smCoefIter <- array(0, dim=c(m, length(outcomeModBeta), numit))
          
        }
        
        if (var==length(partialVars)) {
          #then we have reached end of a cycle
          
          smCoefIter[imp,,cyclenum] <- outcomeModBeta
          
        }
        
        #impute x, either directly where possibly, or using rejection sampling otherwise
        imputationNeeded <- (1:n)[r[,targetCol]==0]
        
        #----------------------------------------------------
        #RUTH: THIS IS THE START OF THE DIRECT SAMPLING PART
        #----------------------------------------------------
        
        if ((method[targetCol]=="logreg") | (method[targetCol]=="podds") | (method[targetCol]=="mlogit")) {  
          
          #directly sample
          if (method[targetCol]=="logreg") {
            numberOutcomes <- 2
            fittedMean <- cbind(1-xfitted, xfitted)
          } else {
            numberOutcomes <- nlevels(imputations[[imp]][,targetCol])
            fittedMean <- xfitted
          }
          
          outcomeDensCovDens = array(dim=c(length(imputationNeeded),numberOutcomes),0)
          
          for (xMisVal in 1:numberOutcomes) {
            if (method[targetCol]=="logreg") {
              if (is.factor(imputations[[imp]][,targetCol])==TRUE) {
                valToImpute <- levels(imputations[[imp]][,targetCol])[xMisVal]
              }else {
                valToImpute <- xMisVal-1
              }
            } else {
              valToImpute <- levels(imputations[[imp]][,targetCol])[xMisVal]
            }
            imputations[[imp]][imputationNeeded,targetCol] <- valToImpute
            
            #update passive variables
            imputations[[imp]] <- updatePassiveVars(imputations[[imp]], method, passiveVars)
            
            #RUTH: generate covariate-by-time interactions in imputations[[imp]]
            if(!is.null(tvarfunc)){
              for(n.explan in 1:length(names.explanCol)){
                eval(parse(text=paste0("imputations[[imp]]$",names.explanCol[n.explan],
                                       ".t<-imputations[[imp]][,names.explanCol[n.explan]]*
                                       tvarfunc(imputations[[imp]][,timeCol])")))}
              }else if(!is.null(spline)){
                temp=lapply(1:length(names.explanCol), function(k)imputations[[imp]][,names.explanCol[k]]*as.matrix(timetransform.sp))
                imputations[[imp]][,names.explanCol.new]=matrix(unlist(temp),ncol=2*(nknots-1))
            }
            
            outmodxb <-  model.matrix(as.formula(smformula3),imputations[[imp]])
            outmodxb <- outmodxb[,2:dim(outmodxb)[2]] %*% outcomeModBeta
            
            #------------
            #RUTH: estimate survivor probs
            #------------
            
            if(!is.null(tvarfunc)){
              tfunc.matrix=as.matrix(cbind(rep(1,length(cut.points.all)),cut.points.tvarfunc))
              for(n.explan in 1:length(names.explanCol)){
                outcomeModBeta.temp=as.matrix(outcomeModBeta[c(names.explanCol[n.explan],names.explanCol.new[n.explan])])
                eval(parse(text=paste0("exp.t.",names.explanCol[n.explan],"=c(exp(t(outcomeModBeta.temp)%*%t(tfunc.matrix)))")))
                eval(parse(text=paste0("outer.exp.t.",names.explanCol[n.explan],
                                       "=outer(exp.t.",names.explanCol[n.explan],",imputations[[imp]][,",explanCol[n.explan],"],FUN=\"^\")")))
              }
            }else if(!is.null(spline)){
              tfunc.matrix=as.matrix(cbind(rep(1,length(cut.points.all)),cut.points.sp))
              for(n.explan in 1:length(names.explanCol)){
                outcomeModBeta.temp=as.matrix(outcomeModBeta[c(names.explanCol[n.explan],names.explanCol.new.mat[n.explan,])])
                eval(parse(text=paste0("exp.t.",names.explanCol[n.explan],"=c(exp(t(outcomeModBeta.temp)%*%t(tfunc.matrix)))")))
                eval(parse(text=paste0("outer.exp.t.",names.explanCol[n.explan],
                                       "=outer(exp.t.",names.explanCol[n.explan],",imputations[[imp]][,",explanCol[n.explan],"],FUN=\"^\")")))
              }
            }
            
            #RUTH: for each individual calculate the sum of basehaz*exp.t^x for all times up to the event/censoring time. This is their cumulative hazard.
            eval(parse(text=paste0("cumhaz.elements=basehaz*",paste0("outer.exp.t.",names.explanCol,collapse="*"))))
            
            #cumhaz.elements are the terms at each event time, which are to be summed up to the event or censoring time of the individual
            list.times=originaldata[,timeCol]
            cumhaz=sapply(1:n,function(x){sum(cumhaz.elements[which(cut.points.all<=list.times[x]),x])})
            survest=exp(-cumhaz)#this is the estimated survival function
            
            #--------------
            
            #outcomeDens <- exp(-H0[imputationNeeded] * exp(outmodxb[imputationNeeded]))* (exp(outmodxb[imputationNeeded])^d[imputationNeeded])
            outcomeDens <- survest[imputationNeeded]* (exp(outmodxb[imputationNeeded]))^d[imputationNeeded]#RUTH: added - not sure about this!
            
            outcomeDensCovDens[,xMisVal] <- outcomeDens * fittedMean[imputationNeeded,xMisVal]
          }#end of xMisVal loop
          directImpProbs = outcomeDensCovDens / rowSums(outcomeDensCovDens)
          
          if (method[targetCol]=="logreg") {
            directImpProbs = directImpProbs[,2]
            if (is.factor(imputations[[imp]][,targetCol])==TRUE) {
              imputations[[imp]][imputationNeeded,targetCol] <- levels(imputations[[imp]][,targetCol])[1]
              imputations[[imp]][imputationNeeded,targetCol][rbinom(length(imputationNeeded),1,directImpProbs)==1] <- levels(imputations[[imp]][,targetCol])[2]
            }
            else {
              imputations[[imp]][imputationNeeded,targetCol] <- rbinom(length(imputationNeeded),1,directImpProbs)
            }
          }
          else {
            imputations[[imp]][imputationNeeded,targetCol] <- levels(imputations[[imp]][,targetCol])[apply(directImpProbs, 1, catdraw)]
          }
          
          #update passive variables
          imputations[[imp]] <- updatePassiveVars(imputations[[imp]], method, passiveVars)
        }
        #RUTH: THIS IS THE END OF THE DIRECT SAMPLING PART
        #----------------------------------------------------
        #----------------------------------------------------
        #RUTH: THIS IS THE START OF THE REJECTION SAMPLING PART
        #----------------------------------------------------
        else {
          #use rejection sampling
          #first draw for all subjects who need imputing, using a small number of attempts
          firstTryLimit <- 25
          j <- 1
          #print(firstTryLimit)
          #print(length(imputationNeeded))
          
          while ((length(imputationNeeded)>0) & (j<firstTryLimit)) {
            #print(j)
            #sample from covariate model
            if ((method[targetCol]=="norm") | (method[targetCol]=="latnorm")) {
              imputations[[imp]][imputationNeeded,targetCol] <- rnorm(length(imputationNeeded),xfitted[imputationNeeded],newsigmasq^0.5)
            } else if (method[targetCol]=="poisson") {
              imputations[[imp]][imputationNeeded,targetCol] <- rpois(length(imputationNeeded),xfitted[imputationNeeded])
            }
            
            #update passive variables
            imputations[[imp]] <- updatePassiveVars(imputations[[imp]], method, passiveVars)
            
            #accept/reject
            uDraw <- runif(length(imputationNeeded))
            
            # imputations[[imp]]$x1.miss.t=imputations[[imp]]$x1.miss*tvarfunc(imputations[[imp]][,timeCol])
            # imputations[[imp]]$x2.miss.t=imputations[[imp]]$x2.miss*tvarfunc(imputations[[imp]][,timeCol])
            
            # #RUTH: generate covariate-by-time interactions in imputations[[imp]]
            if(!is.null(tvarfunc)){
              for(n.explan in 1:length(names.explanCol)){
                eval(parse(text=paste0("imputations[[imp]]$",names.explanCol[n.explan],
                                       ".t<-imputations[[imp]][,names.explanCol[n.explan]]*
                                       tvarfunc(imputations[[imp]][,timeCol])")))}
              }else if(!is.null(spline)){
                temp=lapply(1:length(names.explanCol), function(k)imputations[[imp]][,names.explanCol[k]]*as.matrix(timetransform.sp))
                imputations[[imp]][,names.explanCol.new]=matrix(unlist(temp),ncol=2*(nknots-1))
            }
            
            outmodxb<-model.matrix(as.formula(smformula3),imputations[[imp]])
            
            outmodxb <- outmodxb[,2:dim(outmodxb)[2]] %*% outcomeModBeta
            
            #RUTH: for each individual calculate the sum of basehaz*exp.t^x for all times up to the event/censoring time. This is their cumulative hazard.
            if(!is.null(tvarfunc)){
              for(n.explan in 1:length(names.explanCol)){
                outcomeModBeta.temp=as.matrix(outcomeModBeta[c(names.explanCol[n.explan],names.explanCol.new[n.explan])])
                eval(parse(text=paste0("exp.t.",names.explanCol[n.explan],"=c(exp(t(outcomeModBeta.temp)%*%t(tfunc.matrix)))")))
                eval(parse(text=paste0("outer.exp.t.",names.explanCol[n.explan],
                                       "=outer(exp.t.",names.explanCol[n.explan],",imputations[[imp]][,",explanCol[n.explan],"],FUN=\"^\")")))
              }
            }else if(!is.null(spline)){
              for(n.explan in 1:length(names.explanCol)){
                outcomeModBeta.temp=as.matrix(outcomeModBeta[c(names.explanCol[n.explan],names.explanCol.new.mat[n.explan,])])
                eval(parse(text=paste0("exp.t.",names.explanCol[n.explan],"=c(exp(t(outcomeModBeta.temp)%*%t(tfunc.matrix)))")))
                eval(parse(text=paste0("outer.exp.t.",names.explanCol[n.explan],
                                       "=outer(exp.t.",names.explanCol[n.explan],",imputations[[imp]][,",explanCol[n.explan],"],FUN=\"^\")")))
              }
            }
            eval(parse(text=paste0("cumhaz.elements=basehaz*",paste0("outer.exp.t.",names.explanCol,collapse="*"))))
            #cumhaz.elements are the terms at each event time, which are to be summed up to the event or censoring time of the individual
            cumhaz=sapply(1:n,function(x){sum(cumhaz.elements[which(cut.points.all<=list.times[x]),x])})
            survest=exp(-cumhaz)#this is the estimated survival function
            #----------------
            
            #s_t = exp(-H0[imputationNeeded]* exp(outmodxb[imputationNeeded]))
            #prob = exp(1 + outmodxb[imputationNeeded] - (H0[imputationNeeded]* exp(outmodxb[imputationNeeded])) ) * H0[imputationNeeded]
            #prob = d[imputationNeeded]*prob + (1-d[imputationNeeded])*s_t
            prob.d0<-survest[imputationNeeded] #RUTH
            prob.d1<-basehaz.all[imputationNeeded]*exp(1+outmodxb[imputationNeeded])*survest[imputationNeeded] #RUTH
            prob = d[imputationNeeded]*prob.d1 + (1-d[imputationNeeded])*prob.d0 #RUTH
            reject = 1*(uDraw > prob )
            
            imputationNeeded <- imputationNeeded[reject==1]
            
            j <- j+1
          }
          #print(length(imputationNeeded))
          #now, for those remaining, who must have low acceptance probabilities, sample by subject
          for (i in imputationNeeded) {
            #print("Imputation needed for:")
            #print(i)
            
            tempData <- imputations[[imp]][i,]
            tempData <- tempData[rep(1,rjlimit),]
            if (method[targetCol]=="norm") {
              tempData[,targetCol] <- rnorm(rjlimit,xfitted[i],newsigmasq^0.5)
            } else if (method[targetCol]=="logreg") {
              tempData[,targetCol] <- rbinom(rjlimit,size=1,xfitted[i])
            } else if (method[targetCol]=="poisson") {
              tempData[,targetCol] <- rpois(rjlimit,xfitted[i])
            } else if (method[targetCol]=="latnorm") {
              tempData[,targetCol] <- rnorm(rjlimit,xfitted[i],newsigmasq[i]^0.5)
            }
            
            #passively impute
            tempData <- updatePassiveVars(tempData, method, passiveVars)
            
            #accept reject
            uDraw <- runif(rjlimit)
            
            #RUTH: generate covariate-by-time interactions in tempData
            
            if(!is.null(tvarfunc)){
              for(n.explan in 1:length(names.explanCol)){
                eval(parse(text=paste0("tempData$",names.explanCol[n.explan],
                                       ".t<-tempData[,names.explanCol[n.explan]]*tvarfunc(tempData[,timeCol])")))}
            }else if(!is.null(spline) & !is.null(nknots)){
              temp=lapply(1:length(names.explanCol), function(k)tempData[,names.explanCol[k]]*as.matrix(rcspline.eval(tempData[,timeCol],knots=knots.sp,inclx=T)))
              tempData[,names.explanCol.new]=matrix(unlist(temp),ncol=2*(nknots-1))
            }
            
            outmodxb <-  model.matrix(as.formula(smformula3),tempData)
            
            outmodxb <- outmodxb[,2:dim(outmodxb)[2]] %*% outcomeModBeta
            
            #RUTH: for each individual calculate the sum of basehaz*exp.t^x for all times up to the event/censoring time. This is their cumulative hazard.
            
            if(!is.null(tvarfunc)){
              for(n.explan in 1:length(names.explanCol)){
                outcomeModBeta.temp=as.matrix(outcomeModBeta[c(names.explanCol[n.explan],names.explanCol.new[n.explan])])
                eval(parse(text=paste0("exp.t.",names.explanCol[n.explan],"=c(exp(t(outcomeModBeta.temp)%*%t(tfunc.matrix)))")))
                eval(parse(text=paste0("outer.exp.t.",names.explanCol[n.explan],
                                       "=outer(exp.t.",names.explanCol[n.explan],",tempData[,",explanCol[n.explan],"],FUN=\"^\")")))
              }
            }else if(!is.null(spline)){
              for(n.explan in 1:length(names.explanCol)){
                outcomeModBeta.temp=as.matrix(outcomeModBeta[c(names.explanCol[n.explan],names.explanCol.new.mat[n.explan,])])
                eval(parse(text=paste0("exp.t.",names.explanCol[n.explan],"=c(exp(t(outcomeModBeta.temp)%*%t(tfunc.matrix)))")))
                eval(parse(text=paste0("outer.exp.t.",names.explanCol[n.explan],
                                       "=outer(exp.t.",names.explanCol[n.explan],",tempData[,",explanCol[n.explan],"],FUN=\"^\")")))
              }
            }
            eval(parse(text=paste0("cumhaz.elements=basehaz*",paste0("outer.exp.t.",names.explanCol,collapse="*"))))
            #cumhaz.elements are the terms at each event time, which are to be summed up to the event or censoring time of the individual
            cumhaz=sapply(1:rjlimit,function(x){sum(cumhaz.elements[which(cut.points.all<=tempData[,timeCol][x]),x])})
            survest.rej=exp(-cumhaz)#this is the estimated survival function
            
            prob.d0<-survest.rej #RUTH
            prob.d1<-basehaz.all[i]*exp(1+outmodxb)*survest.rej #RUTH
            prob = d[i]*prob.d1 + (1-d[i])*prob.d0 #RUTH
            reject = 1*(uDraw > prob )
            
            if (sum(reject)<rjlimit) {
              imputations[[imp]][i,targetCol] <- tempData[reject==0,targetCol][1]
            } else {
              #print("Rejection sampling has failed for one record. You may want to increase the rejecton sampling limit.")
            }
          }
          #update passive variables
          imputations[[imp]] <- updatePassiveVars(imputations[[imp]], method, passiveVars)
        }
        
        #RUTH: THIS IS THE END OF THE REJECTION SAMPLING PART
        #-----------------
        
        
        } #end of loop over variables with missing data
      
      }#end of cyclenum loop
    
  }#end of imputations loop
  
  
  list(impDatasets=imputations, smCoefIter=smCoefIter)
  
}#end of function 

updatePassiveVars <- function(data, method, passivecols) {
  for (i in passivecols) {
    data[,i] <- with(data, eval(parse(text=method[i])))
  }
  data
}

expit <- function(x) {
  exp(x)/(1+exp(x))
}

sumna <- function(x) {
  sum(is.na(x)==FALSE)
}

#returns first non missing entry of x
firstnonna <- function(x) {
  x[is.na(x)==FALSE][1]
}

catdraw <- function(prob) {
  (1:length(prob))[rmultinom(1,size=1,prob=prob)==1]
}

modPostDraw <- function(modobj) {
  beta <- modobj$coef
  varcov <- vcov(modobj)
  beta + MASS::mvrnorm(1, mu=rep(0,ncol(varcov)), Sigma=varcov)
}