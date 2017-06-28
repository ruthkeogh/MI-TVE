smcfcs.tve <- function(originaldata,smformula.timefixed,method,predictorMatrix=NULL,m=5,numit=10,rjlimit=1000,noisy=FALSE,
                            errorProneMatrix=NULL) {
  
  stopifnot(is.data.frame(originaldata))
  if (ncol(originaldata)!=length(method)) stop("Method argument must have the same length as the number of columns in the data frame.")
  
  n <- dim(originaldata)[1]
  smformula=smformula.timefixed
  
  #------------------
  #find column numbers of partially observed, fully observed variables, outcome, and explanatory variables collectively
  
  timeCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[2]])]
  dCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[3]])]
  outcomeCol <- c(timeCol, dCol)
  
  originaldata<-originaldata[order(originaldata[,timeCol]),] #order data by event/censoring times - this is required for calculating the estimated survival probs
  
  d <- originaldata[,dCol]
  explanCol<-(1:dim(originaldata)[2])[colnames(originaldata) %in% unlist(strsplit(toString(as.formula(smformula)[[3]]),", "))] 
  names.explanCol=names(originaldata[,explanCol])
  
  #------------------
  #name the time and event indicator columns "t" and "d"
  
  colnames(originaldata)[timeCol]="t"
  colnames(originaldata)[dCol]="d"
  
  #------------------
  #cutspoints to be used later in survSplit
  cut.points=unique(originaldata[,timeCol][originaldata[,dCol]==1])
  cut.points.sp5.1=sapply(cut.points,function(t)v(t,knots.5[4],knots.5[5],knots.5[1]))
  cut.points.sp5.2=sapply(cut.points,function(t)v(t,knots.5[4],knots.5[5],knots.5[2]))
  cut.points.sp5.3=sapply(cut.points,function(t)v(t,knots.5[4],knots.5[5],knots.5[3]))
  
  #------------------
  #TRANSFORMED TIMES TO BE USED LATER
  
  timetransform.sp5.1=sapply(originaldata[,timeCol],function(t)v(t,knots.5[4],knots.5[5],knots.5[1]))
  timetransform.sp5.2=sapply(originaldata[,timeCol],function(t)v(t,knots.5[4],knots.5[5],knots.5[2]))
  timetransform.sp5.3=sapply(originaldata[,timeCol],function(t)v(t,knots.5[4],knots.5[5],knots.5[3]))
  
  originaldata.split<- survSplit(Surv(t,d)~.,data = originaldata, cut = cut.points, end = "t",start = "t0", event = "d") 
  split.times.t=originaldata.split$t
  split.times.t.sp5.1=sapply(originaldata.split$t,function(t)v(t,knots.5[4],knots.5[5],knots.5[1]))
  split.times.t.sp5.2=sapply(originaldata.split$t,function(t)v(t,knots.5[4],knots.5[5],knots.5[2]))
  split.times.t.sp5.3=sapply(originaldata.split$t,function(t)v(t,knots.5[4],knots.5[5],knots.5[3]))

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
        
        #SPLIT THE DATA ON ALL EVENT AND CENSORING TIMES
        imp.split<- survSplit(Surv(t,d)~.,data = imputations[[imp]], cut = cut.points, end = "t",start = "t0", event = "d") 
        
        #update the smformula to include time-varying terms and generate covariate-by-time interactions
          smformula2<-"Surv (t , d )~ x1 + x2 + x1t + x1.sp5.1 + x1.sp5.2 + x1.sp5.3 + x2t + x2.sp5.1 + x2.sp5.2 + x2.sp5.3" 
          imp.split$x1t<-imp.split[,names.explanCol[1]]*imp.split$t
          imp.split$x2t<-imp.split[,names.explanCol[2]]*imp.split$t

          imp.split$x1.sp5.1<-imp.split[,names.explanCol[1]]*split.times.t.sp5.1
          imp.split$x1.sp5.2<-imp.split[,names.explanCol[1]]*split.times.t.sp5.2
          imp.split$x1.sp5.3<-imp.split[,names.explanCol[1]]*split.times.t.sp5.3
          
          imp.split$x2.sp5.1<-imp.split[,names.explanCol[2]]*split.times.t.sp5.1
          imp.split$x2.sp5.2<-imp.split[,names.explanCol[2]]*split.times.t.sp5.2
          imp.split$x2.sp5.3<-imp.split[,names.explanCol[2]]*split.times.t.sp5.3
        
        smformula3<-smformula2 #this is used later
        smformula2<-"Surv (t0, t , d )~ x1 + x2 + x1t + x1.sp5.1 + x1.sp5.2 + x1.sp5.3 + x2t + x2.sp5.1 + x2.sp5.2 + x2.sp5.3" 
        
        #fit the substantive model with time-varying effects
        ymod <- survival::coxph(as.formula(smformula2), imp.split)
        outcomeModBeta <- modPostDraw(ymod)
        ymod$coefficients <- outcomeModBeta

        #obtain baseline cumulative hazard
        eval(parse(text=paste0("my.newdata=data.frame(t=cut.points,t0=c(0,cut.points[-length(cut.points)]),d=1,",
                               names.explanCol[1],"=0,",names.explanCol[2],"=0,x1t=0,x2t=0,x1.sp5.1=0,x1.sp5.2=0,x1.sp5.3=0,x2.sp5.1=0,x2.sp5.2=0,x2.sp5.3=0,id=1)")))
        basecumhaz=unique(survfit(ymod,newdata=my.newdata,id=id)$cumhaz) #this doesn't depend on the imputed values, so only needs to be calculated once
        if(originaldata$d[1]==0){basehaz=diff(basecumhaz,lag=1)} else{ #this is the increments in the cum baseline hazard
        basehaz=c(basecumhaz[1],diff(basecumhaz,lag=1))}
        #now calculate exp(beta0+beta1*t+splines funcs) for each event time. Separately for x1 and x2. These only need to be calculated once. 
        exp.t.x1=exp(outcomeModBeta[1]+outcomeModBeta[3]*cut.points+outcomeModBeta[4]*cut.points.sp5.1+outcomeModBeta[5]*cut.points.sp5.2+outcomeModBeta[6]*cut.points.sp5.3)
        exp.t.x2=exp(outcomeModBeta[2]+outcomeModBeta[7]*cut.points+outcomeModBeta[8]*cut.points.sp5.1+outcomeModBeta[9]*cut.points.sp5.2+outcomeModBeta[10]*cut.points.sp5.3)
        #for each individual calculate the sum of basehaz*exp.t^x for all times up to the event/censoring time. This is their cumulative hazard.
        cumhaz.elements=basehaz*outer(exp.t.x1,imputations[[imp]][,explanCol[1]],FUN="^")*
          outer(exp.t.x2,imputations[[imp]][,explanCol[2]],FUN="^") 
        #cumhaz.elements are the terms at each event time, which are to be summed up to the event or censoring time of the individual
        list.times=originaldata[,timeCol]
        cumhaz=sapply(1:n,function(x){sum(cumhaz.elements[which(cut.points<=list.times[x]),x])})
        survest=exp(-cumhaz)#this is the estimated survival function
        #now assign the cumulative baseline hazard to each individual
        basecumhaz.all=rep(NA,n)
        if(originaldata$d[1]==0){
        basecumhaz.all[which(imputations[[imp]][,timeCol]%in%cut.points & duplicated(imputations[[imp]][,timeCol])==F)]=  basecumhaz[-1]} else {
          basecumhaz.all[which(imputations[[imp]][,timeCol]%in%cut.points & duplicated(imputations[[imp]][,timeCol])==F)]=  basecumhaz} #the 'duplicated bit here delas with ties
        basecumhaz.all=zoo::na.locf(basecumhaz.all,na.rm=F)
        basecumhaz.all=ifelse(is.na(basecumhaz.all)==1,0,basecumhaz.all)
        
        basehaz.all=rep(NA,n)
        basehaz.all[which(imputations[[imp]][,timeCol]%in%cut.points & duplicated(imputations[[imp]][,timeCol])==F)]=  basehaz #the 'duplicated bit here deals with ties
        basehaz.all=zoo::na.locf(basehaz.all,na.rm=F)
        basehaz.all=ifelse(is.na(basehaz.all)==1,0,basehaz.all)
        
        #end of code to obtain baseline cumulative hazard
        
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
        #DIRECT SAMPLING PART
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
          
          #generate covariate-by-time interactions in imputations[[imp]]
            imputations[[imp]]$x1t<-imputations[[imp]][,names.explanCol[1]]*imputations[[imp]][,timeCol]
            imputations[[imp]]$x2t<-imputations[[imp]][,names.explanCol[2]]*imputations[[imp]][,timeCol]
            
            imputations[[imp]]$x1.sp5.1<-imputations[[imp]][,names.explanCol[1]]*timetransform.sp5.1
            imputations[[imp]]$x1.sp5.2<-imputations[[imp]][,names.explanCol[1]]*timetransform.sp5.2
            imputations[[imp]]$x1.sp5.3<-imputations[[imp]][,names.explanCol[1]]*timetransform.sp5.3
            
            imputations[[imp]]$x2.sp5.1<-imputations[[imp]][,names.explanCol[2]]*timetransform.sp5.1
            imputations[[imp]]$x2.sp5.2<-imputations[[imp]][,names.explanCol[2]]*timetransform.sp5.2
            imputations[[imp]]$x2.sp5.3<-imputations[[imp]][,names.explanCol[2]]*timetransform.sp5.3
          
          outmodxb <-  model.matrix(as.formula(smformula3),imputations[[imp]])
          outmodxb <- outmodxb[,2:dim(outmodxb)[2]] %*% outcomeModBeta
          
          #------------
          #estimate survivor probs
          #------------
          
          #for each individual calculate the sum of basehaz*exp.t^x for all times up to the event/censoring time. This is their cumulative hazard.
          cumhaz.elements=basehaz*outer(exp.t.x1,imputations[[imp]][,explanCol[1]],FUN="^")*
            outer(exp.t.x2,imputations[[imp]][,explanCol[2]],FUN="^") 
          #cumhaz.elements are the terms at each event time, which are to be summed up to the event or censoring time of the individual
          list.times=originaldata[,timeCol]
          cumhaz=sapply(1:n,function(x){sum(cumhaz.elements[which(cut.points<=list.times[x]),x])})
          survest=exp(-cumhaz)#this is the estimated survival function
          
          outcomeDens <- survest[imputationNeeded]* (exp(outmodxb[imputationNeeded]))^d[imputationNeeded]
          
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
        #END OF THE DIRECT SAMPLING PART
        #----------------------------------------------------
        
        #----------------------------------------------------
        #START OF THE REJECTION SAMPLING PART
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
            
            #generate covariate-by-time interactions in imputations[[imp]]
            if(tvar1!="" & tvar2==""){
              imputations[[imp]]$x1t<-imputations[[imp]][,names.explanCol[1]]*imputations[[imp]][,timeCol]
              
              imputations[[imp]]$x1.sp5.1<-imputations[[imp]][,names.explanCol[1]]*timetransform.sp5.1
              imputations[[imp]]$x1.sp5.2<-imputations[[imp]][,names.explanCol[1]]*timetransform.sp5.2
              imputations[[imp]]$x1.sp5.3<-imputations[[imp]][,names.explanCol[1]]*timetransform.sp5.3
            } else if(tvar1!="" & tvar2!=""){
              imputations[[imp]]$x1t<-imputations[[imp]][,names.explanCol[1]]*imputations[[imp]][,timeCol]
              imputations[[imp]]$x2t<-imputations[[imp]][,names.explanCol[2]]*imputations[[imp]][,timeCol]
              
              imputations[[imp]]$x1.sp5.1<-imputations[[imp]][,names.explanCol[1]]*timetransform.sp5.1
              imputations[[imp]]$x1.sp5.2<-imputations[[imp]][,names.explanCol[1]]*timetransform.sp5.2
              imputations[[imp]]$x1.sp5.3<-imputations[[imp]][,names.explanCol[1]]*timetransform.sp5.3
              
              imputations[[imp]]$x2.sp5.1<-imputations[[imp]][,names.explanCol[2]]*timetransform.sp5.1
              imputations[[imp]]$x2.sp5.2<-imputations[[imp]][,names.explanCol[2]]*timetransform.sp5.2
              imputations[[imp]]$x2.sp5.3<-imputations[[imp]][,names.explanCol[2]]*timetransform.sp5.3
            }
            outmodxb <-  model.matrix(as.formula(smformula3),imputations[[imp]])
            outmodxb <- outmodxb[,2:dim(outmodxb)[2]] %*% outcomeModBeta
            
            #for each individual calculate the sum of basehaz*exp.t^x for all times up to the event/censoring time. This is their cumulative hazard.
            cumhaz.elements=basehaz*outer(exp.t.x1,imputations[[imp]][,explanCol[1]],FUN="^")*
              outer(exp.t.x2,imputations[[imp]][,explanCol[2]],FUN="^") 
            cumhaz=sapply(1:n,function(x){sum(cumhaz.elements[which(cut.points<=list.times[x]),x])})
            survest=exp(-cumhaz)#this is the estimated survival function

            prob.d0<-survest[imputationNeeded]
            prob.d1<-basehaz.all[imputationNeeded]*exp(1+outmodxb[imputationNeeded])*survest[imputationNeeded]
            prob = d[imputationNeeded]*prob.d1 + (1-d[imputationNeeded])*prob.d0
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
            
            #generate covariate-by-time interactions in tempData
              tempData$x1t<-tempData[,names.explanCol[1]]*tempData[,timeCol]
              tempData$x2t<-tempData[,names.explanCol[2]]*tempData[,timeCol]
              
              tempData$x1.sp5.1<-tempData[,names.explanCol[1]]*timetransform.sp5.1[i]
              tempData$x1.sp5.2<-tempData[,names.explanCol[1]]*timetransform.sp5.2[i]
              tempData$x1.sp5.3<-tempData[,names.explanCol[1]]*timetransform.sp5.3[i]
              
              tempData$x2.sp5.1<-tempData[,names.explanCol[2]]*timetransform.sp5.1[i]
              tempData$x2.sp5.2<-tempData[,names.explanCol[2]]*timetransform.sp5.2[i]
              tempData$x2.sp5.3<-tempData[,names.explanCol[2]]*timetransform.sp5.3[i]
            
            outmodxb <-  model.matrix(as.formula(smformula3),tempData)
            outmodxb <- outmodxb[,2:dim(outmodxb)[2]] %*% outcomeModBeta
              
            
            #for each individual calculate the sum of basehaz*exp.t^x for all times up to the event/censoring time. This is their cumulative hazard.
            cumhaz.elements=basehaz*outer(exp.t.x1,tempData[,explanCol[1]],FUN="^")*
              outer(exp.t.x2,tempData[,explanCol[2]],FUN="^") 
            cumhaz=sapply(1:rjlimit,function(x){sum(cumhaz.elements[which(cut.points<=tempData[,timeCol][x]),x])})
            survest.rej=exp(-cumhaz)#this is the estimated survival function
            
            prob.d0<-survest.rej
            prob.d1<-basehaz.all[i]*exp(1+outmodxb)*survest.rej
            prob = d[i]*prob.d1 + (1-d[i])*prob.d0
            reject = 1*(uDraw > prob )
            
            if (sum(reject)<rjlimit) {
              imputations[[imp]][i,targetCol] <- tempData[reject==0,targetCol][1]
            } else {
              print("Rejection sampling has failed for one record. You may want to increase the rejecton sampling limit.")
            }
          }
          #update passive variables
          imputations[[imp]] <- updatePassiveVars(imputations[[imp]], method, passiveVars)
        }
        
        #END OF THE REJECTION SAMPLING PART
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