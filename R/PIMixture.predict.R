#cumulative.risk is called in ipw.lc.splines
cumulative.risk2<-function(prevalence,Lambda,cox.risk){
  prevalence+(1-prevalence)*(1-exp(-Lambda*exp(cox.risk)))
}  #end of function cumulative.risk2

#cum.risk can be applied to the result of ipw.lc.splines
cumulative.risk<-function(logit.coef,cox.coef,logit.predictor,cox.predictor=NULL,Lambda.data,CR.time.points,knots,order,spline.para.est,cov.mat=NULL,...){

  if(!is.null(cox.predictor)){
    if (class(logit.predictor)=="matrix"&class(cox.predictor)=="numeric"){
      stop("cox.predictor should be matrix type.")
    }
  }

  if(class(logit.predictor)=="matrix"){
    if(ncol(logit.predictor)!=length(logit.coef)){
      stop("logit.predictor does not match with the logit.coef")
    }
  }
  if(class(logit.predictor)=="numeric"&length(logit.predictor)!=length(logit.coef)){
    stop("logit.predictor does not match with the logit.coef")
  }
  if(!is.null(cox.predictor)){
    if(class(cox.predictor)=="numeric"&length(cox.predictor)!=length(cox.coef)){
      stop("cox.predictor does not match with the cox.coef.")
    }
    if(class(cox.predictor)=="matrix"){
      if(ncol(cox.predictor)!=length(cox.coef)){
        stop("cox.predictor does not match with the cox.coef.")
      }
    }
  }
  if(class(Lambda.data)!="data.frame"){
    cat(Lambda.data)
    stop("Please put \"cum.hazard\" of the output list from ipw.lc.splines or ipw.lc.semipara into Lambda.data ")
  }
  if(sum(as.numeric(c("time", "cum.hazard") %in% names(Lambda.data)))!=2){
    stop("Lambda.data should include time and cum.hazard columns")
  }

  if(is.null(cox.predictor)){
    betax<-logit.predictor%*%cbind(as.numeric(logit.coef))
    gammaz<-0
  } else {
    betax<-  logit.predictor%*%cbind(as.numeric(logit.coef))
    gammaz<- cox.predictor%*%cbind(as.numeric(cox.coef))
  }


  indx<-which(Lambda.data$time %in% CR.time.points)
  base.chazard<-as.numeric(Lambda.data$cum.hazard[indx])

  CR.out<-NULL

  for(i in 1:nrow(betax)){
    prevalence<-exp(betax[i,])/(1+exp(betax[i,]))

    if(is.null(cox.predictor)) {
      cdf<-1-exp(-base.chazard)
    } else {
      cdf<-1-exp(-base.chazard*exp(gammaz[i,]))
    }

    lpred<-paste0(logit.predictor)
    cpred<-paste0(cox.predictor)
    CR=prevalence+(1-prevalence)*cdf
    lpred<-paste(logit.predictor[i,],collapse=",")

    if(!is.null(cox.predictor)){
      cpred<-paste(cox.predictor[i,],collapse=",")
    }else if(is.null(cox.predictor)){
      cpred<-"NA"
    }

    if(!is.null(cov.mat)){
      if(!is.null(cox.predictor)){
        CRSE<-CR.se(knots,order, CR.time.points, cov.mat,spline.para.est,
                    logit.risk.est=as.numeric(betax[i,]),cox.risk.est=as.numeric(gammaz[i,]),
                    logit.predictor=as.numeric(logit.predictor[i,]),cox.predictor=as.numeric(cox.predictor[i,]))
        }else if(is.null(cox.predictor)){
          CRSE<-CR.se.noGamma(knots,order, CR.time.points,cov.mat, spline.para.est,
          logit.risk.est=as.numeric(betax[i,]),as.numeric(logit.predictor[i,]))
        }

        UL95=CR+1.96*CRSE
        LL95=CR-1.96*CRSE
        CR.out1<-data.frame(logit.predictor=lpred,cox.predictor=cpred,
                            time=as.numeric(CR.time.points),CR=as.numeric(CR),
                            CR.se=as.numeric(CRSE),LL95=LL95,UL95=UL95)
    } else if(is.null(cov.mat)){
      CR.out1<-data.frame(logit.predictor=lpred,cox.predictor=cpred, time=as.numeric(CR.time.points),CR=as.numeric(CR))
    }

    CR.out<-rbind(CR.out,CR.out1)
  }
  return(CR.out)
}    #end of function cumulative.risk

CR.se<-function(knots2,order, CR.time.points, cov.mat,spline.para.est, logit.risk.est,cox.risk.est,logit.predictor,cox.predictor){

  ###Spline-based Lambda Calculation #############

  Ibasis<-Ispline(as.numeric(CR.time.points), order = order, knots = knots2)

  basis.mat<-exp(spline.para.est)*Ibasis

  if(length(CR.time.points)==1){
    Lambda<-as.numeric(sum(basis.mat))
  }
  if(length(CR.time.points)>1){
    Lambda<-as.numeric(colSums(basis.mat))
  }


  cons1<-as.numeric(exp(logit.risk.est)/(1+exp(logit.risk.est))^2)
  cons2<-as.numeric(1/(1+exp(logit.risk.est)))
  SS<-exp(-Lambda*exp(cox.risk.est))
  FF<- 1-SS


  cr.beta<-(cons1*as.numeric(logit.predictor))%*%t(SS)
  cr.gamma<- (cons2*exp(cox.risk.est)*cox.predictor) %*%t((SS*Lambda))
  cr.b<- (cons2*exp(cox.risk.est)*SS) *t(basis.mat)


  if(length(CR.time.points)==1){
    cr.delta<-c(cr.beta,cr.gamma,cr.b)
  }
  if(length(CR.time.points)>1){
    cr.delta<-rbind(cr.beta,cr.gamma,t(cr.b))   # number of parameters X number of times
  }

  cr.cov<-t(cr.delta) %*%cov.mat%*%cr.delta
  cr.se<-sqrt(diag(cr.cov))
  return(cr.se)
} #end of function CR.se

CR.se.noGamma<-function(knots2,order, CR.time.points,cov.mat, spline.para.est, logit.risk.est,logit.predictor){

  ###Spline-based Lambda Calculation #############

  Ibasis<-Ispline(as.numeric(CR.time.points), order = order, knots = knots2)

  basis.mat<-exp(spline.para.est)*Ibasis

  if(length(CR.time.points)==1){
    Lambda<-as.numeric(sum(basis.mat))
  }
  if(length(CR.time.points)>1){
    Lambda<-as.numeric(colSums(basis.mat))
  }

  cox.risk.est<-0
  cons1<-as.numeric(exp(logit.risk.est)/(1+exp(logit.risk.est))^2)
  cons2<-as.numeric(1/(1+exp(logit.risk.est)))
  SS<-exp(-Lambda*exp(cox.risk.est))

  FF<- 1-SS

  cr.beta<-(cons1*as.numeric(logit.predictor))%*%t(SS)
  cr.b<-(cons2*exp(cox.risk.est)*SS) *t(basis.mat)

  if(length(CR.time.points)==1){
    cr.delta<-c(cr.beta,cr.b)
  }
  if(length(CR.time.points)>1){
    cr.delta<-rbind(cr.beta,t(cr.b))  # number of parameters X number of times
  }

  cr.cov<-t(cr.delta) %*%cov.mat%*%cr.delta
  cr.se<-sqrt(diag(cr.cov))
  return(cr.se)
}  #end of function CR.se.noGamma

##################################################
# get design matrix for a model, the
# model: the model going to be used
# data:  the dataset which include the variables
##################################################

get.design.mat<- function(data, model) {
  if(is.null(model)) {
    return(NULL)
  } else {
    bb<-strsplit(model,split="~")
    p.model<-paste0("~",bb[[1]][2])
    fml<-as.formula(p.model)

    mf <- model.frame(formula=fml, data=data)
    design.mat <- model.matrix(attr(mf, "terms"), data=mf)
    return(design.mat)
  }
}

##################################################
# PIMixture prediction for parametric model
##################################################
PIMixture.predict.param <- function(x,times,mat1=NULL,mat2=NULL){
   requireNamespace("survival", quietly = TRUE)
   requireNamespace("flexsurv", quietly = TRUE)
   requireNamespace("interval", quietly = TRUE)

   x <- as.list(x)
   if(!"model"%in%names(x)) {
     stop("x must have a \'model\' element")
   }

   #########################################################
   ## cumulative risk prediction for non-parametric models
   #########################################################
   if (x$model == "non-parametric") {
     nonparam_matrix <- matrix(NA,length(times),4)
     nonparam_matrix[,1] <- times
     for (i in 1:length(times))
     {
       nonparam_matrix[i,2] <- cumsum(x$pf)[max(which(x$intmap[2,] <= max(times[i],.01)))]
       if (x$conf.int=="TRUE")
       {
         nonparam_matrix[i,3] <- 1-x$CI$upper[max(which(x$CI$time <= max(times[i],.01)))]
         nonparam_matrix[i,4] <- 1-x$CI$lower[max(which(x$CI$time <= max(times[i],.01)))]
       }
     }

     colnames(nonparam_matrix) <- c("time","CR","LL95","UL95")
     if(x$conf.int=="TRUE") {
       pred<- nonparam_matrix
     } else {
       pred<- matrix(nonparam_matrix[, 1:2], nrow=nrow(nonparam_matrix))
       colnames(pred)<- colnames(nonparam_matrix)[1:2]
     }

     return(pred)
   }

   #########################################################################
   ## The rests are cumulative risk prediction for parametric models
   #########################################################################

   m1 <- 0
   m2 <- 0
   npred <- 1

   if (!is.null(mat1)){
      mat1 <- as.matrix(mat1)
      m1 <- ncol(mat1)
      npred <- nrow(mat1)
   }
   if (!is.null(mat2)){
      mat2 <- as.matrix(mat2)
      m2 <- ncol(mat2)
      if (!is.null(mat1))
      { if (nrow(mat2) != npred) { stop("number of rows for mat1 and mat2 do not agree") }
      } else { npred <- nrow(mat2) }
   }

   n1<- 1
   n2<- 1
   n3<- 1


   if(x$model != "non-parametric") {
     # for parametric model, the parametric estimates are stored in regression.coef
     # also n1, n2, n3 are calculated based on regression.coef

     if ( !c("regression.coef") %in% names(x)) {
       stop("x must include \'regression.coef'\ element")
     }

     x$par<- x$regression.coef[, 3]
     x$hessian<- -x$hessian   #convert the hessian to negative hessian, for covariance calculation

     if (c("data.summary")%in%names(x)) {
       n1 <- x$data.summary[x$data.summary[,1] == "Known prevalent cases", 2] #Prevalent cases
       n2 <- x$data.summary[x$data.summary[,1] == "Interval censoring", 2] #incidence cases (Interval censored)
       n3 <- x$data.summary[x$data.summary[,1] == "Left censoring",     2] #Unknown (left censored)
     }

     if(sum(tolower(x$regression.coef[,1]) == "logit") ==0) {
       #incidence model only
       n2<- 0
       n3<- 0
     }

     if(sum(tolower(x$regression.coef[,1]) == "logit") == length(x$regression.coef[,1])) {
       #prevalent model only
       n1<- 0
       n3<- 0
     }
   }

   logistic_weibull_varEst <- function(t)
   {
      if (n1==0 & n3==0)   #incident disease only
      {
         if (m2==0) {
            s_lp <- x$par[1]
            mat2 <- matrix(1,npred,1)
         }
         else if (m2>0)
         {
            s_lp <- x$par[1]+mat2%*%x$par[2:(m2+1)]
            mat2 <- cbind(rep(1,npred),mat2)
         }
         tau <- x$par[(m2+2)]
         ST <- 1-pweibull(t,1/tau,exp(s_lp))

         grad <- matrix(NA,npred,m2+2)
         for (i in 1:(m2+1)) {grad[,i] <- -(t^(1/tau))*ST / (tau*exp(s_lp/tau)) * mat2[,i] }
         grad[,m2+2] <- (s_lp-log(t))*(t^(1/tau))*ST / ((tau^2)*exp(s_lp/tau))
         varcr <- rep(NA,npred)
         for (i in 1:npred){ varcr[i] <- t(grad[i,])%*%covvar%*%grad[i,] }
      } else if (n2==0 & n3==0)
      {
         if (m1==0)
         {
            p_lp <- x$par[1]
            p <- exp(p_lp) / (1+exp(p_lp))
            varcr <- p**2 / (1+exp(p_lp))**2 * covvar
         } else if (m1==1) {
            p_lp <- x$par[1]+mat1*x$par[2]
            p <- exp(p_lp) / (1+exp(p_lp))
            mat1 <- cbind(rep(1,npred),mat1)

            grad <- matrix(NA,npred,m1+1)
            for (i in 1:(m1+1)) {grad[,i] <- p / (1+exp(p_lp)) * mat1[,i] }
            varcr <- rep(NA,npred)
            for (i in 1:npred){ varcr[i] <- t(grad[i,])%*%covvar%*%grad[i,] }
         } else if (m1>0) {
            p_lp <- x$par[1]+mat1%*%x$par[2:(m1+1)]
            p <- exp(p_lp) / (1+exp(p_lp))
            mat1 <- cbind(rep(1,npred),mat1)

            grad <- matrix(NA,npred,m1+1)
            for (i in 1:(m1+1)) {grad[,i] <- p / (1+exp(p_lp)) * mat1[,i] }
            varcr <- rep(NA,npred)
            for (i in 1:npred){ varcr[i] <- t(grad[i,])%*%covvar%*%grad[i,] }
         }
      } else
      {
         if (m1==0)
         {
            p_lp <- x$par[1]
            p <- exp(p_lp) / (1+exp(p_lp))
            mat1 <- matrix(1,npred,1)
         } else if (m1==1) {
            p_lp <- x$par[1]+mat1*x$par[2]
            p <- exp(p_lp) / (1+exp(p_lp))
            mat1 <- cbind(rep(1,npred),mat1)
         } else if (m1>1) {
            p_lp <- x$par[1]+mat1%*%x$par[2:(m1+1)]
            p <- exp(p_lp) / (1+exp(p_lp))
            mat1 <- cbind(rep(1,npred),mat1)
         }

         if (m2==0) {
            s_lp <- x$par[(m1+2)]
            mat2 <- matrix(1,npred,1)
         } else if (m2>0)
         {
            s_lp <- x$par[(m1+2)]+mat2%*%x$par[(m1+3):(m1+m2+2)]
            mat2 <- cbind(rep(1,npred),mat2)
         }
         tau <- x$par[(m1+m2+3)]
         ST <- 1-pweibull(t,1/tau,exp(s_lp))

         if (t==0)
         {
            grad <- matrix(NA,npred,m1+1)
            for (i in 1:(m1+1)) {grad[,i] <- p*ST / (1+exp(p_lp)) * mat1[,i] }
            varcr <- rep(NA,npred)
            for (i in 1:npred){ varcr[i] <- t(grad[i,])%*%covvar[1:(m1+1),1:(m1+1)]%*%grad[i,] }
         } else if (t>0)
         {
            grad <- matrix(NA,npred,m1+m2+3)
            for (i in 1:(m1+1)) {grad[,i] <- p*ST / (1+exp(p_lp)) * mat1[,i] }
            for (i in 1:(m2+1)) {grad[,(m1+1+i)] <- -(t^(1/tau))*ST*(1-p) / (tau*exp(s_lp/tau)) * mat2[,i] }
            grad[,m1+m2+3] <- (s_lp-log(t))*(t^(1/tau))*ST*(1-p) / ((tau^2)*exp(s_lp/tau))
            varcr <- rep(NA,npred)
            for (i in 1:npred){ varcr[i] <- t(grad[i,])%*%covvar%*%grad[i,] }
         }
      }
      names(varcr) <- paste(t)
      return(varcr)
   }

   logistic_exponential_varEst <- function(t)
   {
      if (n1==0 & n3==0)
      {
         if (m2==0) {
            s_lp <- x$par[1]
            mat2 <- matrix(1,npred,1)
            ST <- 1-pweibull(t,1,exp(s_lp))
            grad <- -t*ST / exp(s_lp)
            varcr <- grad^2*covvar
         }
         else if (m2>0)
         {
            s_lp <- x$par[1]+mat2%*%x$par[2:(m2+1)]
            mat2 <- cbind(rep(1,npred),mat2)
            ST <- 1-pweibull(t,1,exp(s_lp))
            grad <- matrix(NA,npred,m2+1)
            for (i in 1:(m2+1)) {grad[,i] <- -t*ST / exp(s_lp) * mat2[,i] }
            varcr <- rep(NA,npred)
            for (i in 1:npred){ varcr[i] <- t(grad[i,])%*%covvar%*%grad[i,] }
         }
      } else if (n2==0 & n3==0)
      {
         if (m1==0)
         {
            p_lp <- x$par[1]
            p <- exp(p_lp) / (1+exp(p_lp))
            varcr <- p**2 / (1+exp(p_lp))**2 * covvar
         } else if (m1>0) {
            p_lp <- x$par[1]+mat1%*%x$par[2:(m1+1)]
            p <- exp(p_lp) / (1+exp(p_lp))
            mat1 <- cbind(rep(1,npred),mat1)

            grad <- matrix(NA,npred,m1+1)
            for (i in 1:(m1+1)) {grad[,i] <- p / (1+exp(p_lp)) * mat1[,i] }
            varcr <- rep(NA,npred)
            for (i in 1:npred){ varcr[i] <- t(grad[,i])%*%covvar%*%grad[,i] }
         }
      } else
      {
         if (m1==0)
         {
            p_lp <- x$par[1]
            p <- exp(p_lp) / (1+exp(p_lp))
            mat1 <- matrix(1,npred,1)
         } else if (m1>0) {
            p_lp <- x$par[1]+mat1%*%x$par[2:(m1+1)]
            p <- exp(p_lp) / (1+exp(p_lp))
            mat1 <- cbind(rep(1,npred),mat1)
         }

         if (m2==0) {
            s_lp <- x$par[(m1+2)]
            mat2 <- matrix(1,npred,1)
         } else if (m2>0)
         {
            s_lp <- x$par[(m1+2)]+mat2%*%x$par[(m1+3):(m1+m2+2)]
            mat2 <- cbind(rep(1,npred),mat2)
         }
         ST <- 1-pweibull(t,1,exp(s_lp))

         if (t==0)
         {
            grad <- matrix(NA,npred,m1+1)
            for (i in 1:(m1+1)) {grad[,i] <- p*ST / (1+exp(p_lp)) * mat1[,i] }
            varcr <- rep(NA,npred)
            for (i in 1:npred){ varcr[i] <- t(grad[i,])%*%covvar[1:(m1+1),1:(m1+1)]%*%grad[i,] }
         } else if (t>0)
         {
            grad <- matrix(NA,npred,m1+m2+2)
            for (i in 1:(m1+1)) {grad[,i] <- p*ST / (1+exp(p_lp)) * mat1[,i] }
            for (i in 1:(m2+1)) {grad[,(m1+1+i)] <- -t*ST*(1-p) / exp(s_lp) * mat2[,i] }
            varcr <- rep(NA,npred)
            for (i in 1:npred){ varcr[i] <- t(grad[i,])%*%covvar%*%grad[i,] }
         }
      }
      names(varcr) <- paste(t)
      return(varcr)
   }

   pllog <- function (q, shape = 1, scale = 1, lower.tail = TRUE, log.p = FALSE)
   {
       Fx <- plogis(log(q), location = scale, scale = shape)
       if (!lower.tail)
           Fx <- 1 - Fx
       if (log.p)
           Fx <- log(Fx)
       return(Fx)
   }

   varcumrisk <- function(y)
   {
      if (x$model=="logistic-Weibull")
      {
         logistic_weibull_varEst(y)
      } else if (x$model=="logistic-exponential")
      {
         logistic_exponential_varEst(y)
      }
   }

   if (x$model=="logistic-Weibull")
   {
      if (n1==0 & n3==0) { if (length(x$par) != (m2+2)) stop("Number of covariates not consistent with the number of parameters")
      } else if (n2==0 & n3==0) { if (length(x$par) != (m1+1)) stop("Number of covariates not consistent with the number of parameters")
      } else { if (length(x$par) != (m1+m2+3)) stop("Number of covariates not consistent with the number of parameters") }

      if (n2==0 & n3==0 & m1==0) { covvar <- 1/x$hessian
      } else { covvar <- solve(x$hessian) }

      if (n1==0 & n3==0)
      {
         if (m2==0) {s_lp <- x$par[1]
         } else if (m2>0) {s_lp <- x$par[1]+mat2%*%x$par[2:(m2+1)] }
         tau <- x$par[(m2+2)]
         predrisk <- function(y) {pweibull(y,1/tau,exp(s_lp)) }
         pred_est <- sapply(times,predrisk)
         if (npred == 1) { names(pred_est) <- paste(times)
         } else { colnames(pred_est) <- paste(times) }
      } else if (n2==0 & n3==0)
      {
         if (m1==0) {p_lp <- x$par[1]
         } else if (m1==1) {p_lp <- x$par[1]+mat1*x$par[2]
         } else if (m1>1) {p_lp <- x$par[1]+mat1%*%x$par[2:(m1+1)] }
         pred_est <- matrix(nrow=npred,ncol=length(times),exp(p_lp) / (1+exp(p_lp)))
         if (npred == 1) { names(pred_est) <- paste(times)
         } else { colnames(pred_est) <- paste(times) }
      } else
      {
         if (m1==0) {p_lp <- x$par[1]
         } else if (m1==1) {p_lp <- x$par[1]+mat1*x$par[2]
         } else if (m1>0) {p_lp <- x$par[1]+mat1%*%x$par[2:(m1+1)] }

         if (m2==0) {s_lp <- x$par[(m1+2)]
         } else if (m2>0) {s_lp <- x$par[(m1+2)]+mat2%*%x$par[(m1+3):(m1+m2+2)] }

         tau <- x$par[(m1+m2+3)]
         p <- exp(p_lp) / (1+exp(p_lp))
         predrisk <- function(y) { p+(1-p)*pweibull(y,1/tau,exp(s_lp)) }
         pred_est <- sapply(times,predrisk)
         if (npred == 1) { names(pred_est) <- paste(times)
         } else { colnames(pred_est) <- paste(times) }
      }

      var_pred <- sapply(times,varcumrisk)
      pred_lcl <- 1-exp(-exp(log(-log(1-pred_est))-1.96*sqrt(var_pred)/sqrt((log(1-pred_est)*(1-pred_est))**2)))
      if (npred == 1) { names(pred_lcl) <- paste(times)
      } else { colnames(pred_lcl) <- paste(times) }
      pred_ucl <- 1-exp(-exp(log(-log(1-pred_est))+1.96*sqrt(var_pred)/sqrt((log(1-pred_est)*(1-pred_est))**2)))
      if (npred == 1) { names(pred_ucl) <- paste(times)
      } else
      {
         colnames(var_pred) <- paste(times)
         rownames(var_pred) <- rownames(mat1)
         colnames(pred_ucl) <- paste(times)
      }
      pred <- list(pred_est,sqrt(var_pred),pred_lcl,pred_ucl)
      names(pred) <- c("CR","CR.se","LL95","UL95")
   } else if (x$model=="logistic-exponential")
   {
      if (n1==0 & n3==0) { if (length(x$par) != (m2+1)) stop("Number of covariates not consistent with the number of parameters")
      } else if (n2==0 & n3==0) { if (length(x$par) != (m2+1)) stop("Number of covariates not consistent with the number of parameters")
      } else { if (length(x$par) != (m1+m2+2)) stop("Number of covariates not consistent with the number of parameters") }

      if ((n1==0 & n3==0 & m2==0) | (n2==0 & n3==0 & m1==0)) { covvar <- 1/x$hessian
      } else { covvar <- solve(x$hessian) }

      if (n1==0 & n3==0)
      {
         if (m2==0) {s_lp <- x$par[1]
         } else if (m2>0) {s_lp <- x$par[1]+mat2%*%x$par[2:(m2+1)] }
         predrisk <- function(y) {pweibull(y,1,exp(s_lp)) }
         pred_est <- sapply(times,predrisk)
         if (npred == 1) { names(pred_est) <- paste(times)
         } else { colnames(pred_est) <- paste(times) }
      } else if (n2==0 & n3==0)
      {
         if (m1==0) {p_lp <- x$par[1]
         } else if (m1>0) {p_lp <- x$par[1]+mat1%*%x$par[2:(m1+1)] }
         pred_est <- exp(p_lp) / (1+exp(p_lp))
         if (npred == 1) { names(pred_est) <- paste(times)
         } else { colnames(pred_est) <- paste(times) }
      } else
      {
         if (m1==0) {p_lp <- x$par[1]
         } else if (m1>0) {p_lp <- x$par[1]+mat1%*%x$par[2:(m1+1)] }

         if (m2==0) {s_lp <- x$par[(m1+2)]
         } else if (m2>0) {s_lp <- x$par[(m1+2)]+mat2%*%x$par[(m1+3):(m1+m2+2)] }

         p <- exp(p_lp) / (1+exp(p_lp))
         predrisk <- function(y) { p+(1-p)*pweibull(y,1,exp(s_lp)) }
         pred_est <- sapply(times,predrisk)
         if (npred == 1) { names(pred_est) <- paste(times)
         } else { colnames(pred_est) <- paste(times) }
      }

      var_pred <- sapply(times,varcumrisk)
      pred_lcl <- 1-exp(-exp(log(-log(1-pred_est))-1.96*sqrt(var_pred)/sqrt((log(1-pred_est)*(1-pred_est))**2)))
      if (npred == 1) { names(pred_lcl) <- paste(times)
      } else { colnames(pred_lcl) <- paste(times) }
      pred_ucl <- 1-exp(-exp(log(-log(1-pred_est))+1.96*sqrt(var_pred)/sqrt((log(1-pred_est)*(1-pred_est))**2)))
      if (npred == 1) { names(pred_ucl) <- paste(times)
      } else
      {
         colnames(var_pred) <- paste(times)
         rownames(var_pred) <- rownames(mat1)
         colnames(pred_ucl) <- paste(times)
      }
      pred <- list(pred_est,sqrt(var_pred),pred_lcl,pred_ucl)
      names(pred) <- c("CR","CR.se","LL95","UL95")
   } else if (x$model=="logistic-loglogistic")
   {
      if (n1==0 & n3==0) { if (length(x$par) != (m2+2)) stop("Number of covariates not consistent with the number of parameters")
      } else if (n2==0 & n3==0) { if (length(x$par) != (m2+1)) stop("Number of covariates not consistent with the number of parameters")
      } else { if (length(x$par) != (m1+m2+3)) stop("Number of covariates not consistent with the number of parameters") }

      if (n2==0 & n3==0 & m1==0) { covvar <- 1/x$hessian
      } else { covvar <- solve(x$hessian) }

      if (n1==0 & n3==0)
      {
         if (m2==0) {s_lp <- x$par[1]
         } else if (m2>0) {s_lp <- x$par[1]+mat2%*%x$par[2:(m2+1)] }
         tau <- x$par[(m2+2)]
         predrisk <- function(y) {pred_est <- pllog(y,tau,s_lp) }
         pred_est <- sapply(times,predrisk)
         if (npred == 1) { names(pred_est) <- paste(times)
         } else { colnames(pred_est) <- paste(times) }
      } else if (n2==0 & n3==0)
      {
         if (m1==0) {p_lp <- x$par[1]
         } else if (m1>0) {p_lp <- x$par[1]+mat1%*%x$par[2:(m1+1)] }
         pred_est <- exp(p_lp) / (1+exp(p_lp))
         if (npred == 1) { names(pred_est) <- paste(times)
         } else { colnames(pred_est) <- paste(times) }
      } else
      {
         if (m1==0) {p_lp <- x$par[1]
         } else if (m1>0) {p_lp <- x$par[1]+mat1%*%x$par[2:(m1+1)] }

         if (m2==0) {s_lp <- x$par[(m1+2)]
         } else if (m2>0) {s_lp <- x$par[(m1+2)]+mat2%*%x$par[(m1+3):(m1+m2+2)] }

         tau <- x$par[(m1+m2+3)]
         p <- exp(p_lp) / (1+exp(p_lp))
         predrisk <- function(y) { p+(1-p)*pllog(y,tau,s_lp) }
         pred_est <- sapply(times,predrisk)
         if (npred == 1) { names(pred_est) <- paste(times)
         } else { colnames(pred_est) <- paste(times) }
      }

      pred<- list("CR"=pred_est)
   } else if (x$model=="logistic-lognormal")
   {
      if (n1==0 & n3==0) { if (length(x$par) != (m2+2)) stop("Number of covariates not consistent with the number of parameters")
      } else if (n2==0 & n3==0) { if (length(x$par) != (m2+1)) stop("Number of covariates not consistent with the number of parameters")
      } else { if (length(x$par) != (m1+m2+3)) stop("Number of covariates not consistent with the number of parameters") }

      if (n2==0 & n3==0 & m1==0) { covvar <- 1/x$hessian
      } else { covvar <- solve(x$hessian) }

      if (n1==0 & n3==0)
      {
         if (m2==0) {s_lp <- x$par[1]
         } else if (m2>0) {s_lp <- x$par[1]+mat2%*%x$par[2:(m2+1)] }
         sigma <- x$par[(m2+2)]
         predrisk <- function(y) {pred_est <- plnorm(y,s_lp,sigma) }
         pred_est <- sapply(times,predrisk)
         if (npred == 1) { names(pred_est) <- paste(times)
         } else { colnames(pred_est) <- paste(times) }
      } else if (n2==0 & n3==0)
      {
         if (m1==0) {p_lp <- x$par[1]
         } else if (m1>0) {p_lp <- x$par[1]+mat1%*%x$par[2:(m1+1)] }
         pred_est <- exp(p_lp) / (1+exp(p_lp))
         if (npred == 1) { names(pred_est) <- paste(times)
         } else { colnames(pred_est) <- paste(times) }
      } else
      {
         if (m1==0) {p_lp <- x$par[1]
         } else if (m1>0) {p_lp <- x$par[1]+mat1%*%x$par[2:(m1+1)] }

         if (m2==0) {s_lp <- x$par[(m1+2)]
         } else if (m2>0) {s_lp <- x$par[(m1+2)]+mat2%*%x$par[(m1+3):(m1+m2+2)] }

         sigma <- x$par[(m1+m2+3)]
         p <- exp(p_lp) / (1+exp(p_lp))
         predrisk <- function(y) { p+(1-p)*plnorm(y,s_lp,sigma) }
         pred_est <- sapply(times,predrisk)
         if (npred == 1) { names(pred_est) <- paste(times)
         } else { colnames(pred_est) <- paste(times) }
      }
     pred<- list("CR"=pred_est)

   } else if (x$model=="logistic-gengamma")
   {
      if (n1==0 & n3==0) { if (length(x$par) != (m2+3)) stop("Number of covariates not consistent with the number of parameters")
      } else if (n2==0 & n3==0) { if (length(x$par) != (m2+1)) stop("Number of covariates not consistent with the number of parameters")
      } else { if (length(x$par) != (m1+m2+4)) stop("Number of covariates not consistent with the number of parameters") }

      if (n2==0 & n3==0 & m1==0) { covvar <- 1/x$hessian
      } else { covvar <- solve(x$hessian) }

      if (n1==0 & n3==0)
      {
         if (m2==0) {s_lp <- x$par[1]
         } else if (m2>0) {s_lp <- x$par[1]+mat2%*%x$par[2:(m2+1)] }
         sigma <- x$par[(m2+2)]
         G <- x$par[(m2+3)]
         predrisk <- function(y) {pgengamma(y,s_lp,sigma,G)}
         pred_est <- sapply(times,predrisk)
         if (npred == 1) { names(pred_est) <- paste(times)
         } else { colnames(pred_est) <- paste(times) }
      } else if (n2==0 & n3==0)
      {
         if (m1==0) {p_lp <- x$par[1]
         } else if (m1>0) {p_lp <- x$par[1]+mat1%*%x$par[2:(m1+1)] }
         pred_est <- exp(p_lp) / (1+exp(p_lp))
         if (npred == 1) { names(pred_est) <- paste(times)
         } else { colnames(pred_est) <- paste(times) }
      } else
      {
         if (m1==0) {p_lp <- x$par[1]
         } else if (m1>0) {p_lp <- x$par[1]+mat1%*%x$par[2:(m1+1)] }

         if (m2==0) {s_lp <- x$par[(m1+2)]
         } else if (m2>0) {s_lp <- x$par[(m1+2)]+mat2%*%x$par[(m1+3):(m1+m2+2)] }

         sigma <- x$par[(m1+m2+3)]
         G <- x$par[(m1+m2+4)]
         p <- exp(p_lp) / (1+exp(p_lp))
         predrisk <- function(y) { p+(1-p)*pgengamma(y,s_lp,sigma,G) }
         pred_est <- sapply(times,predrisk)
         if (npred == 1) { names(pred_est) <- paste(times)
         } else { colnames(pred_est) <- paste(times) }
      }

     pred<- list("CR"=pred_est)
   } else if (x$model=="logistic-gamma")
   {
      if (n1==0 & n3==0) { if (length(x$par) != (m2+2)) stop("Number of covariates not consistent with the number of parameters")
      } else if (n2==0 & n3==0) { if (length(x$par) != (m2+1)) stop("Number of covariates not consistent with the number of parameters")
      } else { if (length(x$par) != m1+m2+3) stop("Number of covariates not consistent with the number of parameters") }

      if (n2==0 & n3==0 & m1==0) { covvar <- 1/x$hessian
      } else { covvar <- solve(x$hessian) }

      if (n1==0 & n3==0)
      {
         if (m2==0) {s_lp <- x$par[1]
         } else if (m2>0) {s_lp <- x$par[1]+mat2%*%x$par[2:(m2+1)] }
         tau <- x$par[(m2+2)]
         predrisk <- function(y) {pgamma(y,tau,exp(-s_lp)*tau) }
         pred_est <- sapply(times,predrisk)
         if (npred == 1) { names(pred_est) <- paste(times)
         } else { colnames(pred_est) <- paste(times) }
      } else if (n2==0 & n3==0)
      {
         if (m1==0) {p_lp <- x$par[1]
         } else if (m1>0) {p_lp <- x$par[1]+mat1%*%x$par[2:(m1+1)] }
         pred_est <- exp(p_lp) / (1+exp(p_lp))
         if (npred == 1) { names(pred_est) <- paste(times)
         } else { colnames(pred_est) <- paste(times) }
      } else
      {
         if (m1==0) {p_lp <- x$par[1]
         } else if (m1>0) {p_lp <- x$par[1]+mat1%*%x$par[2:(m1+1)] }

         if (m2==0) {s_lp <- x$par[(m1+2)]
         } else if (m2>0) {s_lp <- x$par[(m1+2)]+mat2%*%x$par[(m1+3):(m1+m2+2)] }

         tau <- x$par[(m1+m2+3)]
         p <- exp(p_lp) / (1+exp(p_lp))
         predrisk <- function(y) { p+(1-p)*pgamma(y,tau,exp(-s_lp)*tau) }
         pred_est <- sapply(times,predrisk)
         if (npred == 1) { names(pred_est) <- paste(times)
         } else { colnames(pred_est) <- paste(times) }
      }
     pred<- list("CR"=pred_est)
   }

   #####################################################
   # rearrange the pred into matrix format
   #####################################################

   if(!is.null(mat1)) {
     xmat1<- cbind(1, mat1)
   } else {
     xmat1<- 1
   }

   if(!is.null(mat2)) {
     xmat2<- cbind(1, mat2)
   } else {
     xmat2<- 1
   }

   pred.xmat1<-  vector()
   pred.xmat2<-  vector()
   CR.out<-NULL

   if(npred==1) {
      pred.xmat1<- rep(paste(xmat1, collapse=","), length(times))
      pred.xmat2<- rep(paste(xmat2, collapse=","), length(times))
      if("CR.se" %in% names(pred)) {
        CR.out<-data.frame(pred.xmat1, pred.xmat2, times, pred$CR, pred$CR.se, pred$LL95, pred$UL95,
                           check.names=FALSE, stringsAsFactors=FALSE)
        colnames(CR.out)<- c("prev.predictor", "survival.predictor", "time", "CR", "CR.se", "LL95", "UL95")
      } else {
        CR.out<-data.frame(pred.xmat1, pred.xmat2, times, pred$CR, check.names=FALSE, stringsAsFactors=FALSE)
        colnames(CR.out)<- c("prev.predictor", "survival.predictor", "time", "CR")
      }
   } else {
     for(i in 1:nrow(pred$CR)) {
        pred.xmat1<- rep(paste(xmat1[i, ], collapse=","), length(times))
        pred.xmat2<- rep(paste(xmat2[i, ], collapse=","), length(times))
        if("CR.se" %in% names(pred)) {
          out1<-data.frame(pred.xmat1, pred.xmat2, times, pred$CR[i, ], pred$CR.se[i, ],
                           pred$LL95[i, ], pred$UL95[i, ], check.names=FALSE, stringsAsFactors=FALSE)
          colnames(out1)<- c("prev.predictor", "survival.predictor", "time", "CR", "CR.se", "LL95", "UL95")
        } else {
          out1<-data.frame(pred.xmat1, pred.xmat2, times, pred$CR, check.names=FALSE, stringsAsFactors=FALSE)
          colnames(out1)<- c("prev.predictor", "survival.predictor", "time", "CR")
        }
        CR.out<- rbind(CR.out, out1)
      }
   }

   rownames(CR.out)<- NULL
   return(CR.out)
}

#' Prevalence-Incidence Mixture Models Predictions
#'
#' This function produces cumulative risk predictions given an object of class PIMix
#' a vector of times, and a data set of covariates from which to make predictions.
#' In place of a PIMix object, all elements of the fitted model can be specified instead.
#'
#' @import survival interval
#'
#' @param x object of class PIMix.
#' @param data data set of covariates from which to make predictions.
#' @param time.points numeric vector of times points to produce cumulative risk estimates for.
#' @param model Character string indicating the specific member of the Prevalence-Incidence Mixture Model family to be fitted.  Options are:
#'  \itemize{
#'  \item "semi-parametric" Fits logistic regression and proportional hazards model as the prevalence and incidence models, respectively.
#'  The baseline hazard function is estimated using the iterative convex minorant algorithm. Variance estimates are obtained using bootstrap methods.  Can be computationally expensive;
#'  \item "weakly-parametric" Fits logistic regression and proportional hazards model as the prevalence and incidence models, respectively.
#'   The baseline hazard function is approximated using integrated B-splines.
#'  \item "logistic-Weibull" Fits logistic regression and proportional hazards model as the prevalence and incidence models, respectively.
#'   The baseline hazard function is approximated using a Weibull distribution.
#'  \item "logistic-exponential" Fits logistic regression and proportional hazards model as the prevalence and incidence models, respectively.
#'   The baseline hazard function is approximated using a exponential distribution.
#'  \item "logistic-lognormal" Fits logistic regression and lognormal survival as the prevalence and incidence models, respectively.
#'  \item "logistic-loglogistic" Fits logistic regression and loglogistic survival as the prevalence and incidence models, respectively.
#'  \item "logistic-gengamma" Fits logistic regression and generalized-gamma survival as the prevalence and incidence models, respectively.
#'  \item "logistic-gamma" Fits logistic regression and gamma survival as the prevalence and incidence models, respectively.
#'  }
#' @param prev.coef A vector containing coefficients for the prevalence model.
#' @param incid.coef A vector containing coefficients for the incidence model.
#' @param Lambda.data For semi-parametric or weakly-parametric models, this is a data frame containing the times and baseline cumulative hazard.
#' @param knots For weakly-parametric models, this is a numeric vector of starting time points for each exponential spline.
#' @param spline.para.est For weakly-parametric models, this is a numeric vector of coefficients for each exponential spline.
#' @param cov.mat A matrix containing the covariance matrix for the parameters (not required for semi-parametric models).
#'
#' @return A data frame containing the following columns
#'  \itemize{
#'  \item prev.predictor The design matrix for the prevalence model.
#'  \item incid.predictor The design matrix for the incidence model.
#'  \item time Time.
#'  \item CR Cumulative risk at the time specified.
#'  \item CR.se Standard error for the cumulative risk.
#'  \item LL95 Lower 95 percent confidence limit for the cumulative risk.
#'  \item LL95 Upper 95 percent confidence limit for the cumulative risk.
#'  }
#'
#' @author Li C. Cheung, \email{li.cheung@nih.gov}, Noorie Hyun \email{nhyun@mcw.edu} Xiaoqin Xiong, Qing Pan, Hormuzd A. Katki
#'
#' @references
#' \itemize{
#'  \item Cheung LC, Qing P, Hyun N, Schiffman M, Fetterman B, Castle P, Lorey T, Katki H.
#'        Mixture models for undiagnosed prevalent disease and interval-censored incident disease:
#'        Applications to a cohort assembled from electronic health records. Statistics in Medicine 2017; 
#'        36(22):3583-95.
#'  \item Hyun N, Cheung LC, Pan Q, Katki H.
#'        Flexible risk prediction models for left or interval-censored data from electronic health records.
#'        Annals of Applied Statistics 11(2), 1063-1084.
#' }
#'
#' @export
#'
#' @examples
#'
#'  #PIMixture includes "PIdata" RData file, and PIdata includes the two datasets, PIdata1 and PIdata2
#'  data(PIdata)
#'  model<-"C_CIN3PLUS+L_CIN3PLUS+R_CIN3PLUS~RES_HPV16"
#'  fit1<-PIMixture(p.model=model,data=PIdata1, model="logistic-Weibull")
#'  fit2<-PIMixture(p.model=model,data=PIdata1, model="weakly-parametric",n.knots=5,order=4)
#'  fit3<-PIMixture(p.model=model,data=PIdata1, model="semi-parametric")
#'
#'  model2<-"C_CIN3PLUS+L_CIN3PLUS+R_CIN3PLUS~1"
#'  fit4<-PIMixture(p.model=model2,data=PIdata1, model="non-parametric", conf.int=TRUE)
#'  fit5<-PIMixture(p.model=model2,data=PIdata1, model="semi-parametric", conf.int=TRUE)
#'
#'  test.data<- data.frame(rbind(1,0))
#'  names(test.data)<- "RES_HPV16"
#'  time.points=c(0,12,36,60)
#'  predict1<-PIMixture.predict(x=fit1, data=test.data, time.points=time.points)
#'  predict2<-PIMixture.predict(x=fit2, data=test.data, time.points=time.points)
#'  predict3<-PIMixture.predict(x=fit3, data=test.data, time.points=time.points)
#'  predict4<-PIMixture.predict(x=fit4, data=test.data, time.points=time.points)
#'  predict5<-PIMixture.predict(x=fit5, data=test.data, time.points=time.points)
#'  predict1
#'  predict2
#'  predict3
#'
#'  predict4
#'  predict5
#'
#'  #For stratified random samples
#'  model3<-"C+L+R~X1+X2"
#'  output1<-PIMixture(p.model=model3,data=PIdata2, model="semi-parametric",sample.design=1)
#'  output2<-PIMixture(p.model=model3,data=PIdata2, model="weakly-parametric", n.knots=7,order=4,sample.design=1)
#'  test.data<- data.frame(X1=1,X2=0.5)
#'  time.points<-seq(0,10,by=2)
#'  predict6<-PIMixture.predict(x=output1, data=test.data, time.points=time.points)
#'  predict7<-PIMixture.predict(x=output2, data=test.data, time.points=time.points)
#'  predict6
#'  predict7
PIMixture.predict <- function(x=NULL,data,time.points, model="semi-parametric",
                              prev.coef=NULL,incid.coef=NULL,Lambda.data=NULL,knots=NULL,order=NULL,spline.para.est=NULL,
                              cov.mat=NULL,...){

  requireNamespace("survival", quietly = TRUE)
  requireNamespace("optimx", quietly = TRUE)
  requireNamespace("fdrtool", quietly = TRUE)
  requireNamespace("flexsurv", quietly = TRUE)


  output<- matrix()

  #########################################################
  ## If x is not null, then it must be PIMix object
  ##
  ## if x  is a PIMix object, and x$model is non-parametric,
  ## the calculate the cumulative risks prediction for
  ## non-parametric models
  ##if x is a PIMix object and x$model is weakly-parametric,
  ##then the order is calculated.
  #########################################################

  if (!is.null(x)) {
    if (class(x) != "PIMix") {
      stop("if x is not null, x must be a PIMix object")
    }

    if(x$model == "non-parametric") {
      return(PIMixture.predict.param(x,time.points,mat1=NULL, mat2=NULL))
    }
    if(class(x) == "PIMix" & x$model == "weakly-parametric") {
      order=length(x$exp.spline.coeff)-length(x$knots)+2
    }
  }

  #########################################################
  ## The rest are the cumulative risk predictions for
  ## weakly-parametric, semi-parametric and parametric models
  ##
  ##
  ##
  ## Scenario 1: x is null
  ## prev.coef, incid.coef, mat1, mat2 and model are provided.
  ## Cumulative risks are predicted for semi-parametric,
  ## weakly-parametric and parametric models.
  ##############################################################

  if (!is.null(x) &  class(x) != "PIMix") {
    stop("if x is not null, x must be a PIMix object")
  }

  if(is.null(x)) {
    if(is.null(prev.coef) & is.null(incid.coef)) {
      stop("Either \"prev.coef\" or \"incid.coef\" must be non-null.")
    }

    if(!is.null(prev.coef) & is.null(names(prev.coef))) {
      stop("Must provide the variable names for \"prev.coef\". For intercept, just use \"(Intercept)\"")
    }

    if(!is.null(incid.coef) & is.null(names(incid.coef))) {
      stop("Must provide the variable names for \"incid.coef\". For intercept, just use \"(Intercept)\"")
    }

    prev.predictor<- NULL
    incid.predictor<- NULL

    ###############################################
    # We will set up prev.predictor now, and will
    # set incid.predictor later
    ###############################################
    if((!is.null(prev.coef))  & (class(prev.coef) != "numeric")) {
      stop("prev.coef should be numberic")
    }

    if((!is.null(incid.coef)) & (class(incid.coef) != "numeric")) {
      stop("incid.coef should be numberic")
    }

    mat1<- NULL
    mat2<- NULL
    prev.predictor<- NULL
    incid.predictor<- NULL

    ######################################################
    # mat1: this is the design matrix without intercept
    #
    # Assume the first element prev.coef is the
    # intercept coefficients
    #
    # prev.predictor: it is the design matrix which
    # includes the intercept
    ######################################################

    if(is.null(prev.coef)) {
      mat1<- NULL
      prev.predictor<- NULL
    } else if (length(prev.coef)==1){
      mat1<- NULL
      prev.predictor<- matrix(1)
    } else {
      index.1<- match(names(prev.coef), names(data))
      index.1<- index.1[!is.na(index.1)]

      if(length(index.1) == 0) {
        stop("No variables in \"prev.coef\" are matched with prediction data")
        mat1<- NULL
      } else {
        mat1<- matrix(data[, index.1], nrow=nrow(data))
        colnames(mat1)<- names(data)[index.1]
      }

      prev.predictor<- cbind(1, mat1)
    }


    ######################################################
    # mat2: this is the design matrix without intercept
    #
    # Assume the first element prev.coef is the
    # intercept coefficients
    #
    # incid.predictor: it is the design matrix which
    # does not include the intercept
    #######################################################

    if(tolower(model)%in% c("logistic-weibull",  "logistic-exponential",
                            "logistic-lognormal","logistic-loglogistic",
                            "logistic-gengamma", "logistic-gamma")) {
      if(is.null(cov.mat)) {
        stop("cov.mat must be provided.")
      }

      ######################################################
      # mat2: set up the matrix for incidence model
      # Assume the first element is the intercept coefficients
      # for parametric model
      ######################################################
      if (is.null(incid.coef)) {
        mat2<- NULL
      } else {
        index.2<- match(names(incid.coef), names(data))
        index.2<- index.2[!is.na(index.2)]

        if(length(index.2) == 0) {
          #warning("No variables in \"incid.coef\" are matched with prediction data, so no covariates are adjusted for incidence model")
          mat2<- NULL
        } else  {
          mat2<- matrix(data[, index.2], nrow=nrow(data))
          colnames(mat2)<- names(data)[index.2]
        }
      }

      x<- list()
      Model<- c(rep("logit", length(prev.coef)), rep("surv", length(incid.coef)))
      Label<- paste0("V", (length(prev.coef) + length(incid.coef)))
      Coef<- c(prev.coef, incid.coef)

      x$regression.coef<- data.frame(Model, Label, Coef, stringsAsFactors=FALSE, check.names=FALSE)
      x$model<- model
      x$covariance<- cov.mat
      x$hessian<- -solve(cov.mat)   #negative hessian is used for PIMixture.predict.param


      output<- PIMixture.predict.param(x,time.points,mat1, mat2)

    } else if (tolower(model) == "weakly-parametric") {
      if (is.null(incid.coef)) {
        incid.predictor=NULL
      } else {
        index.2<- match(names(incid.coef), names(data))
        index.2<- index.2[!is.na(index.2)]

        if(length(index.2) == 0) {
          stop("the variable names in \"incid.coef\" does not match the ones in prediction data")
        } else {
          incid.predictor<- matrix(data[, index.2], nrow=nrow(data))
          colnames(incid.predictor)<- names(data)[index.2]
        }
      }

      output<- cumulative.risk(logit.coef=prev.coef,
                               cox.coef=incid.coef,
                               logit.predictor=prev.predictor,
                               cox.predictor=incid.predictor,
                               Lambda.data=Lambda.data,
                               CR.time.points=time.points,
                               knots=knots,
                               order=order,
                               spline.para.est=spline.para.est,
                               cov.mat=cov.mat)

    } else if (tolower(model) == "semi-parametric") {

      if (is.null(incid.coef)) {
        incid.predictor=NULL
      } else {
        index.2<- match(names(incid.coef), names(data))
        index.2<- index.2[!is.na(index.2)]

        if(length(index.2) == 0) {
          stop("the variable names in \"incid.coef\" does not match the ones in prediction data")
        } else {
          incid.predictor<- matrix(data[, index.2], nrow=nrow(data))
          colnames(incid.predictor)<- names(data)[index.2]
        }
      }

      output<- cumulative.risk(logit.coef=prev.coef,
                               cox.coef=incid.coef,
                               logit.predictor=prev.predictor,
                               cox.predictor=incid.predictor,
                               Lambda.data=Lambda.data,
                               CR.time.points=time.points)
    } else {
      stop("wrong model option")
    }
  }  #end of if(is.null(x))

  ##############################################################
  ## Scenario 2: x is a PIMix object
  ##
  ## cumulative risks are predicted for
  ## semi-parametric, weakly-parametric and parametric models.
  ##############################################################
  if (class(x) == "PIMix") {

    ##############################################################
    # First get the design matrix for logistic and survival model
    ##############################################################

    p.model<- x$p.model
    i.model<- x$i.model

    bb.p<- strsplit(p.model, "~")
    p.model<-paste0("~",bb.p[[1]][2])
    fml<-as.formula(p.model)
    p.vars<- all.vars(fml)
    p.mat<- get.design.mat(data, p.model)

    bb.i<- strsplit(i.model, "~")
    i.model<-paste0("~", bb.i[[1]][2])
    fml<-as.formula(i.model)
    i.vars<- all.vars(fml)
    i.mat<- get.design.mat(data, i.model)

    if(sum(p.vars%in%names(data)) != length(p.vars)) {
      stop(paste0("For p.model, the predicted dataset does not contain variables ",
                  p.vars[!p.vars%in%names(data)]))
    }

    if( sum(i.vars%in%names(data)) != length(i.vars)) {
      stop(paste0("For i.model, the predicted dataset does not contain variables ",
          i.vars[!i.vars%in%names(data)]))
    }

    if(tolower(x$model)%in% c("logistic-weibull",  "logistic-exponential",
                              "logistic-lognormal","logistic-loglogistic",
                              "logistic-gengamma", "logistic-gamma")) {

      if(is.null(p.mat) | (!is.null(p.mat) & ncol(p.mat) ==1)) {
        mat1<- NULL
      } else {
        mat1<- matrix(p.mat[, -1], nrow=nrow(p.mat), ncol=(ncol(p.mat)-1))
        colnames(mat1)<- colnames(p.mat)[-1]
      }

      if(is.null(i.mat) | (!is.null(i.mat) & ncol(i.mat) ==1)) {
        mat2<- NULL
      } else {
        mat2<- matrix(i.mat[, -1], nrow=nrow(i.mat), ncol=(ncol(i.mat)-1))
        colnames(mat2)<- colnames(i.mat)[-1]
      }

      output<- PIMixture.predict.param(x,time.points,mat1, mat2)

    } else if (tolower(x$model)%in% c("weakly-parametric", "semi-parametric")) {
      prev.coef<- x$regression.coef[x$regression.coef[,1]== "logit", 3]
      incid.coef<- x$regression.coef[tolower(x$regression.coef[,1])== "cox",   3]

      if(is.null(p.mat)) {
        prev.predictor<- NULL
      } else if (ncol(p.mat) ==1) {
        prev.predictor<- as.matrix(1)
      } else {
        prev.predictor<- p.mat
      }

      if(is.null(i.mat) | (!is.null(i.mat) & ncol(i.mat) ==1)) {
        incid.predictor<- NULL
      } else {
        incid.predictor<- matrix(i.mat[, -1], nrow=nrow(i.mat), ncol=(ncol(i.mat)-1))
        colnames(incid.predictor)<- colnames(i.mat)[-1]
      }

      Lambda.data<- x$cum.hazard

      if(tolower(x$model) == "weakly-parametric") {
        knots<- x$knots
        spline.para.est<- x$exp.spline.coeff
        if(!is.null(spline.para.est)) {
         spline.para.est= log(spline.para.est)
        }

        cov.mat=x$covariance
        output<- cumulative.risk(logit.coef=prev.coef,
                                 cox.coef=incid.coef,
                                 logit.predictor=prev.predictor,
                                 cox.predictor=incid.predictor,
                                 Lambda.data=Lambda.data,
                                 CR.time.points=time.points,
                                 knots=knots,
                                 order=order,
                                 spline.para.est=spline.para.est,
                                 cov.mat=cov.mat)

      } else if (tolower(x$model) == "semi-parametric") {
        output<- cumulative.risk(logit.coef=prev.coef,
                                 cox.coef=incid.coef,
                                 logit.predictor=prev.predictor,
                                 cox.predictor=incid.predictor,
                                 Lambda.data=Lambda.data,
                                 CR.time.points=time.points)
      } else {
        stop("x$model is not correct")
      }
    }
  }

  return(output)
}
