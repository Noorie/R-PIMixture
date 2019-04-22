Ispline <-function(x,order,knots){
    # I spline with degree=order, 1 for linear, 2 for quadratic, 3 for cubic, etc.
    # x is a row vector
    # knots are a sequence of increasing points containing the two end points of the interval (a, b), where the curve is to be estimated.
    # The output is a matrix with dimension (length(knots)+order-2) by length(x).

    k=order+1
    m=length(knots)
    n=m-2+k # number of parameters
    t=c(rep(1,k)*knots[1], knots[2:(m-1)], rep(1,k)*knots[m]) # newknots

    yy1=array(rep(0,(n+k-1)*length(x)),dim=c(n+k-1, length(x)))
    for (l in k:n){
      yy1[l,]=(x>=t[l] & x<t[l+1])/(t[l+1]-t[l])
    }

    yytem1=yy1
    for (ii in 1:order){
      yytem2=array(rep(0,(n+k-1-ii)*length(x)),dim=c(n+k-1-ii, length(x)))
      for (i in (k-ii):n){
        yytem2[i,]=(ii+1)*((x-t[i])*yytem1[i,]+(t[i+ii+1]-x)*yytem1[i+1,])/(t[i+ii+1]-t[i])/ii
      }
      yytem1=yytem2
    }

    index=rep(0,length(x))
    for (i in 1:length(x)){
      index[i]=sum(t<=x[i])
    }

    yy=array(rep(0,(n-1)*length(x)),dim=c(n-1,length(x)))

    if (order==1){
      for (i in 2:n){
        yy[i-1,]=(i<index-order+1)+(i==index)*(t[i+order+1]-t[i])*yytem2[i,]/(order+1)
      }
    }else{
      for (j in 1:length(x)){
        for (i in 2:n){
          if (i<(index[j]-order+1)){
            yy[i-1,j]=1
          }else if ((i<=index[j]) && (i>=(index[j]-order+1))){
            yy[i-1,j]=(t[(i+order+1):(index[j]+order+1)]-t[i:index[j]])%*%yytem2[i:index[j],j]/(order+1)
          }else{
            yy[i-1,j]=0
          }
        }
      }
    }
    return(yy)
}

score<-function(para,samp.data,n.beta,n.gamma,design.mat,xmat,n.knots,order2){


  beta<-para[1:n.beta]

  mm1<-n.beta+1
  mm2<-n.beta+n.gamma


  mm3<- n.beta+n.gamma+1
  mm4<-length(para)

  beta.risk<-as.numeric(design.mat%*%beta)
  spline.para<-para[mm3:mm4]

  if(n.gamma==0){
    gamma<-0
    gamma.risk<-0
  }else if(n.gamma>0){
    gamma<-para[mm1:mm2]
    gamma.risk<-as.numeric(xmat%*%gamma)
  }

  n.samp2<-nrow(samp.data)

  ###Spline-based Lambda Calculation #############


  time.list<-sort(unique(c(samp.data$L[samp.data$L!=-999],samp.data$R[samp.data$R!=-999&samp.data$R<Inf])))
  by.knot<-1/(n.knots-1)
  knots2<-as.numeric(quantile(time.list,seq(0,1,by=by.knot)))

  basis.L<-Ispline(as.numeric(samp.data$L), order = order2, knots = knots2)
  basis.R<-Ispline(samp.data$R, order = order2, knots = knots2)

  basis.L[,samp.data$L==-999|samp.data$L==0]<-0
  basis.R[,samp.data$R==-999]<-0
  basis.R[,samp.data$R==Inf]<-Inf


  Lambda.L.mat<-exp(spline.para)*basis.L
  Lambda.R.mat<-exp(spline.para)*basis.R

  LambdaL<-as.numeric(colSums(Lambda.L.mat))
  LambdaR<-as.numeric(colSums(Lambda.R.mat))

  ####################################################################################
  K0<-as.numeric(samp.data$K==0)
  K1<-as.numeric(samp.data$K==1)
  C0<-as.numeric(samp.data$C==0)

  logit.risk1<-exp(beta.risk)/(1+exp(beta.risk))
  logit.risk2<-exp(beta.risk)/(1+exp(beta.risk))^2
  S.R<-exp(-LambdaR*exp(gamma.risk))
  S.L<-exp(-LambdaL*exp(gamma.risk))
  Lam.RX<-LambdaR*exp(gamma.risk)
  Lam.LX<-LambdaL*exp(gamma.risk)


  beta.score1<-as.numeric(K1*samp.data$samp.weight*(samp.data$C-logit.risk1))*design.mat
  beta.score2<-as.numeric(logit.risk1*S.R/(1+exp(beta.risk)-S.R))*K0*samp.data$samp.weight*design.mat
  beta.score2[K0==0,]<-0
  beta.score2[samp.data$R==Inf,]<-0

  beta.score<-as.numeric(colSums(beta.score1+beta.score2))
  beta.score.vec<-(beta.score1+beta.score2)

  comp1<-(S.R*Lam.RX-S.L*Lam.LX)/(S.L-S.R)       #when 0<L<R<Inf
  comp2<-(-S.L*Lam.LX)/S.L                       #When 0<L<R=Inf
  comp22<- (S.R*Lam.RX)/(1-S.R)                  #When 0=L<R<Inf

  comp1[samp.data$L>0&samp.data$L!=-999&samp.data$R==Inf]<-comp2[samp.data$L>0&samp.data$L!=-999&samp.data$R==Inf]
  comp1[samp.data$L==0&samp.data$R<Inf&samp.data$L!=-999]<-comp22[samp.data$L==0&samp.data$R<Inf&samp.data$L!=-999]
  comp1[samp.data$L==0&samp.data$R==Inf]<-0

  comp1[samp.data$C==1]<-0   #When C==1 then L=R=-999, that is, NA

  comp3<-(S.R*Lam.RX)/(1+exp(beta.risk)-S.R)
  comp3[K0==0]<-0
  comp3[samp.data$C==1]<-0
  comp3[samp.data$R==Inf] <-0

  if(n.gamma>0){
    gamma.score1<-as.numeric(K1*samp.data$samp.weight*C0*comp1)*xmat
    gamma.score2<-as.numeric(K0*samp.data$samp.weight*comp3)*xmat

    gamma.score<-as.numeric(colSums((gamma.score1+gamma.score2)))
    gamma.score.vec<-(gamma.score1+gamma.score2)
  }


  comp1.L<-(-S.L*exp(gamma.risk))/(S.L-S.R)       #when 0<L<R<Inf
  comp1.R<-(S.R*exp(gamma.risk))/(S.L-S.R)       #when 0<L<R<Inf


  comp1.L[samp.data$L==0]<-0
  comp1.R[samp.data$R==Inf]<-0

  comp1.L[samp.data$C==1]<-0   #When C==1 then L=R=-999, that is, NA
  comp1.R[samp.data$C==1]<-0   #When C==1 then L=R=-999, that is, NA

  comp33<-(S.R*exp(gamma.risk))/(1+exp(beta.risk)-S.R)
  comp33[K0==0]<-0
  comp33[samp.data$C==1]<-0
  comp33[samp.data$R==Inf] <-0


  b.score11<-matrix(as.numeric(K1*C0*samp.data$samp.weight*comp1.L),dim(basis.L)[1],n.samp2,1) *Lambda.L.mat
  b.score12<-matrix(as.numeric(K1*C0*samp.data$samp.weight*comp1.R),dim(basis.L)[1],n.samp2,1) *Lambda.R.mat
  b.score2<-matrix(as.numeric(K0*comp33*samp.data$samp.weight),dim(basis.L)[1],n.samp2,1)*Lambda.R.mat

  b.score12[Lambda.R.mat==Inf]<-0
  b.score2[Lambda.R.mat==Inf]<-0

  b.score<-rowSums(b.score11+b.score12+b.score2)
  b.score.vec<- (b.score11+b.score12+b.score2)

  out<-list()

  #log-likelihood
  if(n.gamma>0){
    out$score<-c(beta.score,gamma.score,b.score)
    out$score.vec<-cbind(beta.score.vec ,gamma.score.vec,t(b.score.vec))
    out$score.vec.nowt<-cbind(beta.score.vec ,gamma.score.vec,t(b.score.vec))/samp.data$samp.weight
  } else if(n.gamma==0){
    out$score<-c(beta.score,b.score)
    out$score.vec<-cbind(beta.score.vec ,t(b.score.vec))
    out$score.vec.nowt<-cbind(beta.score.vec ,t(b.score.vec))/samp.data$samp.weight
  }

  return(out)
}#the end of the function score

hessian<-function(para,samp.data,n.beta,n.gamma,design.mat,xmat,n.knots,order2){


  beta<-para[1:n.beta]

  mm1<-n.beta+1
  mm2<-n.beta+n.gamma


  mm3<- n.beta+n.gamma+1
  mm4<-length(para)

  beta.risk<-as.numeric(design.mat%*%beta)
  spline.para<-para[mm3:mm4]

  if(n.gamma==0){
    gamma<-0
    gamma.risk<-0
  }else if(n.gamma>0){
    gamma<-para[mm1:mm2]
    gamma.risk<-as.numeric(xmat%*%gamma)
  }


  n.samp2<-nrow(samp.data)

  ###Spline-based Lambda Calculation #############

  time.list<-sort(unique(c(samp.data$L[samp.data$L!=-999],samp.data$R[samp.data$R!=-999&samp.data$R<Inf])))
  by.knot<-1/(n.knots-1)
  knots2<-as.numeric(quantile(time.list,seq(0,1,by=by.knot)))

  basis.L<-Ispline(as.numeric(samp.data$L), order = order2, knots = knots2)
  basis.R<-Ispline(samp.data$R, order = order2, knots = knots2)


  basis.L[,samp.data$L==-999|samp.data$L==0]<-0
  basis.R[,samp.data$R==-999]<-0
  basis.R[,samp.data$R==Inf]<-Inf


  Lambda.L.mat<-exp(spline.para)*basis.L
  Lambda.R.mat<-exp(spline.para)*basis.R

  LambdaL<-as.numeric(colSums(Lambda.L.mat))
  LambdaR<-as.numeric(colSums(Lambda.R.mat))


  ####################################################################################
  K0<-as.numeric(samp.data$K==0)
  K1<-as.numeric(samp.data$K==1)
  C0<-as.numeric(samp.data$C==0)
  IPW<-samp.data$samp.weight

  logit.risk1<-exp(beta.risk)/(1+exp(beta.risk))
  logit.risk2<-exp(beta.risk)/(1+exp(beta.risk))^2
  logit.risk3<-exp(2*beta.risk)/(1+exp(beta.risk))^3
  S.R<-exp(-LambdaR*exp(gamma.risk))
  S.L<-exp(-LambdaL*exp(gamma.risk))
  Lam.RX<-LambdaR*exp(gamma.risk)
  Lam.LX<-LambdaL*exp(gamma.risk)

  ######### derivative with respect to beta, beta

  beta.curv1<-as.numeric(K1*IPW*(-logit.risk2))*design.mat
  beta.curv2<-(S.R*logit.risk1*(1-exp(beta.risk)-S.R/(1+exp(beta.risk))))/(1+exp(beta.risk)-S.R)^2
  beta.curv2[K1==1]<-0              #this component is for K0 and When K=1, L=R=-999, that is, NA
  beta.curv22<-as.numeric(beta.curv2*K0*IPW)*design.mat
  beta.curv3<- (beta.curv22+beta.curv1)

  beta.curv<- as.matrix(t(design.mat)%*%beta.curv3)


  ######### derivative with respect to  beta, b
  beta.b1<- -as.numeric(S.R*IPW*exp(gamma.risk)*exp(beta.risk)/(1+exp(beta.risk)-S.R)^2)
  beta.b1[K0==0]<-0
  beta.b1[samp.data$R==Inf]<-0
  beta.b2<- (beta.b1*K0)*t(Lambda.R.mat)
  beta.b2[samp.data$R==Inf,]<-0

  beta.b<-t(design.mat)%*%beta.b2

  if(n.gamma>0){
    ######### derivative with respect to beta, gamma

    beta.gamma1<- -as.numeric(S.R*IPW*Lam.RX*exp(beta.risk)/(1+exp(beta.risk)-S.R)^2)
    beta.gamma1[K0==0]<-0
    beta.gamma1[samp.data$R==Inf]<-0

    beta.gamma<-t(design.mat)%*%(beta.gamma1*xmat)

    ######### derivative with respect to gamma, gamma
    gamma.curv1<-as.numeric(K1*IPW*C0*(S.L*Lam.LX*(Lam.LX-1)+S.R*Lam.RX*(1-Lam.RX) )/(S.L-S.R))
    gamma.curv1[samp.data$C==1|K0==1]<-0
    gamma.curv1[samp.data$R==Inf&samp.data$L==0]<-0

    gamma.curv11<-as.numeric(K1*IPW*C0*(S.L*Lam.LX*(Lam.LX-1) )/(S.L))
    gamma.curv1[samp.data$R==Inf&samp.data$L>0]<-gamma.curv11[samp.data$R==Inf&samp.data$L>0]

    gamma.curv2<- -as.numeric(K1*IPW*C0*((S.R*Lam.RX-S.L*Lam.LX)/(S.L-S.R))^2)
    gamma.curv2[samp.data$C==1|K0==1]<-0
    gamma.curv2[samp.data$R==Inf&samp.data$L==0]<-0

    gamma.curv21<- -as.numeric(K1*IPW*C0*((-S.L*Lam.LX)/(S.L))^2)
    gamma.curv2[samp.data$R==Inf&samp.data$L>0]<-gamma.curv21[samp.data$R==Inf&samp.data$L>0]


    gamma.curv3<- as.numeric(K0*(S.R*IPW*Lam.RX*(1-Lam.RX)/(1+exp(beta.risk)-S.R)))
    gamma.curv3[samp.data$C==1]<-0
    gamma.curv3[samp.data$R==Inf]<-0

    gamma.curv4<- -as.numeric(K0*IPW*((S.R*Lam.RX)/(1+exp(beta.risk)-S.R))^2)
    gamma.curv4[samp.data$C==1]<-0
    gamma.curv4[samp.data$R==Inf]<-0


    gamma.curv.K1<-(gamma.curv1+gamma.curv2)*xmat
    gamma.curv.K0<-(gamma.curv3+gamma.curv4)*xmat

    gamma.curv<-t(xmat)%*%(gamma.curv.K0+gamma.curv.K1)

    ######### derivative with respect to  gamma, b
    #0<L<R<Inf
    gamma.b.L1<-as.numeric(K1*IPW*C0*(S.L*exp(gamma.risk))*(Lam.LX-1)/(S.L-S.R))
    gamma.b.R1<-as.numeric(K1*IPW*C0*(S.R*exp(gamma.risk))*(1-Lam.RX)/(S.L-S.R))


    gamma.b.L1[samp.data$L==0]<-0
    gamma.b.R1[samp.data$R==Inf]<-0
    gamma.b.L1[samp.data$C==1]<-0
    gamma.b.R1[samp.data$C==1]<-0


    gamma.b.R1.K1<-gamma.b.R1*t(Lambda.R.mat)

    gamma.b.R1.K1[t(Lambda.R.mat)==Inf]<-0

    gamma.b.K1.1<-t(xmat)%*%(gamma.b.L1*t(Lambda.L.mat)+gamma.b.R1.K1)


    gamma.b1<-as.numeric(-1*K1*C0*IPW*(S.R*Lam.RX-S.L*Lam.LX)/(S.L-S.R)^2)
    gamma.b2.R<-as.numeric(S.R*exp(gamma.risk))*t(Lambda.R.mat)
    gamma.b2.L<- -as.numeric(S.L*exp(gamma.risk))*t(Lambda.L.mat)


    gamma.b1[samp.data$C==1]<-0
    gamma.b1[samp.data$L==0&samp.data$R==Inf]<-0


    #0<L<R=Inf
    gamma.b11<-as.numeric(K1*IPW*C0*S.L*Lam.LX/(S.L)^2)
    gamma.b1[samp.data$L>0&samp.data$R==Inf]<-gamma.b11[samp.data$L>0&samp.data$R==Inf]


    gamma.b2.R[samp.data$R==Inf,]<-0
    gamma.b2.L[samp.data$C==1,]<-0
    gamma.b2.R[samp.data$C==1,]<-0

    gamma.b.K1.2<-t(gamma.b1*xmat)%*% (gamma.b2.L +gamma.b2.R)


    gamma.b3<- as.numeric(K0*IPW*S.R*exp(gamma.risk)*( 1 -Lam.RX)/(1+exp(beta.risk)-S.R) )
    gamma.b3[samp.data$R==Inf]<-0
    gamma.b3[samp.data$C==1]<-0

    gamma.b33<-gamma.b3*t(Lambda.R.mat)
    gamma.b33[samp.data$R==Inf,]<-0

    gamma.b.K0.1<- t(xmat)%*%(gamma.b33)

    gamma.b4<- -as.numeric(K0*IPW*S.R^2*Lam.RX*exp(gamma.risk)/(1+exp(beta.risk)-S.R)^2)
    gamma.b4[samp.data$R==Inf]<-0
    gamma.b4[samp.data$C==1]<-0

    gamma.b44<-gamma.b4*t(Lambda.R.mat)
    gamma.b44[samp.data$R==Inf,]<-0

    gamma.b.K0.2<- t(xmat)%*%(gamma.b44)

    gamma.b<-gamma.b.K1.1+gamma.b.K1.2+gamma.b.K0.1+gamma.b.K0.2
  }


  ######### derivative with respect to  b, b
  #0<L<R<Inf
  Lambda.R.mat2<-Lambda.R.mat
  Lambda.R.mat2[Lambda.R.mat==Inf]<-0

  b.b.L1<-as.numeric(K1*C0*IPW*(S.L*exp(gamma.risk))/(S.L-S.R))
  b.b.R1<- -as.numeric(K1*C0*IPW*(S.R*exp(gamma.risk))/(S.L-S.R))

  b.b.L1[samp.data$C==1]<-0
  b.b.R1[samp.data$C==1]<-0
  b.b.R1[samp.data$R==Inf]<-0
  b.b.L1[samp.data$L==0]<-0


  b.b.K1.11<-Lambda.L.mat%*%(b.b.L1*as.numeric(exp(gamma.risk))*t(Lambda.L.mat)) + Lambda.R.mat2%*%(b.b.R1*as.numeric(exp(gamma.risk))*t(Lambda.R.mat2))
  b.b.K1.12<- diag(colSums(-b.b.L1*t(Lambda.L.mat)) + colSums(-b.b.R1*t(Lambda.R.mat2)))
  b.b.K1.1<-b.b.K1.11+b.b.K1.12


  b.b.R2<-  as.numeric(K1*C0*IPW*(S.R*exp(gamma.risk))/(S.L-S.R))
  b.b.L2<-  -as.numeric(K1*C0*IPW*(S.L*exp(gamma.risk))/(S.L-S.R))
  b.b.R2[samp.data$C==1]<-0
  b.b.L2[samp.data$C==1]<-0
  b.b.R2[samp.data$R==Inf]<-0
  b.b.L2[samp.data$L==0]<-0

  b.b.2<-b.b.R2*t(Lambda.R.mat2)+b.b.L2*t(Lambda.L.mat)

  b.b.R22<-  as.numeric(K1*C0*(S.R*exp(gamma.risk))/(S.L-S.R))
  b.b.L22<-  -as.numeric(K1*C0*(S.L*exp(gamma.risk))/(S.L-S.R))
  b.b.R22[samp.data$C==1]<-0
  b.b.L22[samp.data$C==1]<-0
  b.b.R22[samp.data$R==Inf]<-0
  b.b.L22[samp.data$L==0]<-0

  b.b.22<-b.b.R22*t(Lambda.R.mat2)+b.b.L22*t(Lambda.L.mat)

  b.b.K1.2<-t(-b.b.2)%*%(b.b.22)

  b.b.3<- as.numeric(K0*IPW*S.R*exp(gamma.risk)/(1+exp(beta.risk)-S.R))
  b.b.3[samp.data$C==1]<-0
  b.b.3[samp.data$R==Inf]<-0

  b.b.K0.11<- -Lambda.R.mat2%*%(b.b.3*as.numeric(exp(gamma.risk))*t(Lambda.R.mat2))
  b.b.K0.12<- diag(colSums(b.b.3*t(Lambda.R.mat2)))

  b.b.K0.1<-b.b.K0.11+b.b.K0.12

  b.b.4<- -as.numeric(K0*IPW*S.R*exp(gamma.risk)/(1+exp(beta.risk)-S.R))
  b.b.4[samp.data$R==Inf]<-0
  b.b.4[samp.data$C==1]<-0

  b.b.44<- as.numeric(K0*S.R*exp(gamma.risk)/(1+exp(beta.risk)-S.R))
  b.b.44[samp.data$R==Inf]<-0
  b.b.44[samp.data$C==1]<-0

  b.b.K0.2<-t(b.b.4*t(Lambda.R.mat2))%*%(b.b.44*t(Lambda.R.mat2))

  b.curv<-b.b.K1.1+b.b.K1.2+b.b.K0.1+b.b.K0.2

  if(n.gamma>0){
    Hess<-hessian<-rbind(cbind(beta.curv,beta.gamma,beta.b),
                         cbind(t(beta.gamma),gamma.curv,gamma.b),
                         cbind(t(beta.b),t(gamma.b),b.curv))
  } else if(n.gamma==0){
    Hess<-hessian<-rbind(cbind(beta.curv,beta.b),
                         cbind(t(beta.b),b.curv))
  }

  colnames(Hess)<-NULL
  #hessian matrix
  return(Hess)
}#the end of the function hessian




## for Phase 1 - study (based on linearization) ###############
cov.mat.ph1<-function(delta.theta,samp.data){
  n.para<-ncol(delta.theta)
  cov.mat<-0*diag(n.para)

  aaa<-which(table(samp.data$strata)!=1)
  strata.list<-as.numeric(names(aaa))

  for (k in strata.list){
    average<-apply(delta.theta[samp.data$strata==k,],2,mean)
    delta.theta.gp<-delta.theta[samp.data$strata==k,]
    n.gp<-nrow(delta.theta.gp)
    ave.mat<-matrix(average,n.gp,n.para,1)
    c.delta.theta.gp<-delta.theta.gp-ave.mat
    cov.mat<-cov.mat+n.gp/(n.gp-1)*(t(c.delta.theta.gp)%*%c.delta.theta.gp)
  }
  return(cov.mat)
}#the end of the function cov.mat.ph1


############### Phase 2 Design - Covariance ###############################################
cov.mat.ph2<-function(score.nowt,samp.data){
  n.para<-ncol(score.nowt)
  cov.mat<-0*diag(n.para)

  aaa<-which(table(samp.data$strata)!=1)
  strata.list<-as.numeric(names(aaa))

  for (k in strata.list){

    score.by.strata<-score.nowt[samp.data$strata==k,]
    n.gp<-nrow(score.by.strata)
    strata.fraction<-as.numeric(samp.data$strata.frac[samp.data$strata==k])[1]
    samp.data$s.frc<-1/samp.data$samp.weight
    samp.fraction<-as.numeric(samp.data$s.frc[samp.data$strata==k])[1]
    A<- (1/n.gp)*(t(score.by.strata)%*%score.by.strata)
    B<-colMeans(score.by.strata)
    B2<-B%*%t(B)
    cond.var<-A-B2
    cov.mat<-cov.mat+(strata.fraction*(1-samp.fraction)/samp.fraction)*(cond.var)
  }
  return(cov.mat)
}#the end of the function cov.mat.ph2



#################################################################################################



ipw.lc.splines<- function(samp.data, n.beta, n.gamma, design.mat, xmat,
                          sample.design,N,n.knots,order, max.time, reg.initials,
                          convergence.criteria,iteration.limit,
                          time.interval,logit.predictor,cox.predictor,CR.time.points, ...){

  output<- list()



  ################ FUNCTIONS: observed likelihood ###################################
  obs.like<-function(para,samp.data,n.beta,n.gamma,design.mat2,xmat2,n.knots,order){


    beta<-para[1:n.beta]

    mm1<-n.beta+1
    mm2<-n.beta+n.gamma


    mm3<- n.beta+n.gamma+1
    mm4<-length(para)

    if(n.gamma==0){
      gamma<-0
      gamma.risk<-0
    }else if(n.gamma>0){
      gamma<-para[mm1:mm2]
      gamma.risk<-as.numeric(xmat2%*%gamma)
    }

    spline.para<-para[mm3:mm4]
    beta.risk<-as.numeric(design.mat2%*%beta)


    ###Spline-based Lambda Calculation #############


    time.list<-sort(unique(c(samp.data$L[samp.data$L!=-999],samp.data$R[samp.data$R!=-999&samp.data$R<Inf])))
    by.knot<-1/(n.knots-1)
    knots2<-as.numeric(quantile(time.list,seq(0,1,by=by.knot)))

    basis.L<-Ispline(as.numeric(samp.data$L), order = order, knots = knots2)
    basis.R<-Ispline(samp.data$R, order = order, knots = knots2)


    basis.L[,samp.data$L==-999|samp.data$L==0]<-0
    basis.R[,samp.data$R==-999]<-0
    basis.R[,samp.data$R==Inf]<-Inf


    Lambda.L.mat<-exp(spline.para)*basis.L
    Lambda.R.mat<-exp(spline.para)*basis.R

    LambdaL<-as.numeric(colSums(Lambda.L.mat))
    LambdaR<-as.numeric(colSums(Lambda.R.mat))



    B<-log(    exp(-LambdaL*exp(gamma.risk))-exp(-LambdaR*exp(gamma.risk)) )
    B[samp.data$C==1]<-0                                               #When C==1, L=R=999
    B[abs(LambdaL-LambdaR)<1e-10]<-0               #to avoid the value of Inf/-Inf when 1/B and B goes to 0


    B2<-1-exp(-LambdaR*exp(gamma.risk))                   #K=0
    B2[samp.data$K==1]<-0

    ll_K1<-as.numeric(samp.data$K==1)*samp.data$samp.weight*(samp.data$C*beta.risk-log(1+exp(beta.risk)) + B)
    ll_K2<-as.numeric(samp.data$K==0)*log(exp(beta.risk)/(1+exp(beta.risk))+B2/(1+exp(beta.risk)) )*samp.data$samp.weight

    #log-likelihood
    return(sum(ll_K1+ll_K2))
  }

  opt1 <- function(para2){
    out<-obs.like(para2,samp.data,n.beta,n.gamma,design.mat2=design.mat,xmat2=xmat,n.knots,order)
    return(-out)
  }

  gr.opt1 <- function(para2){
    out<-score(para2,samp.data,n.beta,n.gamma,design.mat,xmat,n.knots,order)
    return(-out$score)
  }

  ########### Initial parameters  ################################################

  n.regpara<-n.beta+n.gamma
  n.splines=n.knots+order-2

  n.para1<-n.regpara+1
  n.para2<-n.regpara+n.splines


  if(is.null(reg.initials)){
    if(n.gamma==0){
      initials<-rep(0.05,n.beta)
    } else if (n.gamma>0){
      initials<-c(rep(0.05,n.beta),rep(0,n.gamma))
    }
  } else {
    if(length(reg.initials)!=n.regpara){
      stop("The number of regression parameter initials are not matched with the model.")
    } else {
      initials<-reg.initials
    }
  }


  n.knots22<-n.knots+order-2 -4                #n.splines=n.knots+order-2

  if (n.knots22>0){
    ini.spline.para<-c(rep(-3,4),rep(-2,n.knots22))
  } else if(n.knots22==0){
    ini.spline.para<-rep(-3,4)
  } else if(n.knots22<0){
    n.knots222<-n.knots+order-2
    ini.spline.para<-rep(-3,n.knots222)
  }

  time.list<-sort(unique(c(samp.data$L[samp.data$L!=-999],samp.data$R[samp.data$R!=-999&samp.data$R<Inf])))
  by.knot<-1/(n.knots-1)
  knots2<-as.numeric(quantile(time.list,seq(0,1,by=by.knot)))

  if(missing(max.time)){
    time.points2<-seq(0,max(time.list),by=time.interval)
  }else if(!missing(max.time)){
    time.points2<-seq(0,max.time,by=time.interval)
  }


  old.para<-c(initials,ini.spline.para)

  keep.looping<-1
  ptm<-proc.time()
  iteration<-0

  while (keep.looping==1){
    iteration<-iteration+1
    res1<-optim(old.para,gr=gr.opt1,opt1,method='L-BFGS-B')
    new.para<-res1$par


    diff1<-max(abs(old.para-new.para))

    old.para<-new.para


    conv.marker<-diff1


    if(conv.marker<convergence.criteria|iteration>iteration.limit){
      keep.looping<-0
      para.est<-old.para
      spline.para.est<-para.est[n.para1:n.para2]
    }
  }
  run.time<-((proc.time()-ptm)[3])/60
  remove(diff1)


  basis.Ispline<-Ispline(time.points2,order=order,knots=knots2)
  Lambda.est<-colSums(exp(spline.para.est)*basis.Ispline)

  ## Prevalence/Cumulative Risk Estimation #########################################################
  if(!missing(logit.predictor)&!missing(CR.time.points)){
    logit.risk.est<-as.numeric(logit.predictor%*%para.est[1:n.beta])
    prevalence.est<-as.numeric(exp(logit.risk.est)/(1+exp(logit.risk.est)))
  }

  if(n.gamma>0 & !missing(cox.predictor)&!missing(CR.time.points)){
    cox.risk.est<- as.numeric(cox.predictor%*%para.est[c(1+n.beta):n.regpara])
  } else if(n.gamma==0&!missing(CR.time.points)&!missing(logit.predictor)){
    cox.risk.est<- 0
  }

  if(exists("prevalence.est")==TRUE&exists("cox.risk.est")==TRUE){
    cum.risk.est.full<-cumulative.risk2(prevalence.est, Lambda.est ,cox.risk.est)
    if(!missing(CR.time.points)){
      cum.risk.est<-cumulative.risk2(prevalence.est, Lambda.est[which(time.points2 %in% CR.time.points)] ,cox.risk.est)
    }
  }


  ###################### Variance Calculatoin ######################################################
  H<-hessian(para.est,samp.data,n.beta,n.gamma,design.mat,xmat,n.knots,order)
  det.H<-abs(det(H))

  score.regpara<-score(para.est,samp.data,n.beta,n.gamma,design.mat,xmat,n.knots,order)

  ## when no sampling design ##############
  n.samp<- nrow(samp.data)

  if ((sum(samp.data$samp.weight==1)==n.samp)&(det.H>0)&is.null(sample.design)){
    inv.H<-solve(H,tol=1e-300)
    cov.mat<- -inv.H
    theta.se<-sqrt(diag(cov.mat))

    para.95UL<-as.numeric(para.est[1:n.regpara]+1.96*theta.se[1:n.regpara])
    para.95LL<-as.numeric(para.est[1:n.regpara]-1.96*theta.se[1:n.regpara])
  }

  ## when Phase 1 design ##############
  if(!is.null(sample.design)){
    if(det.H>0 & sample.design==1){
      inv.H<-solve(H,tol=1e-300)

      delta.s<-score.regpara$score.vec

      cov.mat.theta<-cov.mat.ph1(delta.s,samp.data)

      cov.mat<-inv.H%*%cov.mat.theta%*%inv.H
      theta.se<-sqrt(diag(cov.mat))

      para.95UL<-as.numeric(para.est[1:n.regpara]+1.96*theta.se[1:n.regpara])
      para.95LL<-as.numeric(para.est[1:n.regpara]-1.96*theta.se[1:n.regpara])
    }
  }

  ## when phase 2 design ##############
  if(!is.null(sample.design)){
    if((sample.design==2)&(c("strata.frac")%in% names(samp.data))&!missing(N)){

      score.nowt<-score.regpara$score.vec.nowt
      score.wt<-score.regpara$score.vec

      info.mat<- (t(score.wt)%*%score.nowt)
      det.info<-abs(det(info.mat))

      if(det.info>0){

        inv.info<-solve(info.mat,tol=1e-300)

        inv.info.mat<- inv.info*(sum(samp.data$samp.weight)) #because of subgroup analysis, sum of sampling weights replaces N

        cov.mat.phase1<-inv.info.mat/N

        cov.mat.by.strata<-cov.mat.ph2(score.nowt,samp.data)

        cov.mat.phase2<-(inv.info.mat%*%cov.mat.by.strata%*%inv.info.mat)/N

        cov.mat<-cov.mat.phase1+cov.mat.phase2
        theta.se<-sqrt(diag(cov.mat))

        para.95UL<-as.numeric(para.est[1:n.regpara]+1.96*theta.se[1:n.regpara])
        para.95LL<-as.numeric(para.est[1:n.regpara]-1.96*theta.se[1:n.regpara])
      }
    }
  }

  ####### Labeling #######################################################
  if (n.gamma==0){
    para.label<-colnames(design.mat)
    para.label2<-rep("logit",n.beta)
  } else if (n.gamma>0){
    para.label<-c(colnames(design.mat),colnames(xmat))
    para.label2<-c(rep("logit",n.beta),rep("Cox",n.gamma))
  }
  ####### variance for exp(beta) and exp(gamma) ###########################
  if(exists("cov.mat")==TRUE){
    if(n.regpara>1){
      exp.cov.mat<-diag(exp(para.est[1:n.regpara]))%*% cov.mat[1:n.regpara,1:n.regpara] %*%diag(exp(para.est[1:n.regpara])) #by the delta method
      exp.theta.se<-sqrt(diag(exp.cov.mat))
    } else if(n.regpara==1){
      exp.cov.mat<-exp(para.est[1])* cov.mat[1:n.regpara,1:n.regpara] * exp(para.est[1]) #by the delta method
      exp.theta.se<-sqrt(exp.cov.mat)
    }

    exp.para.95UL<-as.numeric(exp(para.est[1:n.regpara])+1.96*exp.theta.se[1:n.regpara])
    exp.para.95LL<-as.numeric(exp(para.est[1:n.regpara])-1.96*exp.theta.se[1:n.regpara])
  }

  if(exists("theta.se")==TRUE){
    output$regression.coef<-data.frame(para.label2,para.label,para.est[1:n.regpara],theta.se[1:n.regpara],
                                       para.95LL, para.95UL, exp(para.est[1:n.regpara]),exp.theta.se[1:n.regpara],
                                       exp.para.95LL,exp.para.95UL)
    colnames(output$regression.coef)<-c('Model','Label','Coef.','SE','95%LL','95%UL','exp(Coef.)','exp.SE','exp.95%LL','exp.95%UL')


    output$exp.spline.coeff<-exp(spline.para.est)
    output$cum.hazard<-data.frame(time=time.points2, cum.hazard=Lambda.est)



    if(exists("cum.risk.est")==TRUE){
      if(n.gamma>0){
        cr.se<-CR.se(knots2, CR.time.points,cov.mat, spline.para.est, logit.risk.est,cox.risk.est,logit.predictor,cox.predictor)

        cr.95UL<-as.numeric(cum.risk.est+1.96*cr.se)
        cr.95LL<-as.numeric(cum.risk.est-1.96*cr.se)

      }else if(n.gamma==0){
        cr.se<-CR.se.noGamma(knots2, CR.time.points,cov.mat, spline.para.est, logit.risk.est,logit.predictor)

        cr.95UL<-as.numeric(cum.risk.est+1.96*cr.se)
        cr.95LL<-as.numeric(cum.risk.est-1.96*cr.se)
      }
      cr.time.labels<-paste0('time=',CR.time.points)
      output$cumrisk.est<-data.frame(cr.time.labels,cum.risk.est,cr.se,cr.95LL,cr.95UL)
      colnames(output$cumrisk.est)<-c("Time","Cumulative.Risk","SE","95%LL","95%UL")

      output$cumrisk.est.full<-data.frame(time.points2,cum.risk.est.full)
      colnames(output$cumrisk.est.full)<-c("Time","Cumulative.Risk")

    }

    output$covariance<-cov.mat
    output$hessian<- H

    output$convergence<-c(run.time=run.time,iteration=iteration,convergence=conv.marker)
    output$loglikelihood<-obs.like(para.est,samp.data,n.beta,n.gamma,design.mat,xmat,n.knots,order)

  }    else if(exists("theta.se")==FALSE){
    warning("Hessian is singular, so standard error is not provided.")


    output$regression.coef<-data.frame(para.label2,para.label,para.est[1:n.regpara],exp(para.est[1:n.regpara]), stringsAsFactors=FALSE)
    colnames(output$regression.coef)<-c('Model','Label','Coef.','exp(Coef.)')
    rownames(output$regression.coef)<-NULL

    output$exp.spline.coeff<-exp(spline.para.est)
    output$cum.hazard<-data.frame(time=time.points2, cum.hazard=Lambda.est)

    if(exists("cum.risk.est")==TRUE){

      cr.time.labels<-paste0('time=',CR.time.points)
      output$cumrisk.est<-data.frame(cr.time.labels,cum.risk.est)
      colnames(output$cumrisk.est)<-c("Time","Cumulative.Risk")

      output$cumrisk.est.full<-data.frame(time.points2,cum.risk.est.full)
      colnames(output$cumrisk.est.full)<-c("Time","Cumulative.Risk")
    }


    output$convergence<-c(run.time=run.time,iteration=iteration,convergence=conv.marker)

    output$loglikelihood<-obs.like(para.est,samp.data,n.beta,n.gamma,design.mat2=design.mat,xmat2=xmat,n.knots,order)

    output$hessian<- H
  }

  #output$logit.design<-design.mat
  #output$cox.design<-xmat
  output$knots<-knots2

  cols.1<- c('Model','Label','Coef.','SE','95%LL','95%UL')
  index.1<- match(tolower(cols.1), tolower(names(output$regression.coef)))
  index.1<- index.1[!is.na(index.1)]
  coef<- output$regression.coef[, index.1]

  cols.2<- c('Model','Label','exp(Coef.)','exp.SE','exp.95%LL','exp.95%UL')
  index.2<- match(tolower(cols.2), tolower(names(output$regression.coef)))
  index.2<- index.2[!is.na(index.2)]
  HR<- output$regression.coef[output$regression.coef[,"Model"]=="Cox",   index.2]
  OR<- output$regression.coef[output$regression.coef[,"Model"]=="logit", index.2]

  #change the colnames
  names(HR)[names(HR)=='exp(Coef.)'] <- "HR"
  names(HR)[names(HR)=='exp.SE']<- "SE"
  names(HR)[names(HR)=='exp.95%LL']<- "95%LL"
  names(HR)[names(HR)=='exp.95%UL']<- "95%UL"

  names(OR)[names(OR)=='exp(Coef.)'] <- "OR"
  names(OR)[names(OR)=='exp.SE']<- "SE"
  names(OR)[names(OR)=='exp.95%LL']<- "95%LL"
  names(OR)[names(OR)=='exp.95%UL']<- "95%UL"

  output$regression.coef<- coef
  output$HR<- HR
  output$OR<- OR

  output$model<- "weakly-parametric"
  return(output)

}#the end of the function ipw.lc.splines



ipw.lc.semipara<-function(samp.data,data,fml,fml2, n.beta, n.gamma, design.mat, xmat,
                          reg.initials,convergence.criteria,iteration.limit,
                          time.interval,logit.predictor,cox.predictor,CR.time.points, ...){

  ######################### semiparametric weighted log-like ########################################################
  wobs.semipara.like<-function(para,n.beta,n.gamma,design.mat,xmat,sdata3){
    #sdata3 should include the variables of LambdaL and LambddaR

    n.para<-n.beta+n.gamma
    mm<-n.beta+1

    beta<-para[1:n.beta]
    if(n.gamma>0){
      gamma<-para[mm:n.para]
      gamma.risk<-xmat%*%gamma
    } else if(n.gamma==0){
      gamma<-0
      gamma.risk<-0
    }

    beta.risk<-design.mat%*%beta


    B<-log(    exp(-sdata3$LambdaL*exp(gamma.risk))-exp(-sdata3$LambdaR*exp(gamma.risk)) )
    B[sdata3$C==1]<-0                                            #When C==1, L=R=999
    B[abs(sdata3$LambdaL-sdata3$LambdaR)<1e-10]<-0               #to avoid the value of Inf/-Inf when 1/B and B goes to 0


    B2<-1-exp(-sdata3$LambdaR*exp(gamma.risk))                   #K=0
    B2[sdata3$K==1]<-0


    wll_K1<-as.numeric(sdata3$K==1)*sdata3$samp.weight*(sdata3$C*beta.risk-log(1+exp(beta.risk)) + B)
    wll_K2<-as.numeric(sdata3$K==0)*sdata3$samp.weight*log(exp(beta.risk)/(1+exp(beta.risk))+B2/(1+exp(beta.risk)) )

    #pseudo-log-likelihood
    return(sum(wll_K1+wll_K2))
  }

  wsemi.opt <- function(para){
    -1*wobs.semipara.like(para,n.beta,n.gamma,design.mat,xmat,sdata3)
  }


  ###############################################################################################
  #time.points2<-time points to be set for cumulative hazard estimate
  #fixed.initials initial parameters for beta and gamma just in case when logistic model doesn't converge.

  n.samp<-nrow(samp.data)
  n.regpara<-n.beta+n.gamma
  mm<-n.beta+1

  iter.limit<-400                          #iteration for cumulative hazard function

  ########### Initial parameters  ################################################

  if(is.null(reg.initials)){

    if(n.gamma==0){
      initials<-rep(0.05,n.beta)
    } else if (n.gamma>0){
      initials<-c(rep(0.05,n.beta),rep(0,n.gamma))
    }
  } else {
    if(length(reg.initials)!=n.regpara){
      stop("The number of regression parameter initials are not matched with the model.")
    } else {
      initials<-reg.initials
    }
  }

  missing.upper<-max(samp.data$R[samp.data$R<Inf],samp.data$L)*1.5+2 #the initial Lambda=time.list*1.5 and 2 is an arbitrary positive number

  samp.data$L[samp.data$L==-999]<-missing.upper
  samp.data$R[samp.data$R==-999]<-missing.upper

  ############## Lambda Initials #################################################

  Jn<-sort(unique(c(samp.data$L[samp.data$L>0&samp.data$L<missing.upper],samp.data$R[samp.data$R<missing.upper])))

  #ordered statistics for time: min{X_i: T_i<=X_i} and max{X_i:X_i<=T_i}
  time.1<-min(samp.data$R[samp.data$R<missing.upper])
  time.2<-max(samp.data$L[samp.data$L>0&samp.data$L<missing.upper])

  if(time.2==0){
    stop("NPMLE condition: the maximum left time point should be positive.")
  }

  time.points<-seq(0,time.2,by=time.interval)

  time.list<-Jn[Jn>=time.1&Jn<=time.2]
  time.95<-quantile(time.list,c(0.95))

  Lambda<-time.list*1.5
  old.Lambda<-Lambda.data<-data.frame(time=time.list,Lambda=Lambda)

  Lambda.index<-rbind(Lambda.data,rep(0,2),rep(Inf,2),rep(missing.upper,2))
  sdata<-samp.data

  sdata2<-merge(sdata,Lambda.index,by.x="L",by.y="time",all.x=TRUE)
  colnames(sdata2)[ncol(sdata2)]<-"LambdaL"
  sdata3<-merge(sdata2,Lambda.index,by.x="R",by.y="time",all.x=TRUE)
  colnames(sdata3)[ncol(sdata3)]<-"LambdaR"

  sdata3$LambdaL[sdata3$L<time.1]<-0
  sdata3$LambdaL[sdata3$L>time.2&sdata3$L<missing.upper]<-max(sdata3$LambdaL[sdata3$LambdaL<missing.upper],na.rm=TRUE)
  sdata3$LambdaR[sdata3$R>time.2&sdata3$R<missing.upper]<-max(sdata3$LambdaL[sdata3$LambdaL<missing.upper],na.rm=TRUE)

  remove(sdata,sdata2)
  #remove(design.mat,xmat)

  ########################### LOOPING ###########################################################
  old.wNR.est<-initials
  iteration<-0
  keep.going<-1
  ptm<-proc.time()
  diff.estimates<-100

  while(keep.going==1){

    iteration<- iteration +1

    ######################   Newton-Raphson               ##############################
    #print(old.wNR.est)
    #print(n.beta)
    #print(n.gamma)
    #print(ncol(design.mat))
    #print(ncol(xmat))
    #print(ncol(sdata3))
    if(n.gamma==0){
      resw <- optim(old.wNR.est,wsemi.opt,method="L-BFGS-B")
    }else if(n.gamma>0){
      resw <- optim(old.wNR.est,wsemi.opt)
    }


    wNR.est<-resw$par

    diff.reg.para<-max(abs(old.wNR.est-wNR.est))

    ######################   Cumulative Hazard Estimate   ##############################

    mf <- model.frame(formula=fml, data=samp.data)
    design.mat <- model.matrix(attr(mf, "terms"), data=mf)

    if(is.null(fml2)){
      if(n.beta>1){
        xmat<-as.matrix(design.mat[,-1])
        colnames(xmat)<-colnames(design.mat)[-1]
      }else if(n.beta==1){
        xmat<-0
      }
    } else if(!is.null(fml2)){
      mf2 <- model.frame(formula=fml2, data=samp.data)
      xmat <- model.matrix(terms(fml2), data=mf2)

      if(ncol(xmat)>1){
        xmat<-as.matrix(xmat[,-1])
      } else if(ncol(xmat)==1){
        xmat<-0
      }
    }

    expplus.est<-as.numeric(1+exp(design.mat%*%wNR.est[1:n.beta]))
    if(n.gamma>0){
      cox.hr.est<-as.numeric(exp(xmat%*%wNR.est[mm:n.regpara]))
    }else {
      cox.hr.est<-1
    }


    iter<-0
    keep.looping<-1

    while(keep.looping==1){

      iter<-iter+1

      Lambda.index<-rbind(Lambda.data,rep(0,2),rep(Inf,2),rep(missing.upper,2))


      sdata<-samp.data

      sdata$WL1<-rep(0,n.samp)
      sdata$WL2<-rep(0,n.samp)
      sdata$WL3<-rep(0,n.samp)
      sdata$HR<-cox.hr.est
      sdata$expplus<-expplus.est


      sdata2<-merge(sdata,Lambda.index,by.x="L",by.y="time",all.x=TRUE)
      colnames(sdata2)[ncol(sdata2)]<-"LambdaL"
      sdata3<-merge(sdata2,Lambda.index,by.x="R",by.y="time",all.x=TRUE)
      colnames(sdata3)[ncol(sdata3)]<-"LambdaR"
      sdata3$LambdaL[sdata3$L<time.1&sdata3$L>=0]<-0
      sdata3$LambdaL[sdata3$L>time.2&sdata3$L<missing.upper]<-max(sdata3$LambdaL[sdata3$LambdaL<missing.upper],na.rm=TRUE)
      sdata3$LambdaR[sdata3$R>time.2&sdata3$R<missing.upper]<-max(sdata3$LambdaL[sdata3$LambdaL<missing.upper],na.rm=TRUE)

      remove(sdata,sdata2)

      exc1<-which(abs(sdata3$LambdaR-sdata3$LambdaL)<1e-10)
      exc2<-which(is.na(sdata3$LambdaL))
      exc3<-which(is.na(sdata3$LambdaR))


      subset1<-which(sdata3$K==0&sdata3$R<missing.upper&sdata3$L==0)                 #K0 and C=-999, 0<=T<R<+Infinity :derivative at right point of interval
      subset2<-which(sdata3$K==1&sdata3$C==0&sdata3$L>0)                              #K1 and C=0, 0<L<T<= +Infinity :derivative at left point of interval
      subset3<-which(sdata3$K==1&sdata3$C==0&sdata3$R<missing.upper)                            #K1 and C=0, R<+Infinity :derivative at right point of interval



      # weighted process
      WL1<-sdata3$samp.weight*exp(-sdata3$LambdaR*sdata3$HR)*sdata3$HR/(sdata3$expplus-exp(-sdata3$LambdaR*sdata3$HR))
      sdata3$WL1[subset1]<-WL1[subset1]

      # weighted process
      WL2<- -exp(-sdata3$LambdaL*sdata3$HR)*sdata3$HR/(exp(-sdata3$LambdaL*sdata3$HR)-exp(-sdata3$LambdaR*sdata3$HR))*sdata3$samp.weight
      sdata3$WL2[subset2]<-WL2[subset2]
      sdata3$WL2[union(exc1,exc2)]<-0

      # weighted process
      WL3<- sdata3$samp.weight*exp(-sdata3$LambdaR*sdata3$HR)*sdata3$HR/(exp(-sdata3$LambdaL*sdata3$HR)-exp(-sdata3$LambdaR*sdata3$HR))
      sdata3$WL3[subset3]<-WL3[subset3]
      sdata3$WL3[union(exc1,exc3)]<-0

      sdata3$WL1.sq<-sdata3$WL1^2
      sdata3$WL2.sq<-sdata3$WL2^2
      sdata3$WL3.sq<-sdata3$WL3^2


      k<-0
      Q_L<-0
      W_L<-0
      G_L<-0
      TL<-nrow(Lambda.data)
      Q_Lambda<-dG_Lambda<-G_Lambda<-W_Lambda<-rep(NA,TL)

      for (ordered.t in Lambda.data$time){
        k<-k+1
        L.process<-as.numeric(sdata3$L<=ordered.t)
        R.process<-as.numeric(sdata3$R<=ordered.t)

        L2.process<-as.numeric(sdata3$L==ordered.t)
        R2.process<-as.numeric(sdata3$R==ordered.t)


        W_Lambda[k]<-sum(R.process*(sdata3$WL1+sdata3$WL3)+L.process*sdata3$WL2,na.rm=TRUE)
        G_Lambda[k]<-sum(R.process*(sdata3$WL1.sq+sdata3$WL3.sq)+L.process*sdata3$WL2.sq,na.rm=TRUE)

        if(k==1){dG_Lambda[k]<-G_Lambda[k]
        } else if(k>1){
          dG_Lambda[k]<-G_Lambda[k]-G_Lambda[k-1]
        }
        Q_L<- Q_L+Lambda.data$Lambda[k]*dG_Lambda[k]
        Q_Lambda[k]<-W_Lambda[k]+Q_L
      }


      no.change<-which(dG_Lambda==0)
      if(length(no.change)>0){
        G_Lambda<-G_Lambda[-no.change]
        Q_Lambda<-Q_Lambda[-no.change]
        Lambda.data<-Lambda.data[-no.change,]
      }


      G_Lambda<-c(0,G_Lambda)
      Q_Lambda<-c(0,Q_Lambda)


      GCM = gcmlcm(G_Lambda, Q_Lambda)

      time.knots<-which(G_Lambda[-1] %in% GCM$x.knots[-1])

      update.Lambda<-data.frame(time=Lambda.data$time[time.knots],Lambda=GCM$slope.knots)
      update.Lambda$Lambda[update.Lambda$Lambda<0]<-0

      sfun.1  <- stepfun(update.Lambda$time, c(update.Lambda$Lambda[1],update.Lambda$Lambda), f = 1)
      new.Lambda<-data.frame(time=time.list,Lambda=sfun.1(time.list))
      diff.Lambda<-max(abs(new.Lambda$Lambda[time.list<time.95]-old.Lambda$Lambda[time.list<time.95]))

      Lambda.data<-old.Lambda<-new.Lambda



      if (diff.Lambda<convergence.criteria|iter>iter.limit){
        keep.looping<-0
      }
    } #while (keep.looping==1)



    diff.estimates<-max(diff.reg.para,diff.Lambda)
    diff.reg.est<-diff.reg.para
    diff.Lambda.est<-diff.Lambda

    old.wNR.est<-wNR.est


    ######################To incorporate the updated Lambda in the regression estimation ##################################
    ########################## design matrix within the loop for cumulative hazard estimate ###############################
    mf <- model.frame(formula=fml, data=sdata3)
    design.mat <- model.matrix(attr(mf, "terms"), data=mf)

    if(is.null(fml2)){
      if(n.beta>1){
        xmat<-as.matrix(design.mat[,-1])
        colnames(xmat)<-colnames(design.mat)[-1]
      }else if(n.beta==1){
        xmat<-0
      }
    } else if(!is.null(fml2)){
      mf2 <- model.frame(formula=fml2, data=sdata3)
      xmat <- model.matrix(terms(fml2), data=mf2)
      if(ncol(xmat2)>1){
        xmat<-as.matrix(xmat[,-1])
      } else if(ncol(xmat)==1){
        xmat<-0
      }
    }

    if(diff.estimates<convergence.criteria|iteration>iteration.limit){
      par.est<-wNR.est
      length.nonpara<-length(GCM$x.knots)
      final.Lambda<-data.frame(Lambda=sfun.1(time.points))
      Lambda.est<-sfun.1(time.points)
      keep.going<-0
    }



  } #while(keep.going==1){

  run.time<-((proc.time()-ptm)[3])/60


  ## Prevalence/Cumulative Risk Estimation #########################################################
  if(!missing(logit.predictor)&!missing(CR.time.points)){
    logit.risk.est<-as.numeric(logit.predictor%*%par.est[1:n.beta])
    prevalence.est<-as.numeric(exp(logit.risk.est)/(1+exp(logit.risk.est)))
  }
  if(!missing(cox.predictor)&!missing(CR.time.points)){
    cox.risk.est<- as.numeric(cox.predictor%*%par.est[mm:n.regpara])
  } else if(n.gamma==0&!missing(CR.time.points)&!missing(logit.predictor)){
    cox.risk.est<- 0
  }

  if(exists("prevalence.est")==TRUE&exists("cox.risk.est")==TRUE){
    cum.risk.est<-cumulative.risk2(prevalence.est, Lambda.est[which(time.points %in% CR.time.points)] ,cox.risk.est)
    cum.risk.est.full<-cumulative.risk2(prevalence.est, Lambda.est ,cox.risk.est)
  }

  ####### Labeling For Regression Parameter Estimates#######################################################
  ## newly labeling #04/22/2019

  length.mf<-length(colnames(mf))
  if (length.mf>1){
    for(j in 2:length.mf){
      if(is.factor(mf[,j])){
        level.list<-levels(mf[,j])

        for (k in 1:length(level.list)) {
          colnames(design.mat)<-gsub(level.list[k],paste0("=",level.list[k]),colnames(design.mat),fixed=T)
        }

        remove(level.list)
      }
    } #for(j in 2:length.mf)
  }#if (length.mf>1)
  remove(length.mf)

  if(!exists("mf2")){
    if(n.beta>1){
      xmat<-matrix(design.mat[,-1],  nrow=nrow(design.mat))
      colnames(xmat)<-colnames(design.mat)[-1]
    }else if(n.beta==1){
      xmat<-0
    }
  } else if(exists("mf2")){
    length.mf2<-length(colnames(mf2))
    if (length.mf2>1){
      for(j in 2:length.mf2){
        if(is.factor(mf2[,j])){
          level.list<-levels(mf2[,j])

          for (k in 1:length(level.list)) {
            colnames(xmat)<-gsub(level.list[k],paste0("=",level.list[k]),colnames(xmat),fixed=T)
          }

          remove(level.list)
        }
      } #for(j in 2:length.mf)
    }#if(length.mf2>1)
    remove(length.mf2)
  }#else if(exists("mf2"))




  if (n.gamma==0){
    para.label<-colnames(design.mat)
    para.label2<-rep("logit",n.beta)
  } else if (n.gamma>0){
    para.label<-c(colnames(design.mat),colnames(xmat))
    para.label2<-c(rep("logit",n.beta),rep("Cox",n.gamma))
  }

  output<- list()
  output$regression.coef<-data.frame(para.label2,para.label,par.est,exp(par.est), stringsAsFactors=FALSE)
  colnames(output$regression.coef)<-c('Model','Label','Coef.','exp(Coef.)')

  output$cum.hazard<-data.frame(time=time.points,cum.hazard=final.Lambda$Lambda)



  if(exists("cum.risk.est")==TRUE){
    cr.time.labels<-paste0('time=',CR.time.points)
    output$cumrisk.est<-data.frame(cr.time.labels,cum.risk.est)
    colnames(output$cumrisk.est)<-c("Time","Cumulative.Risk")

    output$cumrisk.est.full<-data.frame(time.points,cum.risk.est.full)
    colnames(output$cumrisk.est.full)<-c("Time","Cumulative.Risk")
  }


  output$convergence<-c(run.time=run.time,iteration=iteration,regpara.conv=diff.reg.est,Lambda.conv=diff.Lambda.est,no.of.nonpara=length.nonpara)
  output$log.pseudolike<-wobs.semipara.like(par.est,n.beta,n.gamma,design.mat,xmat,sdata3)

  cols.1<- c('Model','Label','Coef.','SE','95%LL','95%UL')
  index.1<- match(tolower(cols.1), tolower(names(output$regression.coef)))
  index.1<- index.1[!is.na(index.1)]
  coef<- output$regression.coef[, index.1]

  cols.2<- c('Model','Label','exp(Coef.)','SE','95%LL','95%UL')
  index.2<- match(tolower(cols.2), tolower(names(output$regression.coef)))
  index.2<- index.2[!is.na(index.2)]
  HR<- output$regression.coef[output$regression.coef[,1]=="Cox",   index.2]
  OR<- output$regression.coef[output$regression.coef[,1]=="logit", index.2]

  names(HR)[names(HR)=='exp(Coef.)'] <- "HR"
  names(OR)[names(OR)=='exp(Coef.)'] <- "OR"
  output$regression.coef<- coef
  output$HR<- HR
  output$OR<- OR

  output$model<- "semi-parametric"

  return(output)
}  # the end of the function ipw.lc.semipara



check.data<- function(samp.data, allvars) {
  n.samp<- nrow(samp.data)

  if(nrow(samp.data)==0){
    stop("All observations have NA values, and subjects with NA values are excluded from analysis")
  }

  if (sum(c("C","L","R") %in% (names(samp.data)))!=3){
    stop("Response variables should include three variables and has the form of \"prevalence indicator + left time point + right time point\"")
  }

  if(sum(samp.data$C %in% c(0,1,-999))!=n.samp){
    stop(paste0("Prevalence indicator,",allvars[1]," should have values of 0 (non-prevalence), 1 (prevalence), -999 only. "))
  }

  if(sum(samp.data$L<0&samp.data$L!=-999)>0){
    stop(paste0("Left time points,",allvars[2]," should have non-negative values except when it is coded as -999 for prevalent cases"))
  }

  if(sum(samp.data$R<0&samp.data$R!=-999)>0){
    stop(paste0("Right time points,",allvars[3]," should have non-negative values except when it is coded as -999 for prevalent cases"))
  }

  if (sum(samp.data$C==0)==n.samp){
    stop(paste0("Prevalence indicator,",allvars[1]," has all 0 (non-prevalence) value, so prevalence-incidence model is not applicable."))
  }

}

get.summary<- function(samp.data) {
  n.samp<-nrow(samp.data)
  L.class<-ifelse(samp.data$L==0,0,
                  ifelse(samp.data$L==-999,-999,1))
  R.class<-ifelse(samp.data$R==-999,-999,
                  ifelse(samp.data$R==Inf,Inf,
                         ifelse(samp.data$R==0,0,1)))
  #no.info.case<-sum(as.numeric(L.class==0&R.class==Inf))
  interval.censoring.case<-sum(as.numeric(L.class==0&R.class==1&samp.data$C==0))+ sum(as.numeric(L.class==1&R.class==1))
  left.censoring.case<-sum(as.numeric(L.class==0&R.class==1&samp.data$C==-999))
  right.censoring.case<-sum(as.numeric(L.class==1&R.class==Inf))
  prevalent.case<-sum(as.numeric(L.class==-999&R.class==-999))
  missing.prev.status<-sum(samp.data$C==-999)
  non.info.interval<-  sum(samp.data$C != 1 & samp.data$L ==0 & samp.data$R == Inf)

  data.summary<- data.frame(c("Included subjects","Known prevalent cases","Interval censoring","Left censoring","Right censoring","Missing prevalent status", "Non-informative interval"),
                            c(n.samp,prevalent.case,interval.censoring.case,left.censoring.case,right.censoring.case,missing.prev.status, non.info.interval), stringsAsFactors=FALSE)

  names(data.summary)<- c("label", "no.obs")

  return(data.summary)
}

###########################################################
## PIMixture estimate for parametric model
## init: the initial values for parameters in the model
##
###########################################################

PIMixtureEst <- function(ci,li,ri,mat1=NULL,mat2=NULL,model="logistic-Weibull",conf.int=FALSE, init=NULL)
{
  requireNamespace("survival", quietly = TRUE)
  requireNamespace("optimx", quietly = TRUE)
  requireNamespace("fdrtool", quietly = TRUE)

  n <- length(ci)
  if (length(li) != n | length(ri) != n) { cat("\nSTOP\nInput vectors are of different lengths\n") }

  if(!is.null(mat1)) {
    mat1<- as.matrix(mat1)
  }

  if(!is.null(mat2)) {
    mat2<- as.matrix(mat2)
  }

  ####### PRELIMINARIES #######
  n1 <- sum(!is.na(ci) & ci==1)  #number of strictly prevalent disease
  if (is.na(n1)) n1 <- 0
  n2 <- sum(!is.na(ci) & ci==0 & ri<Inf)  #number of strictly incident disease
  if (is.na(n2)) n2 <- 0
  n3 <- sum(is.na(ci) & ri<Inf)  #number of disease that can be prevalent or incident
  if (is.na(n3)) n3 <- 0

  xsam <- matrix(NA,n,5)
  xsam[,1] <- !is.na(ci)
  xsam[,2] <- ci
  xsam[,3] <- li
  xsam[,4] <- ri
  xsam[,5] <- ifelse(li>0 | ri<Inf,1,0)
  colnames(xsam) <- c("k","ci","li","ri","k2")
  usable <- sum(xsam[,1]==1 | xsam[,5]==1)
  names(usable) <- "usable"

  m1 <- 0
  m2 <- 0

  if (!is.null(mat1)) {
    names.mat1<- colnames(mat1)
    m1 <- ncol(as.matrix(mat1))
    mat1 <- matrix(mat1[xsam[,1]==1 | xsam[,5]==1,], ncol=ncol(mat1))
    colnames(mat1)<- names.mat1
  }
  if (!is.null(mat2)) {
    names.mat2<- colnames(mat2)
    m2 <- ncol(as.matrix(mat2))
    mat2 <- matrix(mat2[xsam[,1]==1 | xsam[,5]==1,], ncol=ncol(mat2))
    colnames(mat2)<- names.mat2
  }

  xsam <- matrix(xsam[xsam[,1]==1 | xsam[,5]==1,], ncol=ncol(xsam))

  ###############################
  ####### START FUNCTIONS #######
  ###############################

  ####### Logistic-Weibull Functions #######
  #find initial parameter estimates
  logistic_weibull_init <- function()
  {
    # make crude assumptions to get initial parameter estimates
    xsam2 <- xsam
    if (n1>n2)
    {
      xsam2[which(is.na(xsam2[,2])),2] <- 1
    } else if (n1<=n2)
    {
      xsam2[which(is.na(xsam2[,2])),2] <- 0
    }
    xsam2[which(xsam2[,3]==0),3] <- NA
    xsam2[which(xsam2[,4]==Inf),4] <- NA

    if (m1==0)
    {
      logistic_res <- glm(xsam2[,2] ~ 1, family = "binomial")
      param_init <- logistic_res$coefficients
      names(param_init) <- "logit-(Intercept)"
    } else {
      logistic_res <- glm(xsam2[,2] ~ mat1, family = "binomial")
      param_init <- logistic_res$coefficients
      if (is.null(colnames(mat1))) {
        names(param_init)<- paste0("logit-", names(logistic_res$coefficients))
      } else{
        names(param_init)<- paste0("logit-", c("(Intercept)", colnames(mat1)))
      }
    }

    subset<- which(xsam2[,2]==0)
    if (m2==0){
      llog <-     survreg(Surv(xsam2[,3],xsam2[,4],type='interval2')~1, dist="loglogistic", subset=subset)
      weib_res <- survreg(Surv(xsam2[,3],xsam2[,4],type='interval2')~1, dist="weibull",
                          init=c(llog$coefficients,llog$scale), subset=subset)

      param_init2 <- c(weib_res$coefficients, weib_res$scale)
      names(param_init2) <- paste0("Weibull-", c("log(scale)", "1/shape"))

    } else {
      llog <-     survreg(Surv(xsam2[,3],xsam2[,4],type='interval2')~ mat2, dist="loglogistic", subset=subset)
      weib_res <- survreg(Surv(xsam2[,3],xsam2[,4],type='interval2')~ mat2, dist="weibull",
                          init=c(llog$coefficients,llog$scale), subset=subset)
      param_init2 <- c(weib_res$coefficients, weib_res$scale)
      if(is.null(colnames(mat2))) {
        names(param_init2)<- paste0("Weibull-", c(names(weib_res$coefficients), "1/shape"))
      } else {
        names(param_init2)<- paste0("Weibull-", c("(Intercept)", colnames(mat2),"1/shape"))
      }
    }

    param_init <- c(param_init,param_init2)

    return(param_init)
  }


  #full log-likelihood
  logistic_weibull_likfun <- function(parest)
  {
    k <- xsam[,1]
    ci <- xsam[,2]
    li <- xsam[,3]
    ri <- xsam[,4]
    k2 <- xsam[,5]

    if (m1==0) {p_lp <- rep(parest[1],usable)
    } else if (m1==1) { p_lp<-parest[1]+mat1*parest[2]
    } else if (m1>1) {p_lp <- parest[1]+mat1%*%parest[2:(m1+1)] }

    if (m2==0)
    {
      s_lp <- rep(parest[(m1+2)],usable)
      SL <- 1-pweibull(li,1/parest[(m1+3)],exp(s_lp))
      SR <- 1-pweibull(ri,1/parest[(m1+3)],exp(s_lp))
    } else if (m2==1)
    {
      s_lp <- parest[(m1+2)]+mat2*parest[(m1+3)]
      SL <- 1-pweibull(li,1/parest[(m1+m2+3)],exp(s_lp))
      SR <- 1-pweibull(ri,1/parest[(m1+m2+3)],exp(s_lp))
    } else if (m2>1)
    {
      s_lp <- parest[(m1+2)]+mat2%*%parest[(m1+3):(m1+m2+2)]
      SL <- 1-pweibull(li,1/parest[(m1+m2+3)],exp(s_lp))
      SR <- 1-pweibull(ri,1/parest[(m1+m2+3)],exp(s_lp))
    }

    #log-likelihood
    tmp <- rep(NA,usable)
    tmp[which(ci==1 & k==1)] <- I(p_lp-log(1+exp(p_lp)))[which(ci==1 & k==1)]
    tmp[which(ci==0 & k==1)] <- I(log(SL-SR)-log(1+exp(p_lp)))[which(ci==0 & k==1)]
    tmp[which(k==0 & k2==1)] <- I(log(exp(p_lp)+1-SR)-log(1+exp(p_lp)))[which(k==0 & k2==1)]
    tmp[which(k==0 & k2==0)] <- 0
    sum(tmp) #full log-likelihood
  }

  #minimize function to maximize likelihood
  logistic_weibull_opt <- function(y)
  {
    -1*logistic_weibull_likfun(y)
  }

  logistic_weibull_gradem <- function(parinit){
    #code EM algorithm
    estep <- function(parest){
      k <- xsam[,1]
      ci <- xsam[,2]
      ri <- xsam[,4]

      if (m1==0) {p_lp <- rep(parest[1],usable)
      } else if (m1==1) {p_lp <- parest[1]+mat1*parest[2]
      } else if (m1>1) {p_lp <- parest[1]+mat1%*%parest[2:(m1+1)] }

      if (m2==0)
      {
        s_lp <- rep(parest[(m1+2)],usable)
        SR <- 1-pweibull(ri,1/parest[(m1+3)],exp(s_lp))
      } else if (m2==1)
      {
        s_lp <- parest[(m1+2)]+mat2*parest[(m1+3)]
        SR <- 1-pweibull(ri,1/parest[(m1+m2+3)],exp(s_lp))
      } else if (m2>1)
      {
        s_lp <- parest[(m1+2)]+mat2%*%parest[(m1+3):(m1+m2+2)]
        SR <- 1-pweibull(ri,1/parest[(m1+m2+3)],exp(s_lp))
      }

      ci_prime <- rep(NA,usable)
      ci_prime[k==1] <- ci[k==1]
      ci_prime[k==0] <- I(exp(p_lp)/(exp(p_lp)+1-SR))[k==0]
      return(ci_prime)
    }

    mstep <- function(parest,ci_prime){
      ci <- ci_prime
      li <- xsam[,3]
      ri <- xsam[,4]

      if (m1==0)
      {
        p_lp <- rep(parest[1],usable)
        # create design matrices
        mat1 <- matrix(1,usable,1)
      } else if (m1==1)
      {
        p_lp <- parest[1]+mat1*parest[2]
        # create design matrices
        mat1 <- cbind(rep(1,usable),mat1)
      } else if (m1>1)
      {
        p_lp <- parest[1]+mat1%*%parest[2:(m1+1)]
        # create design matrices
        mat1 <- cbind(rep(1,usable),mat1)
      }

      if (m2==0)
      {
        s_lp <- rep(parest[(m1+2)],usable)
        SL <- 1-pweibull(li,1/parest[(m1+3)],exp(s_lp))
        SR <- 1-pweibull(ri,1/parest[(m1+3)],exp(s_lp))
        mat2 <- matrix(1,usable,1)
      } else if (m2==1)
      {
        s_lp <- parest[(m1+2)]+mat2*parest[(m1+3)]
        SL <- 1-pweibull(li,1/parest[(m1+m2+3)],exp(s_lp))
        SR <- 1-pweibull(ri,1/parest[(m1+m2+3)],exp(s_lp))
        # create design matrices
        mat2 <- cbind(rep(1,usable),mat2)
      } else if (m2>1)
      {
        s_lp <- parest[(m1+2)]+mat2%*%parest[(m1+3):(m1+m2+2)]
        SL <- 1-pweibull(li,1/parest[(m1+m2+3)],exp(s_lp))
        SR <- 1-pweibull(ri,1/parest[(m1+m2+3)],exp(s_lp))
        # create design matrices
        mat2 <- cbind(rep(1,usable),mat2)
      }

      p <- exp(p_lp)/(1+exp(p_lp))
      tau <- parest[m1+m2+3]

      #score matrix
      u <- rep(NA,m1+m2+3)

      for (i in 1:(m1+1))
      {
        u[i] <- sum((ci-p)*mat1[,i])
      }

      for (i in 1:(m2+1))
      {
        u_gamma_a1 <- sum(((1-ci)*{li**{1/tau}*SL-ri**{1/tau}*SR} / {tau*exp(s_lp/tau)*(SL-SR)}*mat2[,i])[ri<Inf & ci!=1])
        u_gamma_a2 <- sum(((1-ci)*{li**{1/tau}*SL} / {tau*exp(s_lp/tau)*SL}*mat2[,i])[ri==Inf])
        u[m1+1+i] <- u_gamma_a1 + u_gamma_a2
      }

      u_tau_a1 <- sum(((1-ci)*{(log(li)-s_lp)*li**{1/tau}*SL-(log(ri)-s_lp)*ri**{1/tau}*SR} / {tau**2*exp(s_lp/tau)*(SL-SR)})[li>0&ri<Inf])
      u_tau_a2 <- sum(((1-ci)*{-(log(ri)-s_lp)*ri**{1/tau}*SR} / {tau**2*exp(s_lp/tau)*(SL-SR)})[li==0&ri<Inf&ci!=1])
      u_tau_a3 <- sum(((1-ci)*{(log(li)-s_lp)*li**{1/tau}*SL} / {tau**2*exp(s_lp/tau)*(SL-SR)})[li>0&ri==Inf])
      u_tau_a4 <- 0 #li==0 and ri==Inf
      u[m1+m2+3] <- u_tau_a1 + u_tau_a2 + u_tau_a3 + u_tau_a4


      #hessian matrix
      h <- matrix(NA,m1+m2+3,m1+m2+3)
      for (i in 1:(m1+1))
      {
        for (j in i:(m1+1))
        {
          h[i,j] <- -1*sum(p/{1+exp(p_lp)}*mat1[,i]*mat1[,j])
          h[j,i] <- -1*sum(p/{1+exp(p_lp)}*mat1[,i]*mat1[,j])
        }
        for (j in 1:(m2+1))
        {
          h[i,m1+1+j] <- 0
          h[m1+1+j,i] <- 0
        }
        h[i,m1+m2+3] <- 0
        h[m1+m2+3,i] <- 0
      }

      for (i in 1:(m2+1))
      {
        for (j in i:(m2+1))
        {
          h_gamma_a1 <- sum(((1-ci)*{((li*exp(s_lp))**{1/tau}+(ri*exp(s_lp))**{1/tau}-(li**{1/tau}-ri**{1/tau})**2)*SL*SR-(li*exp(s_lp))**{1/tau}*SL**2-(ri*exp(s_lp))**{1/tau}*SR**2} / {(tau*exp(s_lp/tau)*(SL-SR))**2}*mat2[,i]*mat2[,j])[ri<Inf & ci!=1])
          h_gamma_a2 <- sum(((1-ci)*{-(li*exp(s_lp))**{1/tau}*SL**2} / {(tau*exp(s_lp/tau)*(SL-SR))**2}*mat2[,i]*mat2[,j])[ri==Inf])
          h[m1+1+i,m1+1+j] <- h_gamma_a1 + h_gamma_a2
          h[m1+1+j,m1+1+i] <- h_gamma_a1 + h_gamma_a2
        }

        h_gamma_tau_a1 <- sum(((1-ci)*({((li*exp(s_lp))**{1/tau}+(li*ri)**{1/tau}-li**{2/tau})*(log(li)-s_lp)/tau + (li*exp(s_lp))**{1/tau} + ((ri*exp(s_lp))**{1/tau} + (li*ri)**{1/tau}-ri**{2/tau})*(log(ri)-s_lp)/tau+(ri*exp(s_lp))**{1/tau}}*SL*SR-(1+(log(li)-s_lp)/tau)*(li*exp(s_lp))**{1/tau}*SL**2-(1+(log(ri)-s_lp)/tau)*(ri*exp(s_lp))**{1/tau}*SR**2) / ({tau*exp(s_lp/tau)*(SL-SR)}**2)*mat2[,i])[li>0&ri<Inf])
        h_gamma_tau_a2 <- sum(((1-ci)*({((ri*exp(s_lp))**{1/tau}-ri**{2/tau})*(log(ri)-s_lp)/tau+(ri*exp(s_lp))**{1/tau}}*SL*SR-(1+(log(ri)-s_lp)/tau)*(ri*exp(s_lp))**{1/tau}*SR**2) / ({tau*exp(s_lp/tau)*(SL-SR)}**2)*mat2[,i])[li==0&ri<Inf & ci!=1])
        h_gamma_tau_a3 <- sum(((1-ci)*(-(1+(log(li)-s_lp)/tau)*(li*exp(s_lp))**{1/tau}*SL**2) / ({tau*exp(s_lp/tau)*(SL-SR)}**2)*mat2[,i])[li>0&ri==Inf])
        h_gamma_tau_a4 <- 0 #0 when li==0 & ri == Inf
        h[m1+1+i,m1+m2+3] <- h_gamma_tau_a1 + h_gamma_tau_a2 + h_gamma_tau_a3 + h_gamma_tau_a4
        h[m1+m2+3,m1+1+i] <- h_gamma_tau_a1 + h_gamma_tau_a2 + h_gamma_tau_a3 + h_gamma_tau_a4
      }

      h_tau_a1 <- sum(((1-ci)*({(log(li)-s_lp+2*tau)*(log(li)-s_lp)*(li*exp(s_lp))**{1/tau} + (log(ri)-s_lp+2*tau)*(log(ri)-s_lp)*(ri*exp(s_lp))**{1/tau}-((log(li)-s_lp)*li**{1/tau}-(log(ri)-s_lp)*ri**{1/tau})**2}*SL*SR - (log(ri)-s_lp+2*tau)*(log(ri)-s_lp)*(ri*exp(s_lp))**{1/tau}*SR**2 - (log(li)-s_lp+2*tau)*(log(li)-s_lp)*(li*exp(s_lp))**{1/tau}*SL**2) / ({tau**2*exp(s_lp/tau)*(SL-SR)}**2) )[li>0&ri<Inf])
      h_tau_a2 <- sum(((1-ci)*({(log(ri)-s_lp+2*tau)*(log(ri)-s_lp)*(ri*exp(s_lp))**{1/tau} - ((log(ri)-s_lp)*ri**{1/tau})**2}*SR - (log(ri)-s_lp+2*tau)*(log(ri)-s_lp)*(ri*exp(s_lp))**{1/tau}*SR**2) / ({tau**2*exp(s_lp/tau)*(SL-SR)}**2) )[li==0&ri<Inf&ci!=1])
      h_tau_a3 <- sum(((1-ci)*(-(log(li)-s_lp+2*tau)*(log(li)-s_lp)*(li*exp(s_lp))**{1/tau}*SL**2) / ({tau**2*exp(s_lp/tau)*(SL-SR)}**2) )[li>0&ri==Inf])
      h_tau_a4 <- 0 #0 when li==0 & ri == Inf
      h[m1+m2+3,m1+m2+3] <- h_tau_a1 + h_tau_a2 + h_tau_a3 + h_tau_a4

      hessian <- -1*h

      fisher_info <- solve(hessian)

      parest_next <- parest + fisher_info%*%u

      return(parest_next)
    }

    covvar <- function(parest){
      k <- xsam[,1]
      ci <- xsam[,2]
      li <- xsam[,3]
      ri <- xsam[,4]
      k2 <- xsam[,5]

      if (m1==0)
      {
        p_lp <- rep(parest[1],usable)
        # create design matrices
        mat1 <- matrix(1,usable,1)
      } else if (m1==1)
      {
        p_lp <- parest[1]+mat1*parest[2]
        # create design matrices
        mat1 <- cbind(rep(1,usable),mat1)
      } else if (m1>1)
      {
        p_lp <- parest[1]+mat1%*%parest[2:(m1+1)]
        # create design matrices
        mat1 <- cbind(rep(1,usable),mat1)
      }

      if (m2==0)
      {
        s_lp <- rep(parest[(m1+2)],usable)
        SL <- 1-pweibull(li,1/parest[(m1+3)],exp(s_lp))
        SR <- 1-pweibull(ri,1/parest[(m1+3)],exp(s_lp))
        mat2 <- matrix(1,usable,1)
      } else if (m2==1)
      {
        s_lp <- parest[(m1+2)]+mat2*parest[(m1+3)]
        SL <- 1-pweibull(li,1/parest[(m1+m2+3)],exp(s_lp))
        SR <- 1-pweibull(ri,1/parest[(m1+m2+3)],exp(s_lp))
        # create design matrices
        mat2 <- cbind(rep(1,usable),mat2)
      } else if (m2>1)
      {
        s_lp <- parest[(m1+2)]+mat2%*%parest[(m1+3):(m1+m2+2)]
        SL <- 1-pweibull(li,1/parest[(m1+m2+3)],exp(s_lp))
        SR <- 1-pweibull(ri,1/parest[(m1+m2+3)],exp(s_lp))
        # create design matrices
        mat2 <- cbind(rep(1,usable),mat2)
      }

      p <- exp(p_lp)/(1+exp(p_lp))
      tau <- parest[m1+m2+3]


      #hessian matrix
      h <- matrix(NA,m1+m2+3,m1+m2+3)

      for (i in 1:(m1+1))
      {
        for (j in i:(m1+1))
        {
          h[i,j] <- -1*sum(p/{1+exp(p_lp)}*mat1[,i]*mat1[,j])
          h[j,i] <- -1*sum(p/{1+exp(p_lp)}*mat1[,i]*mat1[,j])
        }
        for (j in 1:(m2+1))
        {
          h[i,m1+1+j] <- 0
          h[m1+1+j,i] <- 0
        }
        h[i,m1+m2+3] <- 0
        h[m1+m2+3,i] <- 0
      }

      for (i in 1:(m1+1))
      {
        for (j in i:(m1+1))
        {

          h_beta_a <- -1*sum((p/{1+exp(p_lp)}*mat1[,i]*mat1[,j])[k==1])
          h_beta_b <- sum({({exp(p_lp)*(1-SR)} / {(exp(p_lp)+1-SR)**{2}} - p / {1+exp(p_lp)})*mat1[,i]*mat1[,j]}[k==0])
          h[i,j] <- h_beta_a + h_beta_b
          h[j,i] <- h_beta_a + h_beta_b
        }

        for (j in 1:(m2+1))
        {
          h[i,m1+1+j] <- sum(({exp(p_lp)*ri**{1/tau}*SR} / {tau*exp(s_lp/tau)*(exp(p_lp)+1-SR)**2}*mat1[,i]*mat2[,j])[k==0&ri<Inf])
          h[m1+1+j,i] <- sum(({exp(p_lp)*ri**{1/tau}*SR} / {tau*exp(s_lp/tau)*(exp(p_lp)+1-SR)**2}*mat1[,i]*mat2[,j])[k==0&ri<Inf])
        }
        h[i,m1+m2+3] <- sum(({exp(p_lp)*(log(ri)-s_lp)*ri**{1/tau}*SR} / {tau**2*exp(s_lp/tau)*(exp(p_lp)+1-SR)**2}*mat1[,i])[k==0&ri<Inf])
        h[m1+m2+3,i] <- sum(({exp(p_lp)*(log(ri)-s_lp)*ri**{1/tau}*SR} / {tau**2*exp(s_lp/tau)*(exp(p_lp)+1-SR)**2}*mat1[,i])[k==0&ri<Inf])
      }

      for (i in 1:(m2+1))
      {
        for (j in i:(m2+1))
        {
          h_gamma_a1 <- sum(((1-ci)*{((li*exp(s_lp))**{1/tau}+(ri*exp(s_lp))**{1/tau}-(li**{1/tau}-ri**{1/tau})**2)*SL*SR-(li*exp(s_lp))**{1/tau}*SL**2-(ri*exp(s_lp))**{1/tau}*SR**2} / {(tau*exp(s_lp/tau)*(SL-SR))**2}*mat2[,i]*mat2[,j])[k==1 & ri<Inf & ci!=1])
          h_gamma_a2 <- sum(((1-ci)*{-(li*exp(s_lp))**{1/tau}*SL**2} / {(tau*exp(s_lp/tau)*(SL-SR))**2}*mat2[,i]*mat2[,j])[k==1 & ri==Inf])
          h_gamma_b <- sum(({((ri*exp(s_lp))**{1/tau}-ri**{2/tau})*SR*(1+exp(p_lp))-(ri*exp(s_lp))**{1/tau}*SR**2} / {(tau*exp(s_lp/tau)*(exp(p_lp)+1-SR))**2}*mat2[,i]*mat2[,j])[k==0&ri<Inf])
          h[m1+1+i,m1+1+j] <- h_gamma_a1 + h_gamma_a2 + h_gamma_b
          h[m1+1+j,m1+1+i] <- h_gamma_a1 + h_gamma_a2 + h_gamma_b
        }

        h_gamma_tau_a1 <- sum(((1-ci)*({((li*exp(s_lp))**{1/tau}+(li*ri)**{1/tau}-li**{2/tau})*(log(li)-s_lp)/tau + (li*exp(s_lp))**{1/tau} + ((ri*exp(s_lp))**{1/tau} + (li*ri)**{1/tau}-ri**{2/tau})*(log(ri)-s_lp)/tau+(ri*exp(s_lp))**{1/tau}}*SL*SR-(1+(log(li)-s_lp)/tau)*(li*exp(s_lp))**{1/tau}*SL**2-(1+(log(ri)-s_lp)/tau)*(ri*exp(s_lp))**{1/tau}*SR**2) / ({tau*exp(s_lp/tau)*(SL-SR)}**2)*mat2[,i])[k==1&li>0&ri<Inf])
        h_gamma_tau_a2 <- sum(((1-ci)*({((ri*exp(s_lp))**{1/tau}-ri**{2/tau})*(log(ri)-s_lp)/tau+(ri*exp(s_lp))**{1/tau}}*SL*SR-(1+(log(ri)-s_lp)/tau)*(ri*exp(s_lp))**{1/tau}*SR**2) / ({tau*exp(s_lp/tau)*(SL-SR)}**2)*mat2[,i])[k==1&li==0&ri<Inf & ci!=1])
        h_gamma_tau_a3 <- sum(((1-ci)*(-(1+(log(li)-s_lp)/tau)*(li*exp(s_lp))**{1/tau}*SL**2) / ({tau*exp(s_lp/tau)*(SL-SR)}**2)*mat2[,i])[k==1&li>0&ri==Inf])
        h_gamma_tau_a4 <- 0 #0 when li==0 & ri == Inf
        h_gamma_tau_b <- sum((((1+(log(ri)-s_lp)/tau)*(ri*exp(s_lp))**{1/tau}*SR*(exp(p_lp)+1-SR)-(log(ri)-s_lp)/tau*ri**{2/tau}*SR*(exp(p_lp)+1)) / {(tau*exp(s_lp/tau)*(exp(p_lp)+1-SR))**2}*mat2[,i])[k==0 & ri<Inf])
        h[m1+1+i,m1+m2+3] <- h_gamma_tau_a1 + h_gamma_tau_a2 + h_gamma_tau_a3 + h_gamma_tau_a4 + h_gamma_tau_b
        h[m1+m2+3,m1+1+i] <- h_gamma_tau_a1 + h_gamma_tau_a2 + h_gamma_tau_a3 + h_gamma_tau_a4 + h_gamma_tau_b
      }

      h_tau_a1 <- sum(((1-ci)*({(log(li)-s_lp+2*tau)*(log(li)-s_lp)*(li*exp(s_lp))**{1/tau} + (log(ri)-s_lp+2*tau)*(log(ri)-s_lp)*(ri*exp(s_lp))**{1/tau}-((log(li)-s_lp)*li**{1/tau}-(log(ri)-s_lp)*ri**{1/tau})**2}*SL*SR - (log(ri)-s_lp+2*tau)*(log(ri)-s_lp)*(ri*exp(s_lp))**{1/tau}*SR**2 - (log(li)-s_lp+2*tau)*(log(li)-s_lp)*(li*exp(s_lp))**{1/tau}*SL**2) / ({tau**2*exp(s_lp/tau)*(SL-SR)}**2) )[k==1&li>0&ri<Inf])
      h_tau_a2 <- sum(((1-ci)*({(log(ri)-s_lp+2*tau)*(log(ri)-s_lp)*(ri*exp(s_lp))**{1/tau} - ((log(ri)-s_lp)*ri**{1/tau})**2}*SR - (log(ri)-s_lp+2*tau)*(log(ri)-s_lp)*(ri*exp(s_lp))**{1/tau}*SR**2) / ({tau**2*exp(s_lp/tau)*(SL-SR)}**2) )[k==1&li==0&ri<Inf&ci!=1])
      h_tau_a3 <- sum(((1-ci)*(-(log(li)-s_lp+2*tau)*(log(li)-s_lp)*(li*exp(s_lp))**{1/tau}*SL**2) / ({tau**2*exp(s_lp/tau)*(SL-SR)}**2) )[k==1&li>0&ri==Inf])
      h_tau_a4 <- 0 #0 when li==0 & ri == Inf
      h_tau_b <- sum((((log(ri)-s_lp)*(log(ri)-s_lp+2*tau)*(ri*exp(s_lp))**{1/tau}*SR*(exp(p_lp)+1-SR)-(log(ri)-s_lp)**2*ri**{2/tau}*SR*(exp(p_lp)+1)) / ({tau**2*exp(s_lp/tau)*(exp(p_lp)+1-SR)}**2))[k==0 & ri<Inf])
      h[m1+m2+3,m1+m2+3] <- h_tau_a1 + h_tau_a2 + h_tau_a3 + h_tau_a4 + h_tau_b

      hessian <- -1*h

      return(hessian)
    }

    likfun <- function(parest)
    {
      k <- xsam[,1]
      ci <- xsam[,2]
      li <- xsam[,3]
      ri <- xsam[,4]
      k2 <- xsam[,5]

      if (m1==0) {p_lp <- rep(parest[1],usable)
      } else if (m1==1) {p_lp <- parest[1]+mat1*parest[2]
      } else if (m1>1) {p_lp <- parest[1]+mat1%*%parest[2:(m1+1)] }

      if (m2==0)
      {
        s_lp <- rep(parest[(m1+2)],usable)
        SL <- 1-pweibull(li,1/parest[(m1+3)],exp(s_lp))
        SR <- 1-pweibull(ri,1/parest[(m1+3)],exp(s_lp))
      } else if (m2==1)
      {
        s_lp <- parest[(m1+2)]+mat2*parest[(m1+3)]
        SL <- 1-pweibull(li,1/parest[(m1+m2+3)],exp(s_lp))
        SR <- 1-pweibull(ri,1/parest[(m1+m2+3)],exp(s_lp))
      } else if (m2>1)
      {
        s_lp <- parest[(m1+2)]+mat2%*%parest[(m1+3):(m1+m2+2)]
        SL <- 1-pweibull(li,1/parest[(m1+m2+3)],exp(s_lp))
        SR <- 1-pweibull(ri,1/parest[(m1+m2+3)],exp(s_lp))
      }

      #log-likelihood
      tmp <- rep(NA,usable)
      tmp[which(ci==1 & k==1)] <- I(p_lp-log(1+exp(p_lp)))[which(ci==1 & k==1)]
      tmp[which(ci==0 & k==1)] <- I(log(SL-SR)-log(1+exp(p_lp)))[which(ci==0 & k==1)]
      tmp[which(k==0 & k2==1)] <- I(log(exp(p_lp)+1-SR)-log(1+exp(p_lp)))[which(k==0 & k2==1)]
      tmp[which(k==0 & k2==0)] <- 0
      sum(tmp) #full log-likelihood
    }


    parest <- parinit
    print(parinit)
    parest_next <- c(rep(0,m1+m2+2),1)
    i <- 1
    while(max(abs(parest_next-parest))>.001){
      if (i>1) {parest <- parest_next}
      ci_prime <- estep(parest)
      parest_next <- mstep(parest,ci_prime)
      i = i+1
      #print(parest_next)
      #print(max(abs(parest_next-parest)))
    }
    fin_res <- list(parest_next,
                    covvar(parest_next),
                    likfun(parest_next))
    names(fin_res) <- c("par","hessian","loglikelihood")
    return(fin_res)
  }

  ####### Logistic-Exponential Functions #######
  #find initial parameter estimates
  logistic_exponential_init <- function()
  {
    # make crude assumptions to get initial parameter estimates
    xsam2 <- xsam
    if (n1>n2)
    {
      xsam2[which(is.na(xsam2[,2])),2] <- 1
    } else if (n1<=n2)
    {
      xsam2[which(is.na(xsam2[,2])),2] <- 0
    }
    xsam2[which(xsam2[,3]==0),3] <- NA
    xsam2[which(xsam2[,4]==Inf),4] <- NA

    if (m1==0)
    {
      logistic_res <- glm(xsam2[,2] ~ 1, family = "binomial")
      param_init <- logistic_res$coefficients[1]
      names(param_init) <- "logit-(Intercept)"
    } else if (m1>0)
    {
      logistic_res <- glm(xsam2[,2] ~ mat1, family = "binomial")
      param_init <- logistic_res$coefficients

      if (is.null(colnames(mat1))) {
        names(param_init)<- paste0("logit-", names(logistic_res$coefficients))
      } else{
        names(param_init)<- paste0("logit-", c("(Intercept)", colnames(mat1)))
      }
    }

    subset<- which(xsam2[,2]==0)
    if (m2==0)
    {
      exp_res <- survreg(Surv(xsam2[,3],xsam2[,4],type='interval2')~1, dist="weibull", scale=1, subset=subset)
      param_init2 <- exp_res$coefficients[1]
      names(param_init2) <- "Exponential-log(scale)"
    } else if (m2>0)
    {
      exp_res <- survreg(Surv(xsam2[,3],xsam2[,4],type='interval2')~ mat2, dist="weibull", scale=1, subset=subset)
      param_init2 <- exp_res$coefficients

      if(is.null(colnames(mat2))) {
        names(param_init2)<- paste0("Exponential-", names(exp_res$coefficients))
      } else {
        names(param_init2)<- paste0("Exponential-", c("(Intercept)", colnames(mat2)))
      }
    }

    param_init <- c(param_init,param_init2)
    return(param_init)
  }

  #full log-likelihood
  logistic_exponential_likfun <- function(parest)
  {
    k <- xsam[,1]
    ci <- xsam[,2]
    li <- xsam[,3]
    ri <- xsam[,4]
    k2 <- xsam[,5]

    if (m1==0) {p_lp <- rep(parest[1],usable)
    } else {p_lp <- parest[1]+mat1%*%parest[2:(m1+1)] }

    if (m2==0)
    {
      s_lp <- rep(parest[(m1+2)],usable)
      SL <- 1-pweibull(li,1,exp(s_lp))
      SR <- 1-pweibull(ri,1,exp(s_lp))
    } else if (m2>0)
    {
      s_lp <- parest[(m1+2)]+mat2%*%parest[(m1+3):(m1+m2+2)]
      SL <- 1-pweibull(li,1,exp(s_lp))
      SR <- 1-pweibull(ri,1,exp(s_lp))
    }

    #log-likelihood
    tmp <- rep(NA,usable)
    tmp[which(ci==1 & k==1)] <- I(p_lp-log(1+exp(p_lp)))[which(ci==1 & k==1)]
    tmp[which(ci==0 & k==1)] <- I(log(SL-SR)-log(1+exp(p_lp)))[which(ci==0 & k==1)]
    tmp[which(k==0 & k2==1)] <- I(log(exp(p_lp)+1-SR)-log(1+exp(p_lp)))[which(k==0 & k2==1)]
    tmp[which(k==0 & k2==0)] <- 0
    sum(tmp) #full log-likelihood
  }

  #minimize function to maximize likelihood
  logistic_exponential_opt <- function(y)
  {
    -1*logistic_exponential_likfun(y)
  }

  logistic_exponential_gradem <- function(parinit){
    #code EM algorithm
    estep <- function(parest){
      k <- xsam[,1]
      ci <- xsam[,2]
      ri <- xsam[,4]

      if (m1==0) {p_lp <- rep(parest[1],usable)
      } else if (m1>0) {p_lp <- parest[1]+mat1%*%parest[2:(m1+1)] }

      if (m2==0)
      {
        s_lp <- rep(parest[(m1+2)],usable)
        SR <- 1-pweibull(ri,1,exp(s_lp))
      } else if (m2>0)
      {
        s_lp <- parest[(m1+2)]+mat2%*%parest[(m1+3):(m1+m2+2)]
        SR <- 1-pweibull(ri,1,exp(s_lp))
      }

      ci_prime <- rep(NA,usable)
      ci_prime[k==1] <- ci[k==1]
      ci_prime[k==0] <- I(exp(p_lp)/(exp(p_lp)+1-SR))[k==0]
      return(ci_prime)
    }

    mstep <- function(parest,ci_prime){
      ci <- ci_prime
      li <- xsam[,3]
      ri <- xsam[,4]

      if (m1==0)
      {
        p_lp <- rep(parest[1],usable)
        # create design matrices
        mat1 <- matrix(1,usable,1)
      } else if (m1>0)
      {
        p_lp <- parest[1]+mat1%*%parest[2:(m1+1)]
        # create design matrices
        mat1 <- cbind(rep(1,usable),mat1)
      }

      if (m2==0)
      {
        s_lp <- rep(parest[(m1+2)],usable)
        SL <- 1-pweibull(li,1,exp(s_lp))
        SR <- 1-pweibull(ri,1,exp(s_lp))
        mat2 <- matrix(1,usable,1)
      } else if (m2>0)
      {
        s_lp <- parest[(m1+2)]+mat2%*%parest[(m1+3):(m1+m2+2)]
        SL <- 1-pweibull(li,1,exp(s_lp))
        SR <- 1-pweibull(ri,1,exp(s_lp))
        # create design matrices
        mat2 <- cbind(rep(1,usable),mat2)
      }

      p <- exp(p_lp)/(1+exp(p_lp))
      tau <- 1

      #score matrix
      u <- rep(NA,m1+m2+2)

      for (i in 1:(m1+1))
      {
        u[i] <- sum((ci-p)*mat1[,i])
      }

      for (i in 1:(m2+1))
      {
        u_gamma_a1 <- sum(((1-ci)*{li**{1/tau}*SL-ri**{1/tau}*SR} / {tau*exp(s_lp/tau)*(SL-SR)}*mat2[,i])[ri<Inf & ci!=1])
        u_gamma_a2 <- sum(((1-ci)*{li**{1/tau}*SL} / {tau*exp(s_lp/tau)*SL}*mat2[,i])[ri==Inf])
        u[m1+1+i] <- u_gamma_a1 + u_gamma_a2
      }


      #hessian matrix
      h <- matrix(NA,m1+m2+2,m1+m2+2)
      for (i in 1:(m1+1))
      {
        for (j in i:(m1+1))
        {
          h[i,j] <- -1*sum(p/{1+exp(p_lp)}*mat1[,i]*mat1[,j])
          h[j,i] <- -1*sum(p/{1+exp(p_lp)}*mat1[,i]*mat1[,j])
        }
        for (j in 1:(m2+1))
        {
          h[i,m1+1+j] <- 0
          h[m1+1+j,i] <- 0
        }
      }

      for (i in 1:(m2+1))
      {
        for (j in i:(m2+1))
        {
          h_gamma_a1 <- sum(((1-ci)*{((li*exp(s_lp))**{1/tau}+(ri*exp(s_lp))**{1/tau}-(li**{1/tau}-ri**{1/tau})**2)*SL*SR-(li*exp(s_lp))**{1/tau}*SL**2-(ri*exp(s_lp))**{1/tau}*SR**2} / {(tau*exp(s_lp/tau)*(SL-SR))**2}*mat2[,i]*mat2[,j])[ri<Inf & ci!=1])
          h_gamma_a2 <- sum(((1-ci)*{-(li*exp(s_lp))**{1/tau}*SL**2} / {(tau*exp(s_lp/tau)*(SL-SR))**2}*mat2[,i]*mat2[,j])[ri==Inf])
          h[m1+1+i,m1+1+j] <- h_gamma_a1 + h_gamma_a2
          h[m1+1+j,m1+1+i] <- h_gamma_a1 + h_gamma_a2
        }
      }

      hessian <- -1*h

      fisher_info <- solve(hessian)

      parest_next <- parest + fisher_info%*%u

      return(parest_next)
    }

    covvar <- function(parest){
      k <- xsam[,1]
      ci <- xsam[,2]
      li <- xsam[,3]
      ri <- xsam[,4]
      k2 <- xsam[,5]

      if (m1==0)
      {
        p_lp <- rep(parest[1],usable)
        # create design matrices
        mat1 <- matrix(1,usable,1)
      } else if (m1>0)
      {
        p_lp <- parest[1]+mat1%*%parest[2:(m1+1)]
        # create design matrices
        mat1 <- cbind(rep(1,usable),mat1)
      }

      if (m2==0)
      {
        s_lp <- rep(parest[(m1+2)],usable)
        SL <- 1-pweibull(li,1,exp(s_lp))
        SR <- 1-pweibull(ri,1,exp(s_lp))
        mat2 <- matrix(1,usable,1)
      } else if (m2>0)
      {
        s_lp <- parest[(m1+2)]+mat2%*%parest[(m1+3):(m1+m2+2)]
        SL <- 1-pweibull(li,1,exp(s_lp))
        SR <- 1-pweibull(ri,1,exp(s_lp))
        # create design matrices
        mat2 <- cbind(rep(1,usable),mat2)
      }

      p <- exp(p_lp)/(1+exp(p_lp))
      tau <- 1

      #hessian matrix
      h <- matrix(NA,m1+m2+2,m1+m2+2)

      for (i in 1:(m1+1))
      {
        for (j in i:(m1+1))
        {
          h_beta_a <- -1*sum((p/{1+exp(p_lp)}*mat1[,i]*mat1[,j])[k==1])
          h_beta_b <- sum(({exp(p_lp)*(1-SR)} / {(exp(p_lp)+1-SR)**{2}} - p / {1+exp(p_lp)}*mat1[,i]*mat1[,j])[k==0])
          h[i,j] <- h_beta_a + h_beta_b
          h[j,i] <- h_beta_a + h_beta_b
        }

        for (j in 1:(m2+1))
        {
          h[i,m1+1+j] <- sum(({exp(p_lp)*ri**{1/tau}*SR} / {tau*exp(s_lp/tau)*(exp(p_lp)+1-SR)**2}*mat1[,i]*mat2[,j])[k==0&ri<Inf])
          h[m1+1+j,i] <- sum(({exp(p_lp)*ri**{1/tau}*SR} / {tau*exp(s_lp/tau)*(exp(p_lp)+1-SR)**2}*mat1[,i]*mat2[,j])[k==0&ri<Inf])
        }
      }

      for (i in 1:(m2+1))
      {
        for (j in i:(m2+1))
        {
          h_gamma_a1 <- sum(((1-ci)*{((li*exp(s_lp))**{1/tau}+(ri*exp(s_lp))**{1/tau}-(li**{1/tau}-ri**{1/tau})**2)*SL*SR-(li*exp(s_lp))**{1/tau}*SL**2-(ri*exp(s_lp))**{1/tau}*SR**2} / {(tau*exp(s_lp/tau)*(SL-SR))**2}*mat2[,i]*mat2[,j])[k==1 & ri<Inf & ci!=1])
          h_gamma_a2 <- sum(((1-ci)*{-(li*exp(s_lp))**{1/tau}*SL**2} / {(tau*exp(s_lp/tau)*(SL-SR))**2}*mat2[,i]*mat2[,j])[k==1 & ri==Inf])
          h_gamma_b <- sum(({((ri*exp(s_lp))**{1/tau}-ri**{2/tau})*SR*(1+exp(p_lp))-(ri*exp(s_lp))**{1/tau}*SR**2} / {(tau*exp(s_lp/tau)*(exp(p_lp)+1-SR))**2}*mat2[,i]*mat2[,j])[k==0&ri<Inf])
          h[m1+1+i,m1+1+j] <- h_gamma_a1 + h_gamma_a2 + h_gamma_b
          h[m1+1+j,m1+1+i] <- h_gamma_a1 + h_gamma_a2 + h_gamma_b
        }
      }
      hessian <- -1*h

      return(hessian)
    }

    likfun <- function(parest)
    {
      k <- xsam[,1]
      ci <- xsam[,2]
      li <- xsam[,3]
      ri <- xsam[,4]
      k2 <- xsam[,5]

      if (m1==0) {p_lp <- rep(parest[1],usable)
      } else {p_lp <- parest[1]+mat1%*%parest[2:(m1+1)] }

      if (m2==0)
      {
        s_lp <- rep(parest[(m1+2)],usable)
        SL <- 1-pweibull(li,1,exp(s_lp))
        SR <- 1-pweibull(ri,1,exp(s_lp))
      } else if (m2>0)
      {
        s_lp <- parest[(m1+2)]+mat2%*%parest[(m1+3):(m1+m2+2)]
        SL <- 1-pweibull(li,1,exp(s_lp))
        SR <- 1-pweibull(ri,1,exp(s_lp))
      }

      #log-likelihood
      tmp <- rep(NA,usable)
      tmp[which(ci==1 & k==1)] <- I(p_lp-log(1+exp(p_lp)))[which(ci==1 & k==1)]
      tmp[which(ci==0 & k==1)] <- I(log(SL-SR)-log(1+exp(p_lp)))[which(ci==0 & k==1)]
      tmp[which(k==0 & k2==1)] <- I(log(exp(p_lp)+1-SR)-log(1+exp(p_lp)))[which(k==0 & k2==1)]
      tmp[which(k==0 & k2==0)] <- 0
      sum(tmp) #full log-likelihood
    }


    parest <- parinit
    print(parinit)
    parest_next <- rep(0,m1+m2+2)
    i <- 1
    while(max(abs(parest_next-parest))>.00001){
      if (i>1) {parest <- parest_next}
      ci_prime <- estep(parest)
      parest_next <- mstep(parest,ci_prime)
      i = i+1
      #print(parest_next)
      #print(max(abs(parest_next-parest)))
    }
    fin_res <- list(parest_next,
                    covvar(parest_next),
                    likfun(parest_next))
    names(fin_res) <- c("par","hessian","loglikelihood")
    return(fin_res)
  }


  ####### Logistic-Loglogistic Functions #######
  #find initial parameter estimates
  pllog <- function (q, shape = 1, scale = 1, lower.tail = TRUE, log.p = FALSE)
  {
    Fx <- plogis(log(q), location = scale, scale = shape)
    if (!lower.tail)
      Fx <- 1 - Fx
    if (log.p)
      Fx <- log(Fx)
    return(Fx)
  }

  logistic_loglogistic_init <- function()
  {
    # make crude assumptions to get initial parameter estimates
    xsam2 <- xsam
    if (n1>n2)
    {
      xsam2[which(is.na(xsam2[,2])),2] <- 1
    } else if (n1<=n2)
    {
      xsam2[which(is.na(xsam2[,2])),2] <- 0
    }
    xsam2[which(xsam2[,3]==0),3] <- NA
    xsam2[which(xsam2[,4]==Inf),4] <- NA

    if (m1==0)
    {
      logistic_res <- glm(xsam2[,2] ~ 1, family = "binomial")
      param_init <- logistic_res$coefficients[1]
      names(param_init) <- "logit-(Intercept)"
    } else if (m1>0)
    {
      logistic_res <- glm(xsam2[,2] ~ mat1, family = "binomial")
      param_init <- logistic_res$coefficients

      if(is.null(colnames(mat1))) {
        names(param_init) <- paste0("logit-",names(logistic_res$coefficients))
      } else {
        names(param_init) <- paste0("logit-",c("(Intercept)", colnames(mat1)))
      }
    }

    subset<- which(xsam2[,2]==0)
    if (m2==0)
    {
      llog_res <- survreg(Surv(xsam2[,3],xsam2[,4],type='interval2')~1, dist="loglogistic", subset=subset)
      param_init2 <- c(llog_res$coefficients[1], llog_res$scale)
      names(param_init2) <- c("Loglogistic-scale","Loglogistic-shape")
    } else if (m2>0)
    {
      llog_res <- survreg(Surv(xsam2[,3],xsam2[,4],type='interval2')~mat2, dist="loglogistic", subset=subset)
      param_init2 <- c(llog_res$coefficients, llog_res$scale)

      if(is.null(colnames(mat2))) {
        names(param_init2)<- paste0("Loglogistic-", c(names(llog_res$coefficients),  "shape"))
      } else {
        names(param_init2)<- paste0("Loglogistic-", c("(Intercept)", colnames(mat2), "shape"))
      }
    }

    param_init <- c(param_init,param_init2)
    return(param_init)
  }

  #full log-likelihood
  logistic_loglogistic_likfun <- function(parest)
  {
    k <- xsam[,1]
    ci <- xsam[,2]
    li <- xsam[,3]
    ri <- xsam[,4]
    k2 <- xsam[,5]

    if (m1==0) {p_lp <- rep(parest[1],usable)
    } else {p_lp <- parest[1]+mat1%*%parest[2:(m1+1)] }

    if (m2==0)
    {
      s_lp <- rep(parest[(m1+2)],usable)
      SL <- 1-pllog(li,parest[(m1+3)],s_lp)
      SR <- 1-pllog(ri,parest[(m1+3)],s_lp)
    } else if (m2>0)
    {
      s_lp <- parest[(m1+2)]+mat2%*%parest[(m1+3):(m1+m2+2)]
      SL <- 1-pllog(li,parest[(m1+m2+3)],s_lp)
      SR <- 1-pllog(ri,parest[(m1+m2+3)],s_lp)
    }

    #log-likelihood
    tmp <- rep(NA,usable)
    tmp[which(ci==1 & k==1)] <- I(p_lp-log(1+exp(p_lp)))[which(ci==1 & k==1)]
    tmp[which(ci==0 & k==1)] <- I(log(SL-SR)-log(1+exp(p_lp)))[which(ci==0 & k==1)]
    tmp[which(k==0 & k2==1)] <- I(log(exp(p_lp)+1-SR)-log(1+exp(p_lp)))[which(k==0 & k2==1)]
    tmp[which(k==0 & k2==0)] <- 0
    sum(tmp) #full log-likelihood
  }

  #minimize function to maximize likelihood
  logistic_loglogistic_opt <- function(y)
  {
    -1*logistic_loglogistic_likfun(y)
  }

  ####### Logistic-Lognormal Functions #######
  logistic_lognormal_init <- function()
  {
    # make crude assumptions to get initial parameter estimates
    xsam2 <- xsam
    if (n1>n2)
    {
      xsam2[which(is.na(xsam2[,2])),2] <- 1
    } else if (n1<=n2)
    {
      xsam2[which(is.na(xsam2[,2])),2] <- 0
    }
    xsam2[which(xsam2[,3]==0),3] <- NA
    xsam2[which(xsam2[,4]==Inf),4] <- NA

    if (m1==0)
    {
      logistic_res <- glm(xsam2[,2] ~ 1, family = "binomial")
      param_init <- logistic_res$coefficients[1]
      names(param_init) <- "logit-(Intercept)"
    } else if (m1>0)
    {
      logistic_res <- glm(xsam2[,2] ~ mat1, family = "binomial")
      param_init <- logistic_res$coefficients

      if(is.null(colnames(mat1))) {
        names(param_init) <- paste0("logit-",names(logistic_res$coefficients))
      } else {
        names(param_init) <- paste0("logit-",c("(Intercept)", colnames(mat1)))
      }
    }

    subset<- which(xsam2[,2]==0)
    if (m2==0)
    {
      lnorm_res <- survreg(Surv(xsam2[,3],xsam2[,4],type='interval2')~1, dist="lognormal", subset=subset)
      param_init2 <- c(lnorm_res$coefficients[1], lnorm_res$scale)
      names(param_init2) <- c("lognormal-meanlog","lognormal-sdlog")
    } else if (m2>0)
    {
      lnorm_res <- survreg(Surv(xsam2[,3],xsam2[,4],type='interval2')~mat2, dist="lognormal", subset=subset)
      param_init2 <- c(lnorm_res$coefficients, lnorm_res$scale)

      if(is.null(colnames(mat2))) {
        names(param_init2)<- paste0("lognormal-", c(names(llog_res$coefficients), "sdlog"))
      } else {
        names(param_init2)<- paste0("lognormal-", c("(Intercept)", colnames(mat2),"sdlog"))
      }
    }

    param_init <- c(param_init,param_init2)
    return(param_init)

  }

  #full log-likelihood
  logistic_lognormal_likfun <- function(parest)
  {
    k <- xsam[,1]
    ci <- xsam[,2]
    li <- xsam[,3]
    ri <- xsam[,4]
    k2 <- xsam[,5]

    if (m1==0) {p_lp <- rep(parest[1],usable)
    } else {p_lp <- parest[1]+mat1%*%parest[2:(m1+1)] }

    if (m2==0)
    {
      s_lp <- rep(parest[(m1+2)],usable)
      SL <- 1-plnorm(li,s_lp,parest[(m1+3)])
      SR <- 1-plnorm(ri,s_lp,parest[(m1+3)])
    } else if (m2>0)
    {
      s_lp <- parest[(m1+2)]+mat2%*%parest[(m1+3):(m1+m2+2)]
      SL <- 1-plnorm(li,s_lp,parest[(m1+m2+3)])
      SR <- 1-plnorm(ri,s_lp,parest[(m1+m2+3)])
    }

    #log-likelihood
    tmp <- rep(NA,usable)
    tmp[which(ci==1 & k==1)] <- I(p_lp-log(1+exp(p_lp)))[which(ci==1 & k==1)]
    tmp[which(ci==0 & k==1)] <- I(log(SL-SR)-log(1+exp(p_lp)))[which(ci==0 & k==1)]
    tmp[which(k==0 & k2==1)] <- I(log(exp(p_lp)+1-SR)-log(1+exp(p_lp)))[which(k==0 & k2==1)]
    tmp[which(k==0 & k2==0)] <- 0
    sum(tmp) #full log-likelihood
  }

  #minimize function to maximize likelihood
  logistic_lognormal_opt <- function(y)
  {
    -1*logistic_lognormal_likfun(y)
  }

  ####### Logistic-Generalized Gamma Functions #######
  #find initial parameter estimates
  logistic_gengamma_init <- function()
  {
    # make crude assumptions to get initial parameter estimates
    xsam2 <- xsam
    if (n1>n2)
    {
      xsam2[which(is.na(xsam2[,2])),2] <- 1
    } else if (n1<=n2)
    {
      xsam2[which(is.na(xsam2[,2])),2] <- 0
    }
    xsam2[which(xsam2[,3]==0),3] <- NA
    xsam2[which(xsam2[,4]==Inf),4] <- NA

    if (m1==0)
    {
      logistic_res <- glm(xsam2[,2] ~ 1, family = "binomial")
      param_init <- logistic_res$coefficients[1]
      names(param_init) <- "logit-(Intercept)"
    } else if (m1>0) {
      logistic_res <- glm(xsam2[,2] ~ mat1, family = "binomial")
      param_init <- logistic_res$coefficients
      if (is.null(colnames(mat1))) {
        names(param_init)<- paste0("logit-", names(logistic_res$coefficients))
      } else{
        names(param_init)<- paste0("logit-", c("(Intercept)", colnames(mat1)))
      }
    }

    subset <- which(xsam2[,2]==0)
    if (m2==0)
    {
      ggamma_res <- survreg(Surv(xsam2[,3],xsam2[,4],type='interval2')~1, dist="weibull", subset=subset)
      param_init2 <- c(ggamma_res$coefficients[1], ggamma_res$scale,1)
      names(param_init2) <- c("Generalized Gamma-mu","Generalized Gamma-sigma","Generalized Gamma-Q")
    } else if (m2>0)
    {
      ggamma_res <- survreg(Surv(xsam2[,3],xsam2[,4],type='interval2')~ mat2, dist="weibull",subset=subset)
      param_init2 <- c(ggamma_res$coefficients, ggamma_res$scale,1)
      if(is.null(colnames(mat2))) {
        names(param_init2) <- c(paste0("Generalized Gamma-", names(ggamma_res$coefficients)),
                                "Generalized Gamma-sigma","Generalized Gamma-Q")
      } else {
        names(param_init2) <- c(paste0("Generalized Gamma-",c("(Intercept)", colnames(mat1))),
                                "Generalized Gamma-sigma", "Generalized Gamma-Q")
      }
    }

    param_init <- c(param_init,param_init2)
    return(param_init)
  }

  #full log-likelihood
  logistic_gengamma_likfun <- function(parest)
  {
    k <- xsam[,1]
    ci <- xsam[,2]
    li <- xsam[,3]
    ri <- xsam[,4]
    k2 <- xsam[,5]

    if (m1==0) {p_lp <- rep(parest[1],usable)
    } else {p_lp <- parest[1]+mat1%*%parest[2:(m1+1)] }

    if (m2==0)
    {
      s_lp <- rep(parest[(m1+2)],usable)
      SL <- 1-pgengamma(li,s_lp,parest[(m1+3)],parest[(m1+4)])
      SR <- 1-pgengamma(ri,s_lp,parest[(m1+3)],parest[(m1+4)])

    } else if (m2>0)
    {
      s_lp <- parest[(m1+2)]+mat2%*%parest[(m1+3):(m1+m2+2)]
      SL <- 1-pgengamma(li,s_lp,parest[(m1+m2+3)],parest[(m1+m2+4)])
      SR <- 1-pgengamma(ri,s_lp,parest[(m1+m2+3)],parest[(m1+m2+4)])
    }

    #log-likelihood
    tmp <- rep(NA,usable)
    tmp[which(ci==1 & k==1)] <- I(p_lp-log(1+exp(p_lp)))[which(ci==1 & k==1)]
    tmp[which(ci==0 & k==1)] <- I(log(SL-SR)-log(1+exp(p_lp)))[which(ci==0 & k==1)]
    tmp[which(k==0 & k2==1)] <- I(log(exp(p_lp)+1-SR)-log(1+exp(p_lp)))[which(k==0 & k2==1)]
    tmp[which(k==0 & k2==0)] <- 0
    sum(tmp) #full log-likelihood
  }

  #minimize function to maximize likelihood
  logistic_gengamma_opt <- function(y)
  {
    -1*logistic_gengamma_likfun(y)
  }

  ####### Logistic-Gamma Functions #######
  #find initial parameter estimates
  logistic_gamma_init <- function()
  {
    # make crude assumptions to get initial parameter estimates
    xsam2 <- xsam
    if (n1>n2)
    {
      xsam2[which(is.na(xsam2[,2])),2] <- 1
    } else if (n1<=n2)
    {
      xsam2[which(is.na(xsam2[,2])),2] <- 0
    }
    xsam2[which(xsam2[,3]==0),3] <- NA
    xsam2[which(xsam2[,4]==Inf),4] <- NA

    if (m1==0)
    {
      logistic_res <- glm(xsam2[,2] ~ 1, family = "binomial")
      param_init <- logistic_res$coefficients[1]
      names(param_init) <- "logit-(Intercept)"
    } else if (m1>0)
    {
      logistic_res <- glm(xsam2[,2] ~ mat1, family = "binomial")
      param_init <- logistic_res$coefficients
      if (is.null(colnames(mat1))) {
        names(param_init)<- paste("logit-", names(logistic_res$coefficients), sep="")
      } else{
        names(param_init)<- paste("logit-", c("(Intercept)", colnames(mat1)), sep="")
      }
    }

    subset<- which(xsam2[,2]==0)

    if (m2==0)
    {
      gamma_res <- survreg(Surv(xsam2[,3],xsam2[,4],type='interval2')~1, dist="weibull", subset=subset)
      param_init2 <- c(gamma_res$coefficients[1], gamma_res$scale**(-2))
      names(param_init2) <- c("Gamma-mu","Gamma-shape")
    } else if (m2>0)
    {
      gamma_res <- survreg(Surv(xsam2[,3],xsam2[,4],type='interval2')~ mat2, dist="weibull", subset=subset)
      param_init2 <- c(gamma_res$coefficients, gamma_res$scale**(-2))

      if(is.null(colnames(mat2))) {
        names(param_init2) <- c(paste0("Gamma-",names(gamma_res$coefficients)), "Gamma-shape")
      } else {
        names(param_init2) <-c(paste0("Gamma-",c("Intercept", colnames(mat2))),"Gamma-shape")
      }
    }

    param_init <- c(param_init,param_init2)
    return(param_init)
  }

  #full log-likelihood
  logistic_gamma_likfun <- function(parest)
  {
    k <- xsam[,1]
    ci <- xsam[,2]
    li <- xsam[,3]
    ri <- xsam[,4]
    k2 <- xsam[,5]

    if (m1==0) {p_lp <- rep(parest[1],usable)
    } else {p_lp <- parest[1]+mat1%*%parest[2:(m1+1)] }

    if (m2==0)
    {
      s_lp <- rep(parest[(m1+2)],usable)
      SL <- 1-pgamma(li,parest[(m1+3)],exp(-1*s_lp)*parest[(m1+3)])
      SR <- 1-pgamma(ri,parest[(m1+3)],exp(-1*s_lp)*parest[(m1+3)])
    } else if (m2>0)
    {
      s_lp <- parest[(m1+2)]+mat2%*%parest[(m1+3):(m1+m2+2)]
      SL <- 1-pgamma(li,parest[(m1+m2+3)],exp(-1*s_lp)*parest[(m1+m2+3)])
      SR <- 1-pgamma(ri,parest[(m1+m2+3)],exp(-1*s_lp)*parest[(m1+m2+3)])
    }

    #log-likelihood
    tmp <- rep(NA,usable)
    tmp[which(ci==1 & k==1)] <- I(p_lp-log(1+exp(p_lp)))[which(ci==1 & k==1)]
    tmp[which(ci==0 & k==1)] <- I(log(SL-SR)-log(1+exp(p_lp)))[which(ci==0 & k==1)]
    tmp[which(k==0 & k2==1)] <- I(log(exp(p_lp)+1-SR)-log(1+exp(p_lp)))[which(k==0 & k2==1)]
    tmp[which(k==0 & k2==0)] <- 0
    sum(tmp) #full log-likelihood
  }

  #minimize function to maximize likelihood
  logistic_gamma_opt <- function(y)
  {
    -1*logistic_gamma_likfun(y)
  }
  ###############################
  ####### END FUNCTIONS #######
  ###############################


  ########################################################
  ## calculate the odds ratio from glm model
  ## coef: coefficients from GLM model
  ## se: standard error of the coef
  ##
  ##
  ########################################################

  get.OR<- function(coef, se, level=0.05) {
    qa<- qnorm(1 - level/2)
    log.or<- coef
    se.or<- exp(log.or)*se
    LL<- exp(log.or - qa*se)
    UL<- exp(log.or + qa*se)

    ##modified to allow using character "-" in coding nominal variables #######

    sep.names<-strsplit(names(log.or), "-")
    name.mat<-plyr::ldply(sep.names, rbind)

    modname<-as.character(name.mat[,1])
    if(ncol(name.mat)>2){
      label2<-apply(name.mat[,-1],1,function(x) paste(x,collapse = "-"))
      label3<-sub("-NA","",label2)
      parname<-label3
    }else{
      parname<-as.character(name.mat[,2])
    }




    #modname <- matrix(unlist(strsplit(names(log.or),"-")),2)[1,]
    #parname <- matrix(unlist(strsplit(names(log.or),"-")),2)[2,]
    or.sum <- data.frame(modname,parname,exp(log.or),se.or,LL,UL)
    colnames(or.sum) <- c("Model","Label","OR","SE","95%LL","95%UL")
    rownames(or.sum) <- NULL
    return(or.sum)
  }

  ########################################################
  ## calculate the HR ratio from the AFT model
  ## mainly for Weibull and exponential survival model
  ## gamma: the coefficients for design matrix
  ## tau:   the scale parameter, for exponential survial, it will be 1
  ##
  ## cov.mat: covariance matrix between gamma and tau, if tau !=1
  ##          if tau==1, then it will be covariance of gamma
  ##
  ########################################################

  get.HR<- function(gamma, tau, cov.mat, level=0.05) {

    if(tau ==1 ) {
      result<- get.OR(-gamma, sqrt(diag(cov.mat))[1:length(gamma)])
      colnames(result)<- c("Model","Label","exp(-coef/scale)", "SE", "95%LL", "95%UL")
      return(result)
    } else {
      qa<-    qnorm(1 - level/2)
      k<-     length(gamma)

      var.tau<- cov.mat[k+1, k+1]
      var.gamma<- diag(as.matrix(cov.mat[1:k, 1:k]))
      cov.gamma.tau<- cov.mat[1:k, k+1]

      log.HR<- -gamma/tau
      HR<- exp(-gamma/tau)
      var.log.HR<- gamma^2/tau^2*(var.gamma/gamma^2 + var.tau/tau^2 - 2*cov.gamma.tau/(gamma*tau))
      se.log.HR<- sqrt(var.log.HR)

      SE.HR<- sqrt(var.log.HR)*HR

      LL<- exp(log.HR - qa*se.log.HR)
      UL<- exp(log.HR + qa*se.log.HR)


      sep.names<-strsplit(names(HR),"-")
      name.mat<-plyr::ldply(sep.names, rbind)
      modname<-as.character(name.mat[,1])

      if(ncol(name.mat)>2){
        label2<-apply(name.mat[,-1],1,function(x) paste(x,collapse = "-"))
        label3<-sub("-NA","",label2)
        parname<-label3
      }else {
        parname<-as.character(name.mat[,2])
      }




      #modname <- matrix(unlist(strsplit(names(HR),"-")),2)[1,]
      #parname <- matrix(unlist(strsplit(names(HR),"-")),2)[2,]
      hr.sum <- data.frame(modname,parname,HR,SE.HR,LL,UL)
      colnames(hr.sum) <- c("Model","Label","HR","SE","95%LL","95%UL")
      rownames(hr.sum) <- NULL
      return(hr.sum)
    }
  }


  res_optim<- list()
  if (model=="non-parametric") {
    ####### Non-parametric #######
    li2 <- ifelse(xsam[,1]==0,0,
                  ifelse(is.na(xsam[,2])==0 & xsam[,2]==1,0,
                         ifelse(xsam[,3]==0,0.01,xsam[,3])))
    ri2 <- ifelse(xsam[,1]==0,xsam[,4],
                  ifelse(is.na(xsam[,2])==0 & xsam[,2]==1,.01,xsam[,4]))
    if (conf.int=="FALSE")
    {
      nonparam <- icfit(li2,ri2,conf.int=FALSE)
    } else if (conf.int=="TRUE")
    {
      nonparam <- icfit(li2,ri2,conf.int=TRUE)
    }

    nonparam$model <- model
    nonparam$conf.int <- conf.int
    return(nonparam)
  } else if (model=="logistic-Weibull")  ####### logistic-Weibull #######
  {
    ### Special cases ###
    ### No prevalent disease ###

    qa<- qnorm(0.975)
    if (n1==0 & n3==0)
    {
      print ("Incident disease only")
      xsam2 <- xsam
      xsam2[which(xsam2[,3]==0),3] <- NA
      xsam2[which(xsam2[,4]==Inf),4] <- NA

      if (m2==0)
      {
        weib_res <- survreg(Surv(xsam2[,3],xsam2[,4],type='interval2')~1, dist="weibull")
        par <- c(weib_res$coefficients[1], weib_res$scale)
        names(par) <- c("Weibull-log(scale)","Weibull-1/shape")
      } else if (m2>0)
      {
        weib_res <-survreg(Surv(xsam2[,3],xsam2[,4],type='interval2')~ mat2, dist="weibull")
        par <- c(weib_res$coefficients, weib_res$scale)

        if(is.null(colnames(mat2))) {
          names(par)<- paste0("Weibull-", c(names(weib_res$coefficients), "1/shape"))
        } else {
          names(par)<- paste0("Weibull-", c("(Intercept)", colnames(mat2),"1/shape"))
        }
      }

      ### ADD HR FUNCTION HERE
      HR<- get.HR(weib_res$coefficients, weib_res$scale, vcov(weib_res))
      OR<- NULL
      loglikelihood <- weib_res$loglik[1]
      hessian <- solve(weib_res$var)

      res_optim <- list(par,HR,OR,loglikelihood,hessian)
      names(res_optim) <- c("par", "HR", "OR", "loglikelihood","hessian")
    } else if (n2==0 & n3==0)   ### No incident disease ###
    {
      print ("Prevalent disease only")
      xsam2 <- xsam
      xsam2[which(is.na(xsam2[,2])),2] <- 1

      if (m1==0)
      {
        logistic_res <- glm(xsam2[,2] ~ 1, family = "binomial")
        par <- logistic_res$coefficients[1]
        names(par) <- "logit-(Intercept)"
        hessian <- 1/summary(logistic_res)$coefficients[2]^2

      } else if (m1>0)
      {
        logistic_res <- glm(xsam2[,2] ~ mat1, family = "binomial")
        par <- logistic_res$coefficients

        if(is.null(colnames(mat1))) {
          names(par)<- paste0("logit-", names(weib_res$coefficients))
        } else {
          names(par)<- paste0("logit-", c("(Intercept)", colnames(mat1)))
        }

        hessian <- solve(vcov(logistic_res))
      }

      ### ADD OR FUNCTION HERE
      HR<- NULL
      OR<- get.OR(par, sqrt(diag(vcov(logistic_res))))

      loglikelihood <- as.numeric(logLik(logistic_res))
      res_optim <- list(par,HR,OR,loglikelihood,hessian)
      names(res_optim) <- c("par","HR", "OR","loglikelihood","hessian")
    } else
    {
      if (is.null(init)) {
        init_param <- logistic_weibull_init()
      } else
      {
        init_param <- init
      }
      res_optim <- tryCatch(optim(init_param, logistic_weibull_opt, method="L-BFGS-B", lower=c(rep(-Inf,m1+m2+2),0.01), upper=c(rep(Inf,m1+m2+2),100), control=list(trace=TRUE), hessian=TRUE),error=function(e) NULL)
      if (class(res_optim)=="list"){
        names(res_optim)[names(res_optim)=="value"] <- "loglikelihood"
        res_optim$loglikelihood <- -1*res_optim$loglikelihood
      } else
      {
        res_optim <- logistic_weibull_gradem(init_param)
        names(res_optim$par)<- names(init_param)
      }

      ### ADD OR and HR FUNCTION HERE USING res_optim$par
      ### Cov-Var: solve(res_optim$hessian)
      ### First get the coefficients and SE from the logstic model
      ### Next  get the coefficients and SE from the weibull model

      index.logit.par  <- 1:(m1+1)
      index.weibull.par<- (m1+2):length(res_optim$par)

      logit.par<-   res_optim$par[index.logit.par]
      weibull.par<- res_optim$par[index.weibull.par]

      cov.mat<-     solve(res_optim$hessian)
      logit.cov<-   matrix(cov.mat[index.logit.par, index.logit.par], nrow=length(index.logit.par))

      #############################################################################
      #
      # For weibull.par, the first n.gamma parameters are gamma_i, and the last
      # parameter is tau. To apply the function get.HR() to get HR, we need log(tau)
      #
      # Covariance matrix includes the covariances for gamma and tau
      #
      # So we need to transfer the covariance matrix into the matrix of log(tau) and gamma
      #############################################################################

      weibull.cov<- matrix(cov.mat[index.weibull.par, index.weibull.par], nrow=length(index.weibull.par))
      n.gamma<- length(weibull.par)-1

      #weibull.cov[n.gamma+1, n.gamma+1]<- 1/(weibull.par[n.gamma+1])^2 * weibull.cov[n.gamma+1, n.gamma+1]
      #weibull.cov[1:n.gamma, n.gamma+1]<- 1/weibull.par[n.gamma+1] * weibull.cov[1:n.gamma, n.gamma+1]
      #weibull.cov[n.gamma+1, 1:n.gamma]<- 1/weibull.par[n.gamma+1] * weibull.cov[n.gamma+1, 1:n.gamma]


      HR<- get.HR(weibull.par[1:n.gamma], weibull.par[n.gamma+1], weibull.cov)
      OR<- get.OR(logit.par, sqrt(diag(logit.cov)))

      res_optim$HR<- HR
      res_optim$OR<- OR

      res_optim$convergence<- list("convergence"=res_optim$convergence, "counts"=res_optim$counts, "message"=res_optim$message)
    }
  } else if (model=="logistic-exponential")  ####### logistic-Exponential #######
  {
    ### Special cases ###
    ### No prevalent disease ###
    if (n1==0 & n3==0)
    {
      print ("Incident disease only")
      xsam2 <- xsam
      xsam2[which(xsam2[,3]==0),3] <- NA
      xsam2[which(xsam2[,4]==Inf),4] <- NA

      if (m2==0)
      {
        exp_res <- survreg(Surv(xsam2[,3],xsam2[,4],type='interval2')~1, dist="weibull",scale=1)
        par <- exp_res$coefficients[1]
        names(par) <- "Exponential-log(scale)"
        hessian <- 1/(exp_res$var[1,1])
      } else if (m2>0)
      {
        exp_res <- survreg(Surv(xsam2[,3],xsam2[,4],type='interval2')~mat2, dist="weibull",scale=1)
        par <- exp_res$coefficients

        if(is.null(colnames(mat2))) {
          names(par)<- paste0("Exponential-", names(weib_res$coefficients))
        } else {
          names(par)<- paste0("Exponential-", c("(Intercept)", colnames(mat2)))
        }

        hessian <- solve(exp_res$var[1:(m2+1),1:(m2+1)])
      }

      loglikelihood <- exp_res$loglik[1]
      HR<- get.HR(par,1, exp_res$var)
      OR<- NULL
      res_optim <- list(par,HR,OR,loglikelihood,hessian)
      names(res_optim) <- c("par", "HR", "OR", "loglikelihood","hessian")
    } else if (n2==0 & n3==0)   ### No incident disease ###
    {
      print ("Prevalent disease only")
      xsam2 <- xsam
      xsam2[which(is.na(xsam2[,2])),2] <- 1

      if (m1==0)
      {
        logistic_res <- glm(xsam2[,2] ~ 1, family = "binomial")
        par <- logistic_res$coefficients[1]
        names(par) <- "logit-(Intercept)"
        hessian <- 1/summary(logistic_res)$coefficients[2]^2
      } else if (m1>0)
      {
        logistic_res <- glm(xsam2[,2] ~ mat1, family = "binomial")
        par <- logistic_res$coefficients

        if(is.null(colnames(mat1))) {
          names(par) <- paste0("logit-",names(logistic_res$coefficients))
        } else {
          names(par) <- paste0("logit-", c("(Intercept)", colnames(mat1)))
        }

        hessian <- solve(vcov(logistic_res))
      }

      loglikelihood <- as.numeric(logLik(logistic_res))
      HR<- NULL
      OR<- get.OR(par, sqrt(diag(vcov(logistic_res))))
      res_optim <- list(par, HR, OR, loglikelihood,hessian)
      names(res_optim) <- c("par","HR", "OR", "loglikelihood","hessian")
    } else
    {
      if (is.null(init)) {
        init_param <- logistic_exponential_init()
      } else
      {
        init_param <- init
      }
      res_optim <- tryCatch(optim(init_param, logistic_exponential_opt, method="BFGS", control=list(trace=TRUE), hessian=TRUE),error=function(e) NULL)
      if (class(res_optim)=="list"){
        names(res_optim)[names(res_optim)=="value"] <- "loglikelihood"
        res_optim$loglikelihood <- -1*res_optim$loglikelihood
      } else
      {
        res_optim <- logistic_exponential_gradem(init_param)
      }

      index.logit.par  <- 1:(m1+1)
      index.exp.par<- (m1+2):length(res_optim$par)

      logit.par<-   res_optim$par[index.logit.par]
      exp.par<-     res_optim$par[index.exp.par]

      cov.mat<-     solve(res_optim$hessian)
      logit.cov<-   matrix(cov.mat[index.logit.par, index.logit.par], nrow=length(index.logit.par))
      exp.cov<-     matrix(cov.mat[index.exp.par,   index.exp.par], nrow=length(index.exp.par))

      HR<- get.HR(exp.par, 1, exp.cov)
      OR<- get.OR(logit.par, sqrt(diag(logit.cov)))

      res_optim$HR<- HR
      res_optim$OR<- OR
      res_optim$convergence<- list("convergence"=res_optim$convergence, "counts"=res_optim$counts, "message"=res_optim$message)

    }
  } else if (model=="logistic-loglogistic")
  {

    ### Special cases ###
    ### No prevalent disease ###
    if (n1==0 & n3==0)
    {
      print ("Incident disease only")
      xsam2 <- xsam
      xsam2[which(xsam2[,3]==0),3] <- NA
      xsam2[which(xsam2[,4]==Inf),4] <- NA

      if (m2==0)
      {
        llog_res <- survreg(Surv(xsam2[,3],xsam2[,4],type='interval2')~1, dist="loglogistic")
        par <- c(llog_res$coefficients[1], llog_res$scale)
        names(par) <- c("Loglogistic-scale","Loglogistic-shape")
      } else if (m2>0)
      {
        llog_res <- survreg(Surv(xsam2[,3],xsam2[,4],type='interval2')~mat2, dist="loglogistic")
        par <- c(llog_res$coefficients, llog_res$scale)

        if(is.null(colnames(mat2))) {
          names(par) <- c(paste0("Loglogistic-",names(llog_res$coefficients)),  "Loglogistic-shape")
        } else {
          names(par) <- c(paste0("Loglogistic-",c("Intercept", colnames(mat2))),"Loglogistic-shape")
        }
      }
      loglikelihood <- llog_res$loglik[1]
      hessian <- solve(llog_res$var)
      res_optim <- list(par,loglikelihood,hessian)
      names(res_optim) <- c("par","loglikelihood","hessian")
    } else if (n2==0 & n3==0)   ### No incident disease ###
    {
      print ("Prevalent disease only")
      xsam2 <- xsam
      xsam2[which(is.na(xsam2[,2])),2] <- 1

      if (m1==0)
      {
        logistic_res <- glm(xsam2[,2] ~ 1, family = "binomial")
        par <- logistic_res$coefficients[1]
        names(par) <- "logit-Intercept"
        hessian <- 1/summary(logistic_res)$coefficients[2]^2
      } else if (m1>0)
      {
        logistic_res <- glm(xsam2[,2] ~ mat1, family = "binomial")
        par <- logistic_res$coefficients

        if(is.null(colnames(mat1))) {
          names(par) <- paste0("logit-",names(logistic_res$coefficients))
        } else {
          names(par) <- paste0("logit-", c("(Intercept)", colnames(mat1)))
        }

        hessian <- solve(vcov(logistic_res))
      }
      loglikelihood <- as.numeric(logLik(logistic_res))
      res_optim <- list(par,loglikelihood,hessian)
      names(res_optim) <- c("par","loglikelihood","hessian")
    } else
    {
      if (is.null(init)) {
        init_param <- logistic_loglogistic_init()
      } else
      {
        init_param <- init
      }
      res_optim <- optim(init_param, logistic_loglogistic_opt, method="L-BFGS-B", lower=c(rep(-Inf,m1+m2+2),0.01), upper=c(rep(Inf,m1+m2+2),100), control=list(trace=TRUE), hessian=TRUE)
      names(res_optim)[names(res_optim)=="value"] <- "loglikelihood"
      res_optim$loglikelihood <- -1*res_optim$loglikelihood
      res_optim$convergence<- list("convergence"=res_optim$convergence, "counts"=res_optim$counts, "message"=res_optim$message)

    }
  } else if (model=="logistic-lognormal")
  {
    ### Special cases ###
    ### No prevalent disease ###
    if (n1==0 & n3==0)
    {
      print ("Incident disease only")
      xsam2 <- xsam
      xsam2[which(xsam2[,3]==0),3] <- NA
      xsam2[which(xsam2[,4]==Inf),4] <- NA

      if (m2==0)
      {
        lnorm_res <- survreg(Surv(xsam2[,3],xsam2[,4],type='interval2')~1, dist="lognormal")
        par <- c(lnorm_res$coefficients[1], lnorm_res$scale)
        names(par) <- c("lognormal-meanlog","lognormal-sdlog")
      } else if (m2>0)
      {
        #subset <- which(xsam2[,2]==0)
        lnorm_res <- survreg(Surv(xsam2[,3],xsam2[,4],type='interval2')~ mat2, dist="lognormal")
        par <- c(lnorm_res$coefficients, lnorm_res$scale)

        if(is.null(colnames(mat2))) {
          names(par) <- c(paste0("lognormal-",names(lnorm_res$coefficients)),"lognormal-sdlog")
        } else {
          names(par) <- c(paste0("lognormal-",c("Intercept", colnames(mat2))),"lognormal-sdlog")
        }

      }
      loglikelihood <- lnorm_res$loglik[1]
      hessian <- solve(lnorm_res$var)
      res_optim <- list(par,loglikelihood,hessian)
      names(res_optim) <- c("par","loglikelihood","hessian")
    } else if (n2==0 & n3==0)   ### No incident disease ###
    {
      print ("Prevalent disease only")
      xsam2 <- xsam
      xsam2[which(is.na(xsam2[,2])),2] <- 1

      if (m1==0)
      {
        logistic_res <- glm(xsam2[,2] ~ 1, family = "binomial")
        par <- logistic_res$coefficients[1]
        names(par) <- "logistic-Intercept"
        hessian <- 1/summary(logistic_res)$coefficients[2]^2
      } else if (m1>0)
      {
        logistic_res <- glm(xsam2[,2] ~ mat1, family = "binomial")
        par <- logistic_res$coefficients

        if(is.null(colnames(mat1))) {
          names(par) <- paste0("logit-",names(logistic_res$coefficients))
        } else {
          names(par) <- paste0("logit-", c("(Intercept)", colnames(mat1)))
        }
        hessian <- solve(vcov(logistic_res))
      }
      loglikelihood <- as.numeric(logLik(logistic_res))
      res_optim <- list(par,loglikelihood,hessian)
      names(res_optim) <- c("par","loglikelihood","hessian")

    } else
    {
      if (is.null(init)) {
        init_param <- logistic_lognormal_init()
      } else
      {
        init_param <- init
      }
      res_optim <- optim(init_param, logistic_lognormal_opt, method="L-BFGS-B", lower=c(rep(-Inf,m1+m2+2),0.01), upper=c(rep(Inf,m1+m2+2),100), control=list(trace=TRUE), hessian=TRUE)
      names(res_optim)[names(res_optim)=="value"] <- "loglikelihood"
      res_optim$loglikelihood <- -1*res_optim$loglikelihood
      res_optim$convergence<- list("convergence"=res_optim$convergence, "counts"=res_optim$counts, "message"=res_optim$message)
    }
  } else if (model=="logistic-gengamma")
  {
    ### Special cases ###
    ### No prevalent disease ###
    if (n1==0 & n3==0)
    {
      print ("Incident disease only")
      xsam2 <- xsam
      xsam2[which(xsam2[,3]==0),3] <- NA
      xsam2[which(xsam2[,4]==Inf),4] <- NA

      if (m2==0)
      {
        ggamma_res <- survreg(Surv(xsam2[,3],xsam2[,4],type='interval2')~1, dist="weibull")
        param_init <- c(ggamma_res$coefficients[1], ggamma_res$scale,1)
        names(param_init) <- c("Generalized Gamma-mu","Generalized Gamma-sigma","Generalized Gamma-Q")
      } else if (m2>0)
      {
        ggamma_res <- survreg(Surv(xsam2[,3],xsam2[,4],type='interval2') ~ mat2, dist="weibull")
        param_init <- c(ggamma_res$coefficients, ggamma_res$scale,1)

        if(is.null(colnames(mat2))) {
          names(param_init) <- c(paste0("Generalized Gamma-",names(ggamma_res$coefficients)),
                                 "Generalized Gamma-sigma","Generalized Gamma-Q")
        } else {
          names(param_init) <- c(paste0("Generalized Gamma-",c("Intercept", colnames(mat2))),
                                 "Generalized Gamma-sigma","Generalized Gamma-Q")
        }

      }

      #full log-likelihood
      logistic_gengamma_likfun2 <- function(parest)
      {
        li <- xsam[,3]
        ri <- xsam[,4]

        if (m2==0)
        {
          s_lp <- rep(parest[1],usable)
          SL <- 1-pgengamma(li,s_lp,parest[2],parest[3])
          SR <- 1-pgengamma(ri,s_lp,parest[2],parest[3])
        } else if (m2>0)
        {
          s_lp <- parest[1]+mat2%*%parest[2:(m2+1)]
          SL <- 1-pgengamma(li,s_lp,parest[(m2+2)],parest[(m2+3)])
          SR <- 1-pgengamma(ri,s_lp,parest[(m2+2)],parest[(m2+3)])
        }

        #log-likelihood
        tmp <- rep(NA,usable)
        tmp <- log(SL-SR)
        sum(tmp) #full log-likelihood
      }

      #minimize function to maximize likelihood
      logistic_gengamma_opt2 <- function(y)
      {
        -1*logistic_gengamma_likfun2(y)
      }

      res_optim <- optim(init_param, logistic_gengamma_opt2, method="L-BFGS-B", lower=c(rep(-Inf,m2+1),0.01,-Inf), upper=c(rep(Inf,m2+1),100,Inf), control=list(trace=TRUE), hessian=TRUE)
      names(res_optim)[names(res_optim)=="value"] <- "loglikelihood"
      res_optim$loglikelihood <- -1*res_optim$loglikelihood
      res_optim$convergence<- list("convergence"=res_optim$convergence, "counts"=res_optim$counts, "message"=res_optim$message)

    } else if (n2==0 & n3==0)   ### No incident disease ###
    {
      print ("Prevalent disease only")
      xsam2 <- xsam
      xsam2[which(is.na(xsam2[,2])),2] <- 1

      if (m1==0)
      {
        logistic_res <- glm(xsam2[,2] ~ 1, family = "binomial")
        par <- logistic_res$coefficients[1]
        names(par) <- "logit-(Intercept)"
        hessian <- 1/summary(logistic_res)$coefficients[2]^2
      } else if (m1>0)
      {
        logistic_res <- glm(xsam2[,2] ~ mat1, family = "binomial")
        par <- logistic_res$coefficients

        if(is.null(colnames(mat1))) {
          names(par) <- paste0("logit-",names(logistic_res$coefficients))
        } else {
          names(par) <- paste0("logit-", c("(Intercept)", colnames(mat1)))
        }

        hessian <- solve(vcov(logistic_res))
      }
      loglikelihood <- as.numeric(logLik(logistic_res))
      res_optim <- list(par,loglikelihood,hessian)
      names(res_optim) <- c("par","loglikelihood","hessian")
    } else
    {
      if (is.null(init)) {
        init_param <- logistic_gengamma_init()
      } else
      {
        init_param <- init
      }
      res_optim <- optim(init_param, logistic_gengamma_opt, method="L-BFGS-B", lower=c(rep(-Inf,m1+m2+2),0.01,-Inf), upper=c(rep(Inf,m1+m2+2),100,Inf), control=list(trace=TRUE), hessian=TRUE)
      names(res_optim)[names(res_optim)=="value"] <- "loglikelihood"
      res_optim$loglikelihood <- -1*res_optim$loglikelihood
      res_optim$convergence<- list("convergence"=res_optim$convergence, "counts"=res_optim$counts, "message"=res_optim$message)
    }
  } else if (model=="logistic-gamma")
  {
    ### Special cases ###
    ### No prevalent disease ###
    if (n1==0 & n3==0)
    {
      print ("Incident disease only")
      xsam2 <- xsam
      xsam2[which(xsam2[,3]==0),3] <- NA
      xsam2[which(xsam2[,4]==Inf),4] <- NA

      if (m2==0)
      {
        gamma_res <- survreg(Surv(xsam2[,3],xsam2[,4],type='interval2')~1, dist="weibull")
        param_init <- c(gamma_res$coefficients[1], gamma_res$scale**(-2))
        names(param_init) <- c("Gamma-mu","Gamma-shape")
      } else if (m2>0)
      {
        gamma_res <- survreg(Surv(xsam2[,3],xsam2[,4],type='interval2')~mat2, dist="weibull")
        param_init <- c(gamma_res$coefficients, gamma_res$scale**(-2))

        if(is.null(colnames(mat2))) {
          names(param_init) <- c(paste0("Gamma-",names(gamma_res$coefficients)),   "Gamma-shape")
        } else {
          names(param_init) <- c(paste0("Gamma-",c("(Intercept)", colnames(mat2))),"Gamma-shape")
        }
      }

      #full log-likelihood
      logistic_gamma_likfun2 <- function(parest)
      {
        li <- xsam[,3]
        ri <- xsam[,4]

        if (m2==0)
        {
          s_lp <- rep(parest[1],usable)
          SL <- 1-pgamma(li,parest[2],exp(-1*s_lp)*parest[2])
          SR <- 1-pgamma(ri,parest[2],exp(-1*s_lp)*parest[2])
        } else if (m2>0)
        {
          s_lp <- parest[1]+mat2%*%parest[2:(m2+1)]
          SL <- 1-pgamma(li,parest[(m2+2)],exp(-1*s_lp)*parest[(m2+2)])
          SR <- 1-pgamma(ri,parest[(m2+2)],exp(-1*s_lp)*parest[(m2+2)])
        }

        #log-likelihood
        tmp <- rep(NA,usable)
        tmp <- log(SL-SR)
        sum(tmp) #full log-likelihood
      }

      #minimize function to maximize likelihood
      logistic_gamma_opt2 <- function(y)
      {
        -1*logistic_gamma_likfun2(y)
      }

      res_optim <- optim(init_param, logistic_gamma_opt, method="L-BFGS-B", lower=c(rep(-Inf,m2+1),0.01), upper=c(rep(Inf,m2+1),Inf), control=list(trace=TRUE), hessian=TRUE)
      names(res_optim)[names(res_optim)=="value"] <- "loglikelihood"
      res_optim$loglikelihood <- -1*res_optim$loglikelihood
      res_optim$convergence<- list("convergence"=res_optim$convergence, "counts"=res_optim$counts, "message"=res_optim$message)

    } else if (n2==0 & n3==0)   ### No incident disease ###
    {
      print ("Prevalent disease only")
      xsam2 <- xsam
      xsam2[which(is.na(xsam2[,2])),2] <- 1

      if (m1==0)
      {
        logistic_res <- glm(xsam2[,2] ~ 1, family = "binomial")
        par <- logistic_res$coefficients[1]
        names(par) <- "logit-(Intercept)"
        hessian <- 1/summary(logistic_res)$coefficients[2]^2
      } else if (m1>0)
      {
        logistic_res <- glm(xsam2[,2] ~ mat1, family = "binomial")
        par <- logistic_res$coefficients

        if(is.null(colnames(mat1))) {
          names(par) <- paste0("logit-",names(logistic_res$coefficients))
        } else {
          names(par) <- paste0("logit-", c("(Intercept)", colnames(mat1)))
        }
        hessian <- solve(vcov(logistic_res))
      }
      loglikelihood <- as.numeric(logLik(logistic_res))
      res_optim <- list(par,loglikelihood,hessian)
      names(res_optim) <- c("par","loglikelihood","hessian")
    } else
    {
      if (is.null(init))
      {
        init_param <- logistic_gamma_init()
      } else
      {
        init_param <- init
      }
      res_optim <- optim(init_param, logistic_gamma_opt, method="L-BFGS-B", lower=c(rep(-Inf,m1+m2+2),0.01), upper=c(rep(Inf,m1+m2+2),Inf), control=list(trace=TRUE), hessian=TRUE)
      names(res_optim)[names(res_optim)=="value"] <- "loglikelihood"
      res_optim$loglikelihood <- -1*res_optim$loglikelihood
      res_optim$convergence<- list("convergence"=res_optim$convergence, "counts"=res_optim$counts, "message"=res_optim$message)
    }
  }

  #res_optim$hessian is the negative hessian, we need to convert it to the real hessian
  res_optim$hessian<- -res_optim$hessian
  res_optim$covariance <- solve(-res_optim$hessian)
  par.se<- sqrt(diag(res_optim$covariance))
  CI.LL <- res_optim$par - qnorm(0.975)*par.se
  CI.UL <- res_optim$par + qnorm(0.975)*par.se

  res_optim$par<-  data.frame("Coef." = res_optim$par, "SE" = par.se, "95%LL"=CI.LL, "95%UL" = CI.UL, check.names=FALSE)
  row.names(res_optim$par) <- names(init_param)
  sep.names<-strsplit(row.names(res_optim$par), "-")

  ##modified to allow using character "-" in coding nominal variables #######

      name.mat<-plyr::ldply(sep.names, rbind)

      if(ncol(name.mat)>2){
        label2<-apply(name.mat[,-1],1,function(x) paste(x,collapse = "-"))
        label3<-sub("-NA","",label2)
        labels<-data.frame(Model=as.character(name.mat[,1]),Label=label3,stringsAsFactors=FALSE)
        remove(label2,label3)
      }else{
        labels<-data.frame(Model=as.character(name.mat[,1]),Label=as.character(name.mat[,2]),stringsAsFactors=FALSE)
      }
      remove(sep.names,name.mat)

  labels$Model[which(labels$Model=="logistic")]<- "logit"
  res_optim$regression.coef<- data.frame(labels, res_optim$par, row.names=NULL, check.names=FALSE,
                                         stringsAsFactors=FALSE)
  res_optim$model<- model

  cols<- c("regression.coef", "HR", "OR", "covariance",  "hessian", "convergence", "loglikelihood", "model")
  index<- match(cols, names(res_optim))
  index<- index[!is.na(index)]

  return(res_optim[index])

}#the end of function PIMixture.Est



#' Prevalence-Incidence Mixture Models
#'
#' This package fits Prevalence Incidence Mixture models to data for which the time to event is interval censored or
#' for which the event is prevalent at the time zero but is only partially observed.
#' Such data often arises in medical screening for asymptomatic disease or disease precursors, such as precancerous lesions.
#' In such data, 1) onset of incident disease occurs between screening visits (interval-censoring),
#' 2) the disease may have already occurred before the initial screen (prevalent at time zero)
#' but may be initially missed and found some time after the initial screen.
#' These models estimates absolute and relative risks.
#' Semi-parametric (iterative convex minorant algorithm by Robertson,Wright and Dykstra, 1988), weakly-parametric (integrated B-splines), and fully parametric members of the Prevalence-Incidence Mixture model family are supported.
#' A non-parametric estimator (Turnbull 1976) is provided and is useful for checking parametric assumptions of fully parametric Prevalence Incidence Mixture models.
#' Only weakly-parametric and semi-parametric models (no variance calculation) currently support stratified random samples in the two viewpoints,
#' a superpopulation and a finite population. The superpopulation is to view the finite population of interest as an independent sample of size N from an infinite superpopulation. A later version will add this functionality for the logistic-Weibull and logistic-exponential models.
#' Semi-parametric, weakly-parametric models, logistic-Weibull, and logistic-exponential uses a logistic regression model as the prevalence model
#' and a proportion hazard survival model as the incidence model.
#' The semi-parametric model makes no assumptions regarding the baseline hazard function.  However, it can be computationally expensive when the unique number of visit times are over hundreds.
#' The weakly-parametric model approximates the bazeline hazard function using integrated B-splines and is faster.
#' When parametric assumptions can be made, the fully parametric models are fastest.
#' The following parametric assumptions are supported: logistic-Weibull, logistic-exponential, logistic-lognormal,
#' logistic-loglogistic, logistic-gengamma, and logistic-gamma.  Variance estimates are available only for the weakly-parametric,
#' logistic-Weibull, logistic-exponential, and non-parametric.  For the non-parametric, this is achieved through boot-strapping by setting
#' the conf.int parameter to "TRUE" and can be very computationally expensive.
#' For identifiability of the mixture model, the data must contained observed prevalent disease and interval-censored incident disease.
#'
#'
#' @import optimx survival fdrtool interval Icens plyr
#'
#' @section Warning:
#'  The model="semi-parametric" option is very computationally expensive when the unique visit times are over hundreds.
#'
#' @param p.model The prevalence model to be fitted, specified using an expression of the form \emph{c+l+r~model}.
#'  Elements in the expression are as followed:
#'  \itemize{
#'  \item c - Numeric variable indicating whether the event was prevalent at time zero, taking values of 1=="Yes", 0="No", -999="Latent";
#'  \item l - Numeric starting time of the interval in which event occurred, with -999 denoting known prevalent events;
#'  \item r - Ending time of the interval in which event occurred, with -999 and Inf denoting known prevalent events and right-censoring, respectively;
#'  \item model - Linear predictor consisting of a series of terms separated by \emph{+} operators.
#'  }
#' @param i.model The incidence model to be fitted, specified using an expression of the form \emph{c+l+r~model} (see p.model). Defaults to p.model.
#' @param data Data used to fit the model containing columns for each term in p.model and i.model expressions.
#'  For stratified random sampling designs, columns denoted samp.weight and strata are expected indicating the sampling weights and sampling strata.
#'  For sample.design="superpopulation" option, an additional column denoted strata.frac is expected indicating the fraction of the population
#'  that consists of each strata.  For example, if in the target population there are three strata that occurs with proportions 0.2, 0.4, and 0.6,
#'  then strata.frac will take values of 0.2, 0.4 or 0.6.
#' @param model Character string indicating the specific member of the Prevalence-Incidence Mixture Model family to be fitted.  Options are:
#'  \itemize{
#'  \item "semi-parametric" Fits logistic regression and proportional hazards model as the prevalence and incidence models, respectively.
#'  The baseline hazard funcion is non-parametrically estimated using the iterative convex minorant algorithm.
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
#'  \item "non-parametric" Provides the non-parametric cumulative risk estimator.  This is akin to the non-parametric estimates provided using
#'  the Turnbull methods.  Covariates are not supported.  Confidence intervals are obtained through bootstrapping but is computationally expensive.
#'  }
#'  Variance estimates are not available for "semi-parametric", logistic-lognormal", "logistic-loglogistic", "logistic-gengamma", or "logistic-gamma",
#'  but they can be obtained using bootstrap methods.
#'  Defaults to "semi-parametric".
#' @param sample.design Sampling design of the NULL="simple random sampling", 1="finite population", 2="superpopulation".  Defaults to NULL.
#'        For "superpopulation", N is required for variance calculation (only provided when model="weakly-parametric" is used)
#' @param N Population size, required for superpopulation.  Defaults to NULL.
#' @param reg.initials Initial parameter estimates.  Defaults to NULL.
#' @param conf.int For non-parametric model option, FALSE="Do not obtain bootstrap confidence intervals", TRUE="Obtain bootstrap confidence intervals.
#'        Defaults to FALSE.
#' @param convergence.criteria Convergence of models occurs when reduction in the objective is within this convergence tolerance.  Defaults to 0.001.
#' @param iteration.limit Maximum number of iterations allowed to achieve convergence.  Defaults to 250.
#' @param time.interval Define time intervals to output baseline hazards for.  Defaults to .01.
#' @param n.knots Number of knots for splines for "weakly-parametric" model.  Defaults to 5.
#' @param order Degree of splines for "weakly-parametric" model.  Defaults to 4 (cubic splines).
#' @param max.time Define maximum time to output baseline hazards for.  Defaults to the largest finite start/end time given in the data.
#' @param design.out Option to include the design matrix of data used for model fitting in the output.  Defaults to TRUE.
#'
#' @return The output is a list of class PIMix which contains the following elements.
#'  \itemize{
#'  \item data.summary A data frame containing the following:
#'    Num. of observation - total number of observations in data set;
#'    Included subjects - number of observations used in fitting model;
#'    Known prevalent cases - the number of events known to be prevalent at time zero;
#'    Interval censoring - the number of event times occuring in the interval (L>0,R<Inf];
#'    Left censoring - the number of event times known to occur by R<Inf, but can also have been prevalent at time zero;
#'    Right censoring - the number of observations right-censored with event time occurring in the interval (L>0,Inf);
#'    Missing prevalent status - the number of observations where it is unknown whether the event was prevalent at time zero;
#'    Non-informative intervals - the number of observations with intervals (0,Inf) or [0,Inf) (denoting missing prevalent status).
#'  \item regression.coef A data frame summarizing parameter values, standard errors, and 95 percent confidence intervals.
#'  \item OR A data frame summary odds ratios, , standard errors, and 95 percent confidence intervals.
#'  \item HR A data frame summary hazard ratios, , standard errors, and 95 percent confidence intervals.
#'  \item knots If model="weakly-parametric" is specified, this is a numeric vector of starting time points for each exponential spline.
#'  \item exp.spline.coeff If model="weakly-parametric" is specified, this is a numeric vector of coefficients for each exponential spline.
#'  \item cum.hazard If model="semi-parametric" or model="weakly-parametric" is specified, this is a data frame containing the baseline cumulative hazard.
#'  \item covariance A matrix containing the covariance matrix for the parameters (not produced for model="semi-parametric").
#'  \item hessian A matrix containing the hessian matrix for the parameters (not produced for model="semi-parametric").
#'  \item model Character string indicating the specific member of the Prevalence-Incidence Mixture Model family fitted.
#'  \item p.model The prevalence model.
#'  \item prev.design The design matrix for the prevalence model.
#'  \item i.model The incidence model.
#'  \item incid.design The design matrix for the incidence model.
#'  \item loglikelihood For random samples, this is the log-likelihood of the fitted model.
#'        For stratified random samples, the weighted-likelihood approach is used and a log-pseudolikelihood (weighted log-likelihood) is reported.
#'  \item convergence Convergence statistics.
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
#'  \item Turnbull BW (1976). The empirical distribution with arbitrary grouped censored and truncated
#'          data. Journal of the Royal Statistical Society - Series B (Statistical Methodology) 38, 290-295.
#'  \item Robertson T, Wright FT, and Dykstra RL (1988). Order Restricted Statistical Inference. Wiley.
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
#'
#'
#'  #For stratified random samples
#'  model3<-"C+L+R~X1+X2"

#'  #sample.design=1 indicates the target population is a finite population, and the variance is design-based.
#'  fit5<-PIMixture(p.model=model3,data=PIdata2, model="weakly-parametric",n.knots=7,order=4,sample.design=1)
#'
#'  #sample.design=2 indicates the target population is a superpopulation, and the variance consists of
#'  #design-based and model-based variances. Generally, the Variance in the superpopulation frame is slightly larger than the variance in the finite population frame.
#'  fit6<-PIMixture(p.model=model3,data=PIdata2, model="weakly-parametric",n.knots=7,order=4,sample.design=2,N=10000)
#'
#'  fit7<-PIMixture(p.model=model3,data=PIdata2, model="semi-parametric",sample.design=1)
#'
#'
PIMixture<-function(p.model,i.model,data,model="semi-parametric",reg.initials=NULL,conf.int=FALSE,convergence.criteria=0.001,
                    iteration.limit=250,time.interval=1e-2,design.out=TRUE,sample.design=NULL,N=NULL,n.knots=5,order=4,
                    max.time,...){
  requireNamespace("survival", quietly = TRUE)
  requireNamespace("optimx", quietly = TRUE)
  requireNamespace("fdrtool", quietly = TRUE)
  requireNamespace("interval",quietly = TRUE)
  requireNamespace("flexsurv", quietly = TRUE)
  requireNamespace("plyr",quietly = TRUE)

  if (missing(p.model)){
    stop("Regression model is missing")
  }

  if (missing(data)) {
    stop("data is not provided")
  }

  if (!model%in%c("semi-parametric",   "weakly-parametric", "logistic-Weibull", "logistic-exponential",
                  "logistic-lognormal","logistic-loglogistic", "logistic-gengamma","logistic-gamma", "non-parametric")) {
    stop("the model argument is missing, model can be chosen from \"semi-parametric\",   \"weakly-parametric\", \"logistic-Weibull\", \"logistic-exponential\",\"logistic-lognormal\",\"logistic-loglogistic\",\"logistic-gengamma\",\"logistic-gamma\", \"non-parametric\", ")
  }

  if(!is.null(sample.design)){
    if(missing(N)&sample.design==2&model=="weakly-parametric"){
      warning("Variance estimation for superpopulation is not available because population size N is missing")
    }
  }

  #change variable names for outcomes, C, L, R
  fml<-as.formula(p.model)

  if(attr(terms(fml),"response")!=1){
    stop("No response variable in p.model")
  }
  allvars<-all.vars(fml)

  ######################################################
  # save p.model in the output, since the response variable in
  # p.model will be modified in the next step
  ######################################################
  output<-list()
  output$p.model<- p.model

  ######################################################
  # Formulate p.model
  # change the response variable in p.model into C+L+R
  ######################################################
  bb<-strsplit(p.model,split="~")
  p.model<-paste0("C+L+R~",bb[[1]][2])
  fml<-as.formula(p.model)

  remove(bb)

  ######################################################
  # Formulate i.model
  # Assume the response variable name in i.model
  # is the same as the response variable in p.model
  #
  ######################################################

  if(!missing(i.model)){
    bb<-strsplit(i.model,split="~")
    i.model<-paste0("~",bb[[1]][2])

    ####################################################
    # 1) adding the covariates in i.model into allvars
    # 2) ignore the responses in i.model, assuming the
    # response in i.model is the same as that of p.model
    #####################################################
    fml2<-as.formula(i.model)
    allvars2<-all.vars(fml2)
    allvars<-union(allvars,allvars2)

    output$i.model<- i.model

    #######################################################
    #Adding back the response variable "C+L+R" into i.model
    #######################################################

    i.model<- paste0("C+L+R", i.model)
    fml2<-as.formula(i.model)

    remove(bb)
  } else {
    fml2<- NULL
    bb<-strsplit(p.model,split="~")
    output$i.model<- paste0("~",bb[[1]][2])
  }

  #####################################
  # sampling weight options
  #####################################

  if(!is.null(sample.design)) {
    if(!sample.design %in%c(1, 2)) {
      warning("\"sample.design\" is out of range, it should be either NULL, 1 or 2. So NULL sample.design will be used")
      sample.design<- NULL
    }
  }

  if(c("samp.weight") %in% names(data)) {
    allvars<-c(allvars,"samp.weight")
    if(sum(data$samp.weight<0)>0){
      stop("\"samp.weight\" should be positive.")
    }
  } else if(!(c("samp.weight") %in% names(data))){
    allvars<-c(allvars,"samp.weight")
    data$samp.weight<-1
    if(!is.null(sample.design)) {
      warning("\"samp.weight\" is not included in the data, so we run no weighted model.")
    }
  }

  #########################################################
  # strata and phase 2 sampling, mainly for weakly-parametric
  ##########################################################

  if(c("strata") %in% names(data)) {
    allvars<-c(allvars,"strata")
  } else if(!(c("strata") %in% names(data))){
    allvars<-c(allvars,"strata")
    data$strata<-1
    if(!is.null(sample.design)){
      warning("\"strata\" is not included in the data, so we consider no stratified sample.")
    }
  }

  if(!is.null(sample.design)) {
    if((sample.design==2)&(c("strata.frac") %in% names(data))){
      allvars<-c(allvars,"strata.frac")
    } else if((sample.design==2)&!(c("strata.frac") %in% names(data))){
      warning("\"strata.frac\" are needed to cacluate variance estimation for superpopulation.")
    }
  }

  ########################################################
  # set samp.data
  #######################################################
  if(sum(allvars%in%names(data)) != length(allvars)) {
    stop(paste0("The input dataset does not contain variables ", allvars[!allvars%in%names(data)]))
  }

  samp.data2<-data[,allvars]
  names(samp.data2)[names(samp.data2)==allvars[1]]<-"C"
  names(samp.data2)[names(samp.data2)==allvars[2]]<-"L"
  names(samp.data2)[names(samp.data2)==allvars[3]]<-"R"

  samp.data2[is.na(samp.data2$C), "C"]<- -999   #ci with "NA" values are treated as "unknown"
  samp.data2[samp.data2$C == 1,   "L"]<- -999   #For prevalent cases L & R = -999
  samp.data2[samp.data2$C == 1,   "R"]<- -999

  ########################################################
  # Remove all the missing covariates
  ########################################################
  samp.data<-samp.data2[complete.cases(samp.data2),]  #excluding cases with NA
  check.data(samp.data, allvars)

  ############# VALIDATE SAMPLE DATA ###########################################################

  ############# GET THE DATA SUMMARY ###########################################################
  n.obs <-nrow(samp.data2)       #number of observations
  temp.summary<- get.summary(samp.data)
  output$data.summary<- data.frame(label= c("Num. of observations", temp.summary[,1]),
                                   no.obs=c(n.obs, temp.summary[,2]), check.names=FALSE)
  remove(n.obs)

  ###################### design matrix ##############################################################
  mf <- model.frame(formula=fml, data=samp.data)
  design.mat <- model.matrix(attr(mf, "terms"), data=mf)
  n.beta<-ncol(design.mat)

  length.mf<-length(colnames(mf))
  if (length.mf>1){
    for(j in 2:length.mf){
      if(is.factor(mf[,j])){
        level.list<-levels(mf[,j])

        for (k in 1:length(level.list)) {
        colnames(design.mat)<-gsub(level.list[k],paste0("=",level.list[k]),colnames(design.mat),fixed=T)
        }
        remove(level.list)
      }
    } #for(j in 2:length.mf)
  }#if (length.mf>1)
  remove(length.mf)

  if(missing(i.model)){
    if(n.beta>1){
      xmat<-matrix(design.mat[,-1],  nrow=nrow(design.mat))
      colnames(xmat)<-colnames(design.mat)[-1]
    }else if(n.beta==1){
      xmat<-0
    }
    n.gamma<-n.beta-1
  } else if(!missing(i.model)){
    mf2 <- model.frame(formula=fml2, data=samp.data)
    xmat <- model.matrix(terms(fml2), data=mf2)
    length.mf2<-length(colnames(mf2))
    if (length.mf2>1){
        for(j in 2:length.mf2){
          if(is.factor(mf2[,j])){
            level.list<-levels(mf2[,j])

            for (k in 1:length(level.list)) {
              colnames(xmat)<-gsub(level.list[k],paste0("=",level.list[k]),colnames(xmat),fixed=T)
            }

            remove(level.list)
          }
        } #for(j in 2:length.mf)
    }#if(length.mf2>1)
    remove(length.mf2)

    if(ncol(xmat)>1){
      temp<- colnames(xmat)
      xmat<- matrix(xmat[,-1], nrow=nrow(xmat))
      colnames(xmat)<- temp[-1]
      n.gamma<-ncol(xmat)
    } else if(ncol(xmat)==1){
      xmat<-0
      n.gamma<-0
    }
  }

  samp.data$K<-as.numeric(samp.data$C!=-999)
  #n.regpara<-n.beta+n.gamma


  #########################################################################
  # The above is to validate all the parameters. This part is to do
  # the real calculation
  #########################################################################

  output2<- list()
  if (model== "semi-parametric") {
    output2<- ipw.lc.semipara(samp.data,data,fml,fml2,n.beta,n.gamma, design.mat, xmat, reg.initials,
                              convergence.criteria,iteration.limit,time.interval)

  } else if (model == "weakly-parametric") {
    output2<- ipw.lc.splines(samp.data, n.beta, n.gamma, design.mat, xmat,
                             sample.design,N,n.knots,order,max.time,reg.initials,
                             convergence.criteria,iteration.limit,time.interval)

  } else if(model%in% c("logistic-Weibull", "logistic-exponential",
                        "logistic-lognormal","logistic-loglogistic","logistic-gengamma","logistic-gamma", "non-parametric")) {

    ci<- samp.data$C
    li<- samp.data$L
    ri<- samp.data$R

    ci[which(ci==-999)]<- NA
    li[which(ci==1)] <- 0
    ri[which(ci==1)] <- 0

    if(ncol(design.mat) ==1) {
      mat1<- NULL
    } else {
      mat1<- matrix(design.mat[, -1], nrow=nrow(design.mat))
      colnames(mat1)<- colnames(design.mat)[-1]
    }

    mat2<- xmat
    if(class(mat2) == "numeric") {
      mat2<- NULL
      xmat<- cbind(design.mat[,"(Intercept)"])
      colnames(xmat)<- "(Intercept)"
    } else {
      xmat<- cbind(1, xmat)
      colnames(xmat) <- c("(Intercept)", colnames(mat2))
    }

    output2<- PIMixtureEst(ci=ci,li=li, ri=ri, mat1=mat1, mat2=mat2, model=model, conf.int=conf.int, init=reg.initials)

  } else {
    stop("wrong \"model\" option")
  }

  output<- c(output, output2)
  if(design.out) {
    output$prev.design<-design.mat
    output$incid.design<- xmat
  }

  output$model<- model
  class(output)<- "PIMix"

  ##################################################################
  # The output of non-parametric model is different from those
  # of parametric and semi-parametric/weakly parametric model
  ##################################################################
  if (model == "non-parametric") {
    return(output)
  } else {
    cols<- c("data.summary", "regression.coef", "exp.spline.coeff", "HR", "OR", "cum.hazard",
             "cumrisk.est", "cumrisk.est.full","covariance", "hessian",  "convergence", "loglikelihood",
             "knots", "model", "p.model", "i.model")

    if(design.out) {
      cols<- c(cols, "prev.design", "incid.design")
    }

    index<- match(cols, names(output))
    index<- index[!is.na(index)]

    if(design.out) {
      cols<- c(cols, "prev.design", "incid.design")
    }

    output<- output[index]
    class(output)<- "PIMix"
    return(output)
  }
}

