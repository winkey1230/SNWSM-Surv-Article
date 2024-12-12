################################################################################
######### build functions for Shape-constrained Nonlinear Weighted-Sum Model ###
######## standardized method: (x-x_min)/(x_max-x_min)                 ##########
######## date 2024/07/23 #######################################################
#################################################################################
library(fda)
library(survival)
library(doSNOW)
require(splines)
library(boot)


#' @description This function implemented the shape-constrained semi-parametric model in the 
#' framework of Cox proportional hazard regression. 
#' @param time A vector indicating the survival time
#' @param censor If censor, the value is 0, otherwise , 1
#' @param x A matrix or data.frame indicating the exposure variables
#' @param z A vector or matrix indicating the confounding variables
#' @param strata A vector indicating the strata variable
#' @param nknots A integer
#' @param beta_ini The initial values for beta and gamma
#' @param loss_inc_stop A logical value indicating whether the iteration process stops 
#' when the loss value emerge an increase for the first time.
snwsm_surv <- function(time,censor,x,z = NULL,strata = NULL,nknots = NULL,max.iter = 100,knots = NULL,loss_inc_stop = T,
                      degree = 3,eps = 1e-4,trace = FALSE,input_save = F,weightfun = function(x) x/sum(x),boot_for_beta = T){
  if(input_save){
    call_params <- formals(sys.function(sys.parent()))
    input_argue <- sapply(names(call_params),get,envir = environment())
  } else input_argue <- NULL
  if(!is.null(knots) & !is.null(nknots)) stop("knots and nknots exists simultaneously")
  if(anyNA(time)|anyNA(censor)|anyNA(x)) stop("Please remove the NA values before building the model")
  if(is.null(strata)) strata <- rep(1,length(time))
  if(!is.null(z)){
    if(anyNA(z)) stop("Please remove the NA values in Z before building the model")
    if("data.frame" %in% class(z)) z <- as.matrix(model.matrix(~.,z)[,-1])
    else z <- as.matrix(model.matrix(~ z)[,-1])  
    pz <- dim(z)[2] 
  } else pz <- 0 
  x <- as.matrix(x)
  if(dim(x)[2] == 1) stop("x must be a matrix with the number of columns larger than 1")
  origin_x <- x
  range_x <- apply(x, 2, function(a) max(a)-min(a))
  min_x <- apply(x, 2, min)
  
  x <- apply(x, 2, function(a) (a-min(a))/(max(a)-min(a)))
  n <- dim(x)[1]
  p <- dim(x)[2]
  nevent <- sum(1-censor)
  norder <- degree + 1 # the norder for bsplineS function
  
  xz <- cbind(x,z)
  fit_w <- coxph(Surv(time,1-censor)~xz + strata(strata))
  w_summaryNA <- summary(fit_w)$coefficients[1:p,-2]
  w_summaryNA[,] <- NA
  w_summaryNA[,1] <- 0
  
  ab <- summary(fit_w)$coefficients[1:p,c(1,5)]
  a <- ab[,1]
  b <- ab[,2]
  while(any(a<0)){
    e <- match(min(b[a<0]),b)
    a[e] <- 0
    xz <- cbind(x[,which(a != 0)],z)
    fit_w <- suppressWarnings(coxph(Surv(time,1-censor)~xz + strata(strata)))
    b[a != 0] <- summary(fit_w)$coefficients[1:sum(a != 0),5]
    a[a != 0] <- summary(fit_w)$coefficients[1:sum(a != 0),1]
  }
  w_ini <- a 
  if(!all(w_ini == 0)) w_ini <- weightfun(w_ini)
  lossvalue_ini <- -fit_w$loglik[2]/nevent
  if(trace) cat(0,": ",lossvalue_ini,w_ini,"\n")
  
  fit_beta <- NULL
 
  if(is.null(knots)) knots <- seq(0, 1, length = nknots + norder - 1)
  else {
    nknots <- length(knots)
    knots <- c(0,knots,1)
  }
  len.B <- length(knots) + norder - 3
  len.w <- length(w_ini)
  
  SIGMA <- matrix(-1,len.B,len.B)
  for (i in 1:len.B) {
    for (j in 1:len.B){
      if(i < j) SIGMA[i,j] <- 0 else SIGMA[i,j] <- 1
    }
  }
  control.optim <- list()
  
  w <- w_ini
  beta.B <- rep(0,len.B)
  beta.z <- if(all(w==0)) coef(fit_w) else coef(fit_w)[-(1:sum(w!=0))]
  s0 <- matrix(0,ncol = len.B,nrow = len.B)
  lossvalue <- lossvalue_ini
  w_summary <- w_summaryNA
  
  if(all(w==0)){
    warning("All the w is zero")
    conv <- T
  } else{
    conv <- FALSE
    iter <- 0
    while (iter < max.iter & conv == FALSE) {
      iter <- iter + 1
      
      w_summary <- w_summaryNA;
      w_summary[w != 0,] <- summary(fit_w)$coefficients[1:sum(w != 0),-2]
      w_summary[,1] <- weightfun(w_summary[,1])
      w_summary[,2] <- w_summary[,1]/w_summary[,3]
      
      u <- pu(x, w)
      eta <- u$u
      B <- bsplineS(eta, breaks = quantile(eta, knots),norder = norder)[,-1] # intercept is false
      TransB <- B%*%SIGMA
      
      ## 限制单调递增
      BZ <- cbind(TransB,z)
      fit_beta <- suppressWarnings(coxph(Surv(time,1-censor)~BZ + strata(strata)))
      ab <- summary(fit_beta)$coefficients[1:len.B,c(1,5)]
      a <- ab[,1]
      b <- ab[,2]
      while(any(a<0)){
        e <- match(min(b[a<0]),b)
        a[e] <- 0
        BZ <- cbind(TransB[,which(a != 0)],z)
        fit_beta <- suppressWarnings(coxph(Surv(time,1-censor)~BZ + strata(strata)))
        b[a != 0] <- summary(fit_beta)$coefficients[1:sum(a != 0),5]
        a[a != 0] <- summary(fit_beta)$coefficients[1:sum(a != 0),1]
      }
      lossvalue_new <- -fit_beta$loglik[2]/nevent
      loss_dec_continue <- ifelse(loss_inc_stop,lossvalue > lossvalue_new,T)  
      beta.B_new <- drop(SIGMA%*%a)
      beta.z_new <- coef(fit_beta)[-(1:sum(a!=0))]
      s0_new <- matrix(0,ncol = len.B,nrow = len.B) # s0为SIGMA转换前的样条基函数估计的系数
      s0_new[which(a!=0),which(a!=0)] <- vcov(fit_beta)[1:sum(a!=0),1:sum(a!=0)]
      if(all(beta.B_new == 0)) break
      if(loss_dec_continue){
        if(trace) cat(iter,": ",lossvalue_new,w,"\n")
        lossvalue <- lossvalue_new
        beta.B <- beta.B_new
        beta.z <- beta.z_new
        s0 <- s0_new
        final_w_summary <- w_summary
        final_w <- w
      } else break
      
      B_deriv <- bsplineS(eta, breaks = quantile(eta, knots),norder = norder, nderiv = 1)[,-1]
      newx <- drop((B_deriv * (u$deriv)) %*% beta.B) * x
      off <- B %*% beta.B - newx %*% w
      xz <- cbind(newx,z)
      fit_w <- suppressWarnings(coxph(Surv(time,1-censor)~xz + strata(strata)+offset(off)))
      ab <- summary(fit_w)$coefficients[1:p,c(1,5)]
      a <- ab[,1]
      b <- ab[,2]
      while(any(a<0)){
        e <- match(min(b[a<0]),b)
        a[e] <- 0
        xz <- cbind(newx[,which(a != 0)],z)
        fit_w <- (coxph(Surv(time,1-censor)~xz + strata(strata) + offset(off)))
        b[a != 0] <- summary(fit_w)$coefficients[1:sum(a != 0),5]
        a[a != 0] <- summary(fit_w)$coefficients[1:sum(a != 0),1]
      }
      w_new <- a
      if(sum(w_new) != 0) w_new <- weightfun(w_new) else break # 所有w都为0时，停止计算
      conv <- (max(abs(w - w_new)) < eps)
      w <- w_new
    }
    w <- final_w
    w_summary <- final_w_summary
  }
  
  beta.B.t <- solve(SIGMA) %*% beta.B # 重参数化后的估计值
  if(any(beta.B.t == 0) & boot_for_beta & !all(w == 0)){
    u <- pu(x, w)
    eta <- u$u
    B <- bsplineS(eta, breaks = quantile(eta, knots),norder = norder)[,-1] # intercept is false
    TransB <- B%*%SIGMA
    BZ <- cbind(TransB,z)
    data_boot <- cbind(time, censor,strata,BZ)
    betaBfun <- function(data,i){
      time <- data[i,1]
      censor <- data[i,2]
      strata <- data[i,3]
      XB <- data[i,4:(3+len.B)]
      if(is.null(z))
        fit_beta <- suppressWarnings(coxph(Surv(time,1-censor)~XB + strata(strata)))
      else{
        offsetZ <- drop(data[i,-c(1:(3+len.B))] %*% beta.z) 
        fit_beta <- suppressWarnings(coxph(Surv(time,1-censor)~XB + strata(strata) + offset(offsetZ)))
      }
      ab <- summary(fit_beta)$coefficients[1:len.B,c(1,5)]
      a <- ab[,1]
      b <- ab[,2]
      while(any(a<0)){
        e <- match(min(b[a<0]),b)
        a[e] <- 0
        if(all(a == 0)) break
        XB <- data[i,4:(3+len.B)][,which(a != 0)]
        if(is.null(z))
          fit_beta <- suppressWarnings(coxph(Surv(time,1-censor)~XB + strata(strata)))
        else{
          fit_beta <- suppressWarnings(coxph(Surv(time,1-censor)~XB + strata(strata) + offset(offsetZ)))
        }
        b[a != 0] <- summary(fit_beta)$coefficients[1:sum(a != 0),5]
        a[a != 0] <- summary(fit_beta)$coefficients[1:sum(a != 0),1]
      }
      a
    }
    set.seed(2)
    fitboot <- boot(data_boot,betaBfun,R = 99,stype = "i",parallel = "snow",strata = data_boot[,3])
    s0 <- cov(fitboot$t) 
  }
  
  u <- pu(x, w)
  eta <- u$u
  B <- bsplineS(eta, breaks = quantile(eta, knots),norder = norder)[,-1] # intercept is false
  fit.score <- drop(x%*%w) 
  fit.hazard.log <- drop(B%*%beta.B)
 
  beta.x <- w/range_x
  df <- sum(w != 0) + sum(beta.B != 0) + pz - 1
  AIC <- nevent*lossvalue*2 + 2*df
  BIC <- nevent*lossvalue*2 + log(nevent)*df
  attr(AIC,"df") <- attr(BIC,"df") <- df
  range.x <- range_x
  names(w) <- names(beta.x) <- names(range.x) <- colnames(origin_x)
  names(beta.z) <- colnames(z)
  list(w = w,beta.z = beta.z,beta.x = beta.x,beta.B = beta.B,fit.score = fit.score,
       fit.hazard.log = fit.hazard.log,w_summary = w_summary,
       range.x = range.x,converge = conv, AIC = AIC, BIC = BIC,knots = knots,x.normalized = x,
       input_argue = input_argue,min.x = min_x,df = df,s0 = s0,SIGMA = SIGMA)
} 

# tool function
pu <- function (x, beta) { 
  x <- as.matrix(x)
  d <- (length(beta) + 1)/2
  v <- x %*% beta
  alpha <- quantile(sqrt(rowSums(x * x)), 0.95)
  t <- (v + alpha)/(2 * alpha)
  u <- pbeta(t, d, d)
  deriv <- dbeta(t, d, d)/(2 * alpha)
  return(list(u = u, deriv = drop(deriv)))
}
score2u <- function(x,beta,score){
  x <- as.matrix(x)
  d <- (length(beta) + 1)/2
  v <- score
  alpha <- quantile(sqrt(rowSums(x * x)), 0.95)
  t <- (v + alpha)/(2 * alpha)
  u <- pbeta(t, d, d)
  deriv <- dbeta(t, d, d)/(2 * alpha)
  names(u) <- names(deriv) <- names(v)
  return(list(u = u, deriv = drop(deriv)))
}

#' @description predictive log(hazard) for a given score or newx
#' @param object an object from snwsm_surv
#' @param newx a vector or a matrix indicating the exposure values to be predicted
#' @param score a vector or a scalar indicating the score to be predicted
#' @param cenx a vector indicating the centred exposure values
predict.hazard <- function(object,newx = NULL,score = NULL,cenx = NULL,CI = 0.95){
  if((!is.null(newx) & !is.null(score)))
    stop("newx and score cannot be nonnull simutaneously")
  if(is.null(newx) & is.null(score)) score <- sort(unique(object$fit.score)) 
  eta <- pu(object$x.normalized,object$w)$u
  degree <- object$input_argue$degree
  if(is.null(cenx)){
    cenx <- apply(object$input_argue$x,2,median)
    cat("The cenx has been set as the median of x")
  }  
  cenx <- t(as.matrix(cenx))
  cenx.norm <- (t(cenx) - object$min.x)/object$range.x
  cenx.norm <- t(as.matrix(cenx.norm))
  cen.score <- drop(cenx.norm%*%object$w) 
  cen.u <- score2u(object$x.normalized,object$w,cen.score)$u
  cen.B <- splines::bs(cen.u,knots = quantile(eta, object$knots[-c(1,length(object$knots))]),
                       Boundary.knots = range(eta),degree = degree,intercept = F,warn.outside = F)
  cen.hazard.log <- drop(cen.B%*%object$beta.B) 
  s <- object$s0
  zalpha <- qnorm(CI/2+0.5)
  range.score <- range(fit$fit.score)
  # browser()
  if(!is.null(newx)){
    if(is.vector(newx)) newx <- t(as.matrix(newx))
    else newx <- as.matrix(newx)
    newx.norm <- (t(newx) - object$min.x)/object$range.x
    newx.norm <- t(as.matrix(newx.norm))
    pred.score <- drop(newx.norm%*%object$w)
    index.min <- pred.score < range.score[1]
    index.max <- pred.score > range.score[2]
    if(any(index.min)){
      warning(paste("This is",sum(index.min),"pred.score less than the minimal score, which has been set as the minimal score") )
      pred.score[pred.score < range.score[1]] <- range.score[1]
    } 
    if(any(index.max)){
      warning(paste("This is",sum(index.max),"pred.score larger than the maximal score, which has been set as the maximal score") )
      pred.score[pred.score > range.score[2]] <- range.score[2]
    } 
    u <- score2u(object$x.normalized,object$w,pred.score)$u
    # B <- bsplineS(u, breaks = quantile(eta, object$knots),norder = 3)[,-1] # intercept is false
    B <- splines::bs(u,knots = quantile(eta, object$knots[-c(1,length(object$knots))]),
                     Boundary.knots = range(eta),degree = degree,intercept = F,warn.outside = F)
    pred.HR <- exp(drop(B%*%object$beta.B) - cen.hazard.log)
    difB <- t(t(B) - c(cen.B))
    sd <- apply(difB%*%object$SIGMA, 1, function(x) sqrt(x%*%s%*%x))
    pred.HR_lower <- exp(log(pred.HR)-zalpha*sd)
    pred.HR_upper <- exp(log(pred.HR)+zalpha*sd)
    list(pred.score = pred.score,pred.HR = pred.HR,pred.HR_lower = pred.HR_lower,pred.HR_upper = pred.HR_upper,
         cen.hazard.log = cen.hazard.log,cen.score = cen.score)
  } else {
    u <- score2u(object$x.normalized,object$w,score)$u
    # B <- bsplineS(u, breaks = quantile(eta, object$knots),norder = 3)[,-1] # intercept is false
    B <- splines::bs(u,knots = quantile(eta, object$knots[-c(1,length(object$knots))]),
                     Boundary.knots = range(eta),degree = degree,intercept = F,warn.outside = F)
    pred.HR <- exp(drop(B%*%object$beta.B) - cen.hazard.log) 
    difB <- t(t(B) - c(cen.B))
    sd <- apply(difB%*%object$SIGMA, 1, function(x) sqrt(x%*%s%*%x))
    pred.HR_lower <- exp(log(pred.HR)-zalpha*sd)
    pred.HR_upper <- exp(log(pred.HR)+zalpha*sd)
    
    list(pred.score = score,pred.HR = pred.HR,pred.HR_lower = pred.HR_lower,pred.HR_upper = pred.HR_upper,
         cen.hazard.log = cen.hazard.log,cen.score = cen.score)
  }
}


#' @description using bootstrap to calculate 95% confidence interval
#' @param object an object from snwsm_surv
snwsm_bootstrap <- function(object,R = 100,ncore = 4,seed = 1,max.iter = 10){
  if(is.null(object$input_argue)) stop("The input arguement parameters for snwsm_surv cannot be Null")
  if(is.null(max.iter)) max.iter <- object$input_argue$max.iter
  n <- length(object$input_argue$time)
  set.seed(seed = 1)
  indices.boot <- matrix(sample(1:n,size = n*R,replace = T),ncol = R) 
  
  cl <- makeCluster(ncore)
  registerDoSNOW(cl)
  
  clusterExport(cl = cl,c("snwsm_surv","pu"))
  clusterExport(cl = cl,"bsplineS",envir = asNamespace("fda"))
  clusterExport(cl = cl,c("strata","coxph","Surv","coxph.control"),envir = asNamespace("survival"))
  
  pb <- txtProgressBar(max = R, style = 3)
  progress <- function(nn) setTxtProgressBar(pb, nn)
  opts <- list(progress = progress)
  
  fit.boot <- foreach(i = 1:R,.combine = "list",.multicombine = T,
                      .options.snow = opts,.maxcombine = R) %dopar%{
                        indices <- indices.boot[,i]
                        time <- object$input_argue$time[indices]
                        x <- object$input_argue$x[indices,]
                        z <- if(is.null(object$input_argue$z)) NULL 
                        else if(is.vector(object$input_argue$z)) object$input_argue$z[indices]
                        else object$input_argue$z[indices,]
                        strata <- if(is.null(object$input_argue$strata)) NULL 
                        else object$input_argue$strata[indices]
                        censor <- object$input_argue$censor[indices]
                        fit.R <- snwsm_surv(time = time,censor = censor,strata = strata,x = x,z = z,degree = object$input_argue$degree,
                                             nknots = object$input_argue$nknots,loss_inc_stop = object$input_argue$loss_inc_stop,
                                             max.iter = max.iter,eps = object$input_argue$eps,knots = object$input_argue$knots,
                                             trace = F,input_save = F,boot_for_beta = F)
                        fit.R$fit.hazard.log <- fit.R$fit.score <- NULL
                        fit.R
                      }
  close(pb)
  stopCluster(cl)
  fit.boot
}

snwsm_CI <- function(object,object.boot,newx = NULL,cenx = NULL,prob.ci = 0.95,length = 100,cenfun = median){
  if(is.null(newx)) {
    rangex <- apply(object$input_argue$x, 2, range)
    newx <- apply(rangex, 2, function(x) seq(x[1],x[2],length = length))
  } 
  if(is.vector(newx)) newx <- t(matrix(newx))
  if(is.null(cenx)) cenx <- apply(newx,2,cenfun)
  prob.par <- c(0.5-prob.ci/2,0.5+prob.ci/2)
  newx <- unique(newx)
  if(is.null(row.names(newx))) row.names(newx) <- 1:dim(newx)[1]
  
  w.boot <- t(sapply(object.boot, function(x) x$w))
  w.ci <- t(apply(w.boot, 2, quantile,prob = c(prob.par,0.5)))
  w.point <- object$w
  w.mean <- apply(w.boot, 2, mean)
  w.pvalue <- apply(w.boot, 2,function(x) mean(x <= 0))
  w <- cbind(w.point,w.mean,w.ci,w.pvalue)
  colnames(w) <- c("w.point","w.mean",paste0("w.",colnames(w.ci)),"w.pvalue")
  
  predict.hazard.boot <- function(objectx,newx,cenx,degree){
    object <- objectx
    eta <- pu(object$x.normalized,object$w)$u
    cenx <- t(as.matrix(cenx))
    cenx.norm <- (t(cenx) - object$min.x)/object$range.x
    cenx.norm <- t(as.matrix(cenx.norm))
    cen.score <- drop(cenx.norm%*%object$w) 
    cen.u <- score2u(object$x.normalized,object$w,cen.score)$u
    cen.B <- splines::bs(cen.u,knots = quantile(eta, object$knots[-c(1,length(object$knots))]),
                         Boundary.knots = range(eta),degree = degree,intercept = F,warn.outside = F)
    cen.hazard.log <- drop(cen.B%*%object$beta.B) 
    if(is.vector(newx)) newx <- t(as.matrix(newx))
    else newx <- as.matrix(newx)
    newx.norm <- (t(newx) - object$min.x)/object$range.x
    newx.norm <- t(as.matrix(newx.norm))
    pred.score <- drop(newx.norm%*%object$w) 
    u <- score2u(object$x.normalized,object$w,pred.score)$u
    B <- splines::bs(u,knots = quantile(eta, object$knots[-c(1,length(object$knots))]),
                     Boundary.knots = range(eta),degree = degree,intercept = F,warn.outside = F)
    pred.HR <- exp(drop(B%*%object$beta.B) - cen.hazard.log)
    pred.HR
  }
  
  HR.boot <- sapply(object.boot, function(x) predict.hazard.boot(x,newx = newx,cenx = cenx,degree = object$input_argue$degree))
  HR.ci <- t(apply(HR.boot, 1, quantile,prob = c(prob.par,0.5)))
  HR.point <- predict.hazard(object,newx = newx,cenx = cenx)
  HR <- cbind(HR.point$pred.score,HR.point$pred.HR,HR.ci)
  colnames(HR) <- c("score","HR.point",paste0("HR.",colnames(HR.ci)))
  
  rownames(HR) <- rownames(newx)
  list(HR = as.data.frame(HR),w = as.data.frame(w),newx = newx,cenx = cenx,cen.score = HR.point$cen.score,n.boot = nrow(w.boot))
}

#' @description This function calculate the instantaneous benefit, i.e., the benefit under
#'  a given small reduction in a specific air pollutant
#' @param object
#' @param x A vector or matrix indicating the exposure variables at which the derivation is conducted
#' @param percent logical. If true, indicating the percent decrease
snwsm_derivation <- function(object,x,unit = rep(1,length(object$w)),percent = F){
  if(is.vector(x)) origin_x <- t(as.matrix(x))
  else origin_x <- as.matrix(x)
  
  x <- t((t(origin_x) - object$min.x)/object$range.x) # 标化
  
  eta <- pu(object$x.normalized,object$w)$u
  x.score0 <- drop(x%*%object$w)
  range_score <- range(object$fit.score)
  is_score_out <- (range_score[2] < x.score0) | (range_score[1] > x.score0)
  x.score <- x.score0
  x.score[range_score[2] < x.score0] <- range_score[2]
  x.score[range_score[1] > x.score0] <- range_score[1]
  x.u <- score2u(object$x.normalized,object$w,x.score)
  B_deriv <- bsplineS(x.u$u, breaks = quantile(eta, object$knots),norder = object$input_argue$degree + 1,nderiv = 1)[,-1]
  aa <- object$w/object$range.x * unit
  if(percent) {
    derivs <- t(t(origin_x) * aa)/100
  } else{
    derivs <- matrix(rep(aa,nrow(x)),ncol = ncol(x),byrow = T) 
  }
  derivs <- drop((B_deriv * (x.u$deriv)) %*% object$beta.B) * derivs
  derivs <- exp(derivs)
  colnames(derivs) <- colnames(x)
  rownames(derivs) <- rownames(x)
  list(derivs = drop(derivs),score = x.score0,is_score_out = is_score_out)
}



#' @description This function implements the estimation of survival model with non-negative constraint
snwsm_surv_line <- function(time,censor,x,z = NULL,strata = NULL,input_save = F){
  if(input_save){
    call_params <- formals(sys.function(sys.parent()))
    # call_params <- as.list(match.call())[-1] # 这个不返回默认参数
    input_argue <- sapply(names(call_params),get,envir = environment())
  } else input_argue <- NULL
  if(anyNA(time)|anyNA(censor)|anyNA(x)) stop("Please remove the NA values before building the model")
  if(is.null(strata)) strata <- rep(1,length(time))
  if(!is.null(z)){
    if(anyNA(z)) stop("Please remove the NA values in Z before building the model")
    if("data.frame" %in% class(z)) z <- as.matrix(model.matrix(~.,z)[,-1])
    else z <- as.matrix(model.matrix(~ z)[,-1])  
    pz <- dim(z)[2] 
  } else pz <- 0 
  x <- as.matrix(x)
  if(dim(x)[2] == 1) stop("x must be a matrix with the number of columns larger than 1")
  origin_x <- x
  range_x <- apply(x, 2, function(a) max(a)-min(a))
  min_x <- apply(x, 2, min)
  
  x <- apply(x, 2, function(a) (a-min(a))/(max(a)-min(a)))
  n <- dim(x)[1]
  p <- dim(x)[2]
  nevent <- sum(1-censor)
  # browser()
  
  # beta_ini 
  xz <- cbind(x,z)
  fit_w <- coxph(Surv(time,1-censor)~xz + strata(strata))
  ab <- summary(fit_w)$coefficients[1:p,c(1,5)]
  a <- ab[,1]
  b <- ab[,2]
  while(any(a<0)){
    e <- match(min(b[a<0]),b)
    a[e] <- 0
    xz <- cbind(x[,which(a != 0)],z)
    fit_w <- coxph(Surv(time,1-censor)~xz + strata(strata))
    b[a != 0] <- summary(fit_w)$coefficients[1:sum(a != 0),5]
    a[a != 0] <- summary(fit_w)$coefficients[1:sum(a != 0),1]
  }
  
  beta_ini <- a
  w <- beta_ini[1:p]
  if(!all(w == 0)) w <- w/sum(w)
  
  beta.x <- w/range_x
  df <- sum(w != 0) + pz
  AIC <- -logLik(fit_w)*2 + 2*df
  BIC <- -logLik(fit_w)*2 + log(nevent)*df
  attr(AIC,"df") <- attr(BIC,"df") <- df
  range.x <- range_x
  names(w) <- names(beta.x) <- names(range.x) <- colnames(origin_x)
  list(w = w,beta.x = beta.x,fit.w.final = fit_w,
       range.x = range.x,AIC = AIC, BIC = BIC,x.normalized = x,
       input_argue = input_argue,min.x = min_x,df = df)
} 

#' @description This function plot the nonlinear nomogran according to snwsm_surv
#' @param object An object from snwsm_surv
#' @param filename A character with a suffix of ".tiff"
#' @param cenx The referenced values of x when calculating HR
#' @param varnames A vector of character
#' @param unitnames A vector of character
#' @param x_separate A matrix used to present the breaks for the monogram. 
#' The first row is the minimum label, the second row is the maximal label,
#' and the third row is the step. Each column corresponds to a factor.
#' @param total_score_lim The limited values for the total scores (total points)
#' @param hr_yaxis A list including three parameters to set the yaxis of the nonlinear curve. 
#' hr_lim: The maximal and minimal values of HR when plot the nonlinear association;
#' hr_labels: A vector of numeric. The labels needing to be presented at the hr axis;
#' hr_ticks: A vector of numeric. The ticks needing to be presented at the hr axis.
#' @param attri_num_plot A list including three parameters which indicates whether the attributable numbers of event were presented.
#' event_num_at_cenx: A numeric. The total numbers of events at cenx;
#' attri_num_labels: A vector of numeric. The lables for attributable numbers needing to be presented; 
#' attri_num_ticks: A vector of numeric. The ticks for for attributable numbers needing to be presented. 
#' attri_num_yname: A characher indicating the name of yaxis for the attributable events.
#' @param marks_x A list including four parameters which set vertical lines at the nonlinear curve.
#' marks_xvalue: A vector  or a matrix with a row indicating a vector of x needing to be marked.
#' marks_xname: A vector indicating the names of x needing to be marked
#' marks_xname_at_hr: A vector of numeric, indicating the position of xname needing to be marked,
#' marks_xcol: A vector indicating the colors of the vertical lines.
nnomogram.snwsm <- function(object,filename,cenx,varnames,unitnames,x_separate,total_score_lim,
                          segment_width,hr_yaxis,marks_x,attri_num_plot){
  # 命名下面向量
  w_names <- if(is.null(names(object$w))) paste0("w",1:length(object$w)) else  names(object$w)
  if(missing(varnames)) varnames <- w_names
  names(varnames) <- w_names
  if(!missing(unitnames)) names(unitnames) <- w_names
  
  x_min <- object$min.x
  x_max <- object$range.x + x_min
  nw <- length(object$w)
  range_x <- object$range.x
  
  score_max0 <- object$beta.x * object$range.x # 放大前最大贡献的wx
  magnify_score <- 100/max(score_max0) # 将wx放大的倍数
  score_max1 <- score_max0 * magnify_score # 放大后最大贡献的score/points
  range_score0 <- range(fit$fit.score) # 放大前的观测wx范围
  range_score1 <- range_score0 * magnify_score # 放大后的观测score范围
  
  
  if(!missing(marks_x)) {
    marks_x$marks_xvalue <- matrix(marks_x$marks_xvalue,ncol = nw)
    marks_wx <- apply(marks_x$marks_xvalue, 1, function(x){
      predict.hazard(object,score = range_score0[1],cenx = x)$cen.score
    })
  }
  
  predhr <- predict.hazard(object,score = seq(range_score0[1],range_score0[2],length = 100),cenx = cenx)
  # cen_wx <- predhr$cen.score
  
  # 将下面向量顺序化
  wnames_order <- names(sort(object$w,decreasing = T))
  point_var_names <- c("Point",varnames[wnames_order])
  if(!missing(unitnames)) unitnames <- unitnames[wnames_order]
  pvalues_names <- object$w_summary[paste0("xz",wnames_order),4] %>% round(3) %>% na.omit()
  pvalues_names <- sapply(pvalues_names, function(x){
    if(x < 0.001) substitute(italic(P)*" < "*xx,list(xx=0.001))
    else substitute(italic(P)*" = "*xx,list(xx=format(x,nsamll = 3))) 
  })
  pvalues_names <- do.call("expression", pvalues_names)
  
  x_min <- x_min[wnames_order] %>% round(2)
  x_max <- x_max[wnames_order] %>% round(2)
  ws <- fit$w[wnames_order]
  
  # 设置x_separate
  if(missing(x_separate)){
    a1 <- pretty(c(x_min[1],x_max[1]),n = 12)
    a2 <- length(a1)
    a1 <- a1[a1 >= x_min[1] & a1 <= x_max[1]]
    x_separate <- c(min(a1),max(a1),a1[2] - a1[1])
    for (i in 2:length(ws)) {
      a1 <- pretty(c(x_min[i],x_max[i]),n = round(a2*ws[i]/ws[1]))
      a1 <- a1[a1 >= x_min[i] & a1 <= x_max[i]]
      x_separate <- cbind(x_separate,c(min(a1),max(a1),a1[2] - a1[1]))
    }
    colnames(x_separate) <- wnames_order
  } else{
    colnames(x_separate) <- w_names
    x_separate <- x_separate[,wnames_order]
  }  
  
  ymax <- length(point_var_names)+2
  ys <- seq(ymax,1,length = length(point_var_names))
  cex.axis0 = 1.5
  tck = 0.01
  cex.axis1 <- 1.2
  
  
  
  # 根据range(c(predhr$pred.HR_lower,predhr$pred.HR_upper))确定下面参数
  if(missing(hr_yaxis)){
    hr_lim <- range(c(predhr$pred.HR_lower,predhr$pred.HR_upper))
    hr_labels <- pretty(hr_lim)
    hr_lim <- c(hr_lim[1]-diff(hr_lim)/20,hr_lim[2]) # 将下限适当放宽一点
    hr_ticks <- seq(min(hr_labels),max(hr_labels),by = (hr_labels[2] - hr_labels[1])/5)
    hr_labels <- hr_labels[hr_labels >= hr_lim[1] & hr_labels <= hr_lim[2]]
    hr_ticks <- hr_ticks[hr_ticks >= hr_lim[1] & hr_ticks <= hr_lim[2]]
  } else{
    hr_labels <- hr_yaxis$hr_labels
    hr_lim <- hr_yaxis$hr_lim
    hr_ticks <- hr_yaxis$hr_ticks
  }
  
  if(!missing(attri_num_plot)){
    attri_labels_hr <- attri_num_plot$attri_num_labels/attri_num_plot$event_num_at_cenx + 1
    attri_labels_tick_hr <- attri_num_plot$attri_num_ticks/attri_num_plot$event_num_at_cenx + 1
  }
  hr_plot <- do.call("cbind",predhr[1:4]) %>% as.data.frame()
  hr_point <- predhr$pred.HR
  hr_upper <- predhr$pred.HR_upper
  hr_lower <- predhr$pred.HR_lower
  curve_high <- 5 # 为曲线图分配的高度
  if(missing(total_score_lim)) total_score_lim <- c(0,ceiling(range_score1[2])) # 根据range_score1确定 @@@@@@@@@@
  total_score_labels <- seq(total_score_lim[1],total_score_lim[2],by = 10) # 根据total_score_lim确定 @@@@@@@@@@
  total_score_ticks <- seq(total_score_lim[1],total_score_lim[2],by = 2) # 根据total_score_lim确定 @@@@@@@@@@
  
  tiff(filename,width = 15,height = 20,units = "cm",pointsize = 6,res = 300)
  
  # plot 列线图
  plot(x = c(-10,110),y = c(-10,ymax),type = "n",axes = F,xlab = "",ylab = "")
  segments(x0 = seq(0,100,by = 2), y0 = max(ys), x1 = seq(0,100,by = 2), y1 = min(ys), col = "lightgray", lty = "dotted") #添加网格线
  text(-13,ys,labels = point_var_names,cex = 2,adj = 0,col = c("red",rep("black",length(point_var_names)-1)))
  if(!missing(unitnames)) text(-13,ys[-1]-0.5,labels = unitnames,cex = 1.2,adj = 0)
  text(101,ys[(1:length(pvalues_names))+1],labels = pvalues_names,cex = 1.1,adj = 0,col = "black", family = "sans")
  # padj控制label到坐标的上下距离，hadj控制左右距离
  axis(side = 3,at = seq(0,100,10),tick = T,pos = ymax,tck = -tck,cex.axis = cex.axis0,col = "red",padj = 3)
  axis(side = 3,at = seq(0,100,2),labels = F,tick = T,pos = ymax,tck = -tck/2,col = "red")
  
  for (i in 1:length(varnames)) {
    if(ws[i] == 0){
      text(x = 0,y = ys[i+1],labels = "No adverse effect was observed",cex = 1.5,pos = 4)
    }else{
      a <- x_min[i];b <- x_max[i]; d <- x_separate[,i];e <- ws[i]
      at_labels <- seq(d[1],d[2],by = d[3])
      at_score <- (at_labels-a)/(b-a)*e*magnify_score
      at_score_tick <- (seq(d[1],d[2],by = d[3]/5)-a)/(b-a)*e*magnify_score
      at_score_bound <- (c(a,b)-a)/(b-a)*e*magnify_score
      axis(side = 3,at = at_score_bound,tick = T,pos = ys[i+1],tck = -tck,lwd = 0.8,cex.axis = cex.axis1,labels = c(a,b),padj = 0,col.ticks = "blue",col.axis = "blue") 
      axis(side = 3,at = at_score,tick = T,pos = ys[i+1],tck = -tck,lwd = 0.8,cex.axis = cex.axis1,labels = at_labels,padj = 3.5) 
      axis(side = 3,at = at_score_tick,tick = T,pos = ys[i+1],tck = -tck/2,cex.axis = cex.axis1,labels = F,lwd = 0.8)
    }
  }
  
  
  # 下面画曲线
  total_score2at <- function(x) (x-total_score_lim[1])*100/(total_score_lim[2]-total_score_lim[1])
  hr2at <- function(x) min(ys)-1.3-(hr_lim[2]-x)*curve_high/(hr_lim[2]-hr_lim[1]) 
  
  if(missing(segment_width)) segment_width <- c(diff(total_score_labels[1:2])/2,diff(hr_labels[1:2])/2)
  a <- seq(0,total_score_lim[2],by = segment_width[1])
  b <- seq(min(hr_labels)-diff(hr_labels[1:2]),hr_lim[2],by = segment_width[2])
  b <- b[b > hr_lim[1] & b < hr_lim[2]]
  segments(x0 = total_score2at(a), y0 = hr2at(hr_lim[1]), x1 = total_score2at(a), y1 = hr2at(hr_lim[2]), col = "lightgray", lty = "dotted") #添加网格线
  segments(x0 = 0, y0 = hr2at(b), x1 = total_score2at(total_score_lim[2]), y1 = hr2at(b), col = "lightgray", lty = "dotted") #添加网格线
  
  at_hr <- hr2at(hr_labels)
  at_hr_tick <- hr2at(hr_ticks)
  at_total_score <- total_score2at(total_score_labels)
  at_total_score_tick <- total_score2at(total_score_ticks)
  
  axis(side = 2,at = at_hr,tick = T,pos = 0,tck = tck,labels = hr_labels,cex.axis = 1.2)
  axis(side = 2,at = at_hr_tick,tick = T,pos = 0,tck = tck/2,labels = F)
  arrows(x0 = 0,y0 = hr2at(hr_lim[1]),x1 = 0,y1 = hr2at(hr_lim[2])+curve_high/15,length = 0.05)
  
  axis(side = 3,at = at_total_score,tick = T,pos = hr2at(hr_lim[1]),tck = tck,labels = total_score_labels,cex.axis = 1.2,padj = 4)
  axis(side = 3,at = at_total_score_tick,tick = T,pos = hr2at(hr_lim[1]),tck = tck/2,labels = F,cex.axis = 1)
  arrows(x0 = 0,y0 = hr2at(hr_lim[1]),x1 = 102,y1 = hr2at(hr_lim[1]),length = 0.05) 
  
  if(!missing(attri_num_plot)){
    text(108,hr2at(mean(hr_lim)),labels = attri_num_plot$attri_num_yname,cex = 1.5,srt = 90)
    axis(side = 2,at = hr2at(attri_labels_hr),tick = T,pos = 100,tck = tck,labels = attri_num_plot$attri_num_labels,cex.axis = 1.2,padj = 4)
    axis(side = 2,at = hr2at(attri_labels_tick_hr),tick = T,pos = 100,tck = tck/2,labels = F)
    arrows(x0 = 100,y0 = hr2at(hr_lim[1]),x1 = 100,y1 = hr2at(hr_lim[2])+0.5,length = 0.05)
  }
  
  a <- total_score2at(hr_plot$pred.score*magnify_score)
  b <- hr2at(hr_plot$pred.HR)
  b_upper <- hr2at(hr_plot$pred.HR_upper)
  b_upper[b_upper>hr2at(hr_lim[2])] <- hr2at(hr_lim[2])
  b_lower <- hr2at(hr_plot$pred.HR_lower)
  lines(a,b,col = "red")
  redtrans <- rgb(255, 0, 0, 60, maxColorValue=255)
  polygon(c(rev(a),a),c(rev(b_upper),b_lower),col=redtrans, border = NA)
  lines(c(0,100),rep(hr2at(1),2),lty = 2)
  text(50,hr2at(hr_lim[1]-diff(hr_lim)/5),labels = "Total Points",cex = 1.5)  # @@@@@@@@@@@@@@
  text(-8,hr2at(mean(hr_lim)),labels = "HR",cex = 1.5,srt = 90)
  
  if(!missing(marks_x)){
    for (i in 1:length(marks_x$marks_xname)) {
      d <- total_score2at(marks_wx[i]*magnify_score)
      lines(c(d,d),hr2at(hr_lim),lty = 2,col = marks_x$marks_xcol[i])
      b <- paste0(marks_x$marks_xname[i]," (",round(marks_wx[i]*magnify_score,1),")")
      text(d,hr2at(marks_x$marks_xname_at_hr[i]),b,pos = 4,cex = 1.3,col = marks_x$marks_xcol[i])     # @@@@@@@@@@@@
    }
  }
  dev.off()
  RIPs <- ws * magnify_score
  RIPs
}


nnomogram.boot <- function(object,boot_CI,filename,varnames,unitnames,x_separate,total_score_lim,
                               segment_width,hr_yaxis,marks_x,attri_num_plot){
  # 命名下面向量
  w_names <- if(is.null(names(object$w))) paste0("w",1:length(object$w)) else  names(object$w)
  if(missing(varnames)) varnames <- w_names
  names(varnames) <- w_names
  if(!missing(unitnames)) names(unitnames) <- w_names
  
  x_min <- object$min.x
  x_max <- object$range.x + x_min
  nw <- length(object$w)
  range_x <- object$range.x
  
  score_max0 <- object$beta.x * object$range.x # 放大前最大贡献的wx
  magnify_score <- 100/max(score_max0) # 将wx放大的倍数
  score_max1 <- score_max0 * magnify_score # 放大后最大贡献的score/points
  range_score0 <- range(fit$fit.score) # 放大前的观测wx范围
  range_score1 <- range_score0 * magnify_score # 放大后的观测score范围
  
  
  if(!missing(marks_x)) {
    marks_x$marks_xvalue <- matrix(marks_x$marks_xvalue,ncol = nw)
    marks_wx <- apply(marks_x$marks_xvalue, 1, function(x){
      predict.hazard(object,score = range_score0[1],cenx = x)$cen.score
    })
  }
  
  aa <- boot_CI$HR[,c("score","HR.point","HR.2.5%","HR.97.5%")] %>% as.data.frame()
  predhr <- aa[aa$score >= min(object$fit.score) & aa$score <= max(object$fit.score),]
  colnames(predhr) <- c("pred.score","pred.HR","pred.HR_lower","pred.HR_upper")
  
  
  # 将下面向量顺序化
  wnames_order <- names(sort(object$w,decreasing = T))
  ws <- fit$w[wnames_order]
  point_var_names <- c("Point",varnames[wnames_order])
  if(!missing(unitnames)) unitnames <- unitnames[wnames_order]
  pvalues_names <- boot_CI$w[wnames_order,"w.pvalue"] %>% format(nsamll = 3) 
  pvalues_names <- pvalues_names[ws != 0]
  pvalues_names <- sapply(pvalues_names, function(x){
    if(x == 0) substitute(italic(P)*" < "*xx,list(xx=1/boot_CI$n.boot))
    else substitute(italic(P)*" = "*xx,list(xx=x)) 
  })
  pvalues_names <- do.call("expression", pvalues_names)
  
  x_min <- x_min[wnames_order] %>% round(2)
  x_max <- x_max[wnames_order] %>% round(2)
  
  
  # 设置x_separate
  if(missing(x_separate)){
    a1 <- pretty(c(x_min[1],x_max[1]),n = 12)
    a2 <- length(a1)
    a1 <- a1[a1 >= x_min[1] & a1 <= x_max[1]]
    x_separate <- c(min(a1),max(a1),a1[2] - a1[1])
    for (i in 2:length(ws)) {
      a1 <- pretty(c(x_min[i],x_max[i]),n = round(a2*ws[i]/ws[1]))
      a1 <- a1[a1 >= x_min[i] & a1 <= x_max[i]]
      x_separate <- cbind(x_separate,c(min(a1),max(a1),a1[2] - a1[1]))
    }
    colnames(x_separate) <- wnames_order
  } else{
    colnames(x_separate) <- w_names
    x_separate <- x_separate[,wnames_order]
  }  
  
  ymax <- length(point_var_names)+2
  ys <- seq(ymax,1,length = length(point_var_names))
  cex.axis0 = 1.5
  tck = 0.01
  cex.axis1 <- 1.2
  
  
  
  if(missing(hr_yaxis)){
    hr_lim <- range(c(predhr$pred.HR_lower,predhr$pred.HR_upper))
    hr_labels <- pretty(hr_lim)
    hr_lim <- c(hr_lim[1]-diff(hr_lim)/20,hr_lim[2]) # 将下限适当放宽一点
    hr_ticks <- seq(min(hr_labels),max(hr_labels),by = (hr_labels[2] - hr_labels[1])/5)
    hr_labels <- hr_labels[hr_labels >= hr_lim[1] & hr_labels <= hr_lim[2]]
    hr_ticks <- hr_ticks[hr_ticks >= hr_lim[1] & hr_ticks <= hr_lim[2]]
  } else{
    hr_labels <- hr_yaxis$hr_labels
    hr_lim <- hr_yaxis$hr_lim
    hr_ticks <- hr_yaxis$hr_ticks
  }
  
  if(!missing(attri_num_plot)){
    attri_labels_hr <- attri_num_plot$attri_num_labels/attri_num_plot$event_num_at_cenx + 1
    attri_labels_tick_hr <- attri_num_plot$attri_num_ticks/attri_num_plot$event_num_at_cenx + 1
  }
  hr_plot <- predhr
  hr_point <- predhr$pred.HR
  hr_upper <- predhr$pred.HR_upper
  hr_lower <- predhr$pred.HR_lower
  
  curve_high <- 5 # 为曲线图分配的高度
  if(missing(total_score_lim)) total_score_lim <- c(0,ceiling(range_score1[2])) # 根据range_score1确定 @@@@@@@@@@
  total_score_labels <- seq(total_score_lim[1],total_score_lim[2],by = 10) # 根据total_score_lim确定 @@@@@@@@@@
  total_score_ticks <- seq(total_score_lim[1],total_score_lim[2],by = 2) # 根据total_score_lim确定 @@@@@@@@@@
  
  tiff(filename,width = 15,height = 20,units = "cm",pointsize = 6,res = 300)
  
  # plot 列线图
  plot(x = c(-10,110),y = c(-10,ymax),type = "n",axes = F,xlab = "",ylab = "")
  segments(x0 = seq(0,100,by = 2), y0 = max(ys), x1 = seq(0,100,by = 2), y1 = min(ys), col = "lightgray", lty = "dotted") #添加网格线
  text(-13,ys,labels = point_var_names,cex = 2,adj = 0,col = c("red",rep("black",length(point_var_names)-1)))
  if(!missing(unitnames)) text(-13,ys[-1]-0.5,labels = unitnames,cex = 1.2,adj = 0)
  text(101,ys[(1:length(pvalues_names))+1],labels = pvalues_names,cex = 1.1,adj = 0,col = "black", family = "sans")
  # padj控制label到坐标的上下距离，hadj控制左右距离
  axis(side = 3,at = seq(0,100,10),tick = T,pos = ymax,tck = -tck,cex.axis = cex.axis0,col = "red",padj = 3)
  axis(side = 3,at = seq(0,100,2),labels = F,tick = T,pos = ymax,tck = -tck/2,col = "red")
  
  for (i in 1:length(varnames)) {
    if(ws[i] == 0){
      text(x = 0,y = ys[i+1],labels = "No adverse effect was observed",cex = 1.5,pos = 4)
    }else{
      a <- x_min[i];b <- x_max[i]; d <- x_separate[,i];e <- ws[i]
      at_labels <- seq(d[1],d[2],by = d[3])
      at_score <- (at_labels-a)/(b-a)*e*magnify_score
      at_score_tick <- (seq(d[1],d[2],by = d[3]/5)-a)/(b-a)*e*magnify_score
      at_score_bound <- (c(a,b)-a)/(b-a)*e*magnify_score
      axis(side = 3,at = at_score_bound,tick = T,pos = ys[i+1],tck = -tck,lwd = 0.8,cex.axis = cex.axis1,labels = c(a,b),padj = 0,col.ticks = "blue",col.axis = "blue") 
      axis(side = 3,at = at_score,tick = T,pos = ys[i+1],tck = -tck,lwd = 0.8,cex.axis = cex.axis1,labels = at_labels,padj = 3.5) 
      axis(side = 3,at = at_score_tick,tick = T,pos = ys[i+1],tck = -tck/2,cex.axis = cex.axis1,labels = F,lwd = 0.8)
    }
  }
  
  
  # 下面画曲线
  total_score2at <- function(x) (x-total_score_lim[1])*100/(total_score_lim[2]-total_score_lim[1])
  hr2at <- function(x) min(ys)-1.3-(hr_lim[2]-x)*curve_high/(hr_lim[2]-hr_lim[1]) 
  
  if(missing(segment_width)) segment_width <- c(diff(total_score_labels[1:2])/2,diff(hr_labels[1:2])/2)
  a <- seq(0,total_score_lim[2],by = segment_width[1])
  b <- seq(min(hr_labels)-diff(hr_labels[1:2]),hr_lim[2],by = segment_width[2])
  b <- b[b > hr_lim[1] & b < hr_lim[2]]
  segments(x0 = total_score2at(a), y0 = hr2at(hr_lim[1]), x1 = total_score2at(a), y1 = hr2at(hr_lim[2]), col = "lightgray", lty = "dotted") #添加网格线
  segments(x0 = 0, y0 = hr2at(b), x1 = total_score2at(total_score_lim[2]), y1 = hr2at(b), col = "lightgray", lty = "dotted") #添加网格线
  
  at_hr <- hr2at(hr_labels)
  at_hr_tick <- hr2at(hr_ticks)
  at_total_score <- total_score2at(total_score_labels)
  at_total_score_tick <- total_score2at(total_score_ticks)
  
  axis(side = 2,at = at_hr,tick = T,pos = 0,tck = tck,labels = hr_labels,cex.axis = 1.2)
  axis(side = 2,at = at_hr_tick,tick = T,pos = 0,tck = tck/2,labels = F)
  arrows(x0 = 0,y0 = hr2at(hr_lim[1]),x1 = 0,y1 = hr2at(hr_lim[2])+curve_high/15,length = 0.05)
  
  axis(side = 3,at = at_total_score,tick = T,pos = hr2at(hr_lim[1]),tck = tck,labels = total_score_labels,cex.axis = 1.2,padj = 4)
  axis(side = 3,at = at_total_score_tick,tick = T,pos = hr2at(hr_lim[1]),tck = tck/2,labels = F,cex.axis = 1)
  arrows(x0 = 0,y0 = hr2at(hr_lim[1]),x1 = 102,y1 = hr2at(hr_lim[1]),length = 0.05) 
  
  if(!missing(attri_num_plot)){
    text(108,hr2at(mean(hr_lim)),labels = attri_num_plot$attri_num_yname,cex = 1.5,srt = 90)
    axis(side = 2,at = hr2at(attri_labels_hr),tick = T,pos = 100,tck = tck,labels = attri_num_plot$attri_num_labels,cex.axis = 1.2,padj = 4)
    axis(side = 2,at = hr2at(attri_labels_tick_hr),tick = T,pos = 100,tck = tck/2,labels = F)
    arrows(x0 = 100,y0 = hr2at(hr_lim[1]),x1 = 100,y1 = hr2at(hr_lim[2])+0.5,length = 0.05)
  }
  
  a <- total_score2at(hr_plot$pred.score*magnify_score)
  b <- hr2at(hr_plot$pred.HR)
  b_upper <- hr2at(hr_plot$pred.HR_upper)
  b_upper[b_upper>hr2at(hr_lim[2])] <- hr2at(hr_lim[2])
  b_lower <- hr2at(hr_plot$pred.HR_lower)
  lines(a,b,col = "red")
  redtrans <- rgb(255, 0, 0, 60, maxColorValue=255)
  polygon(c(rev(a),a),c(rev(b_upper),b_lower),col=redtrans, border = NA)
  lines(c(0,100),rep(hr2at(1),2),lty = 2)
  text(50,hr2at(hr_lim[1]-diff(hr_lim)/5),labels = "Total Points",cex = 1.5)  # @@@@@@@@@@@@@@
  text(-8,hr2at(mean(hr_lim)),labels = "HR",cex = 1.5,srt = 90)
  
  if(!missing(marks_x)){
    for (i in 1:length(marks_x$marks_xname)) {
      d <- total_score2at(marks_wx[i]*magnify_score)
      lines(c(d,d),hr2at(hr_lim),lty = 2,col = marks_x$marks_xcol[i])
      b <- paste0(marks_x$marks_xname[i]," (",round(marks_wx[i]*magnify_score,1),")")
      text(d,hr2at(marks_x$marks_xname_at_hr[i]),b,pos = 4,cex = 1.3,col = marks_x$marks_xcol[i])     # @@@@@@@@@@@@
    }
  }
  dev.off()
  RIPs <- ws * magnify_score
  RIPs
}

nnomogram <- function(object,boot_CI = NULL,filename,cenx = NULL,varnames,unitnames,x_separate,total_score_lim,
                           segment_width,hr_yaxis,marks_x,attri_num_plot){
  if(is.null(boot_CI)){
    nnomogram.snwsm(object,filename,cenx,varnames,unitnames,x_separate,total_score_lim,
                    segment_width,hr_yaxis,marks_x,attri_num_plot)
  } else{
    nnomogram.boot(object,boot_CI,filename,varnames,unitnames,x_separate,total_score_lim,
                   segment_width,hr_yaxis,marks_x,attri_num_plot)
  }
}