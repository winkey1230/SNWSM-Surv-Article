##############################################################################-#
## The codes give an example for implementing snwsm_surv_fin                  -#
## 2024/12/10                                                                 -#
##############################################################################-#
library(sigmoid)
library(mvtnorm)
source("snwsm_surv_fun.R")
`%>%` <- magrittr::`%>%`
################### generating artificial data------------
set.seed(5)
p <- 5
n <- 1000
# generate X
sigma <- outer(1:p, 1:p, function(i, j){ 2^(-abs(i-j)) })
sigma <- outer(1:p, 1:p, function(i, j){ 2^(-abs((i-j)/5)) })
X <- rmvnorm(n, mean = rep(0,p), sigma = sigma * 0.5)
X <- relu(X + 1) - 1
X <- 1 - relu(1 - X)

sigma <- outer(1:2, 1:2, function(i, j){ 2^(-abs(i-j)) })
Z <- rmvnorm(n, mean = rep(0,2), sigma = sigma * 0.5)
Z <- relu(Z + 1) - 1
Z <- 1 - relu(1 - Z)
cor(X)
# generate the true association
beta1 <- 0:4/sqrt(sum((0:4)^2))
beta1 <- 0:4/sum(0:4)*2
beta2 <- c(0.5,0.5)
mu1 <- X %*% beta1
mu2 <- Z %*% beta2 # confounding effects
ref.hr <- exp(-log(Gompertz(-sum(apply(X,2,median)*beta1),a = 4,b = 1, c = 1)+0.1))
g1 <- -log(Gompertz(-mu1,a = 4,b = 1, c = 1)+0.1); plot(mu1,exp(g1)/ref.hr,type = "p") # the true association
# plot(mu1,g1,type = "p")
g2 <- mu2

# generate y
lambda0 <- 1
lambda <- exp(g1 + g2) * lambda0
S_time <- rexp(n,lambda)
C_time <- runif(n,0,2)
quantile(S_time)
time <- pmin(S_time,C_time)
censor <- as.integer(S_time > C_time) 
data <- data.frame(time,censor)
mean(censor)
xnames <- paste0("X",1:5)
znames <- paste0("Z",1:2)
colnames(X) <- xnames
colnames(Z) <- znames
rm(lambda0,lambda,sigma,C_time,censor,xnames,znames,S_time,time)
################### fit snwsm model -----
# fit model
fit <- snwsm_surv(time = data$time,censor = data$censor,x = X,z = Z,strata = NULL,
                 knots = c(0.1,0.5,0.9),input_save = T,max.iter = 20,degree = 4)
fit$w

# compare the estiamted HR with the true HR
ref <- predict.hazard(fit,newx = apply(X,2,median),cenx = apply(X,2,median))
## The fitted value
plot(mu1,exp(fit$fit.hazard.log-ref$cen.hazard.log),xlab = "w*x",ylab = "HR")
## the true value
points(mu1,exp(g1)/ref.hr,col = "red") 
## The fitted value
plot(fit$fit.score,exp(fit$fit.hazard.log-ref$cen.hazard.log),xlab = "w*x",ylab = "HR") 
## the true value
points(fit$fit.score,exp(g1)/ref.hr,col = "red") 

# using the bootstrap method to estimate 95%CI
cenx <- c(0.5,0.5,0.5,0.5,0.5) # the referenced values
fit.boot <- snwsm_bootstrap(fit,R = 1000,ncore = 4)
res.ci <- snwsm_CI(object = fit,object.boot = fit.boot,cenx = cenx)
hr <- res.ci$HR
hr <- hr[order(hr$score),]
w <- res.ci$w

# plot the association curve between w*x and HR
plot(hr$score,hr$HR.point,type = "l",ylim = range(c(hr$`HR.97.5%`,hr$`HR.2.5%`)))
lines(hr$score,hr$`HR.2.5%`,type = "l",col = "blue")
lines(hr$score,hr$`HR.97.5%`,type = "l",col = "blue")
abline(h = 1)

# plot nomogram
## the reference values of x
varnames <- paste0("X",1:length(fit$w))
marks_x <- list(marks_xvalues = rbind(cenx,rep(0.1,length(fit$w)), rep(-0.5,length(fit$w))),
                marks_xname = c("ref0","ref1","ref2"),
                marks_xname_at_hr = c(0.6,2,4),
                marks_xcols = c("blue","red","orange"))
## the number of event at cenx was set as 2000
attri_num_plot <- list(event_num_at_cenx = 2000,attri_num_labels = seq(-2000,7000,by = 1000),
                       attri_num_ticks = seq(-2000,8000,by = 200),attri_num_yname = "Attributable No.")
total_score_lim = c(0,252)
## plot and get RIPs
RIPs <- nnomogram(object = fit,boot_CI = res.ci,filename = "nomogram.tiff",cenx = cenx,
                       marks_x = marks_x,attri_num_plot = attri_num_plot,varnames = varnames,
                       total_score_lim = total_score_lim)
## plot using the default settings and get RIPS
RIPs <- nnomogram(object = fit,boot_CI = res.ci,filename = "nomogram.tiff",cenx = cenx)



# # using bootstrap to estimate the uncertainty will cost intensive computation sources
# # an alternative method is to use the uncertainty from the sign-constrained Cox model
# # but this method may underestimate the uncertainty
# fit$w_summary
# predhr <- predict.hazard(fit,cenx = cenx)
# plot(predhr$pred.score,predhr$pred.HR,type = "l",xlab = "w*x",ylab = "HR",
#      ylim = range(c(predhr$pred.HR_lower,predhr$pred.HR_upper)))
# lines(predhr$pred.score,predhr$pred.HR_lower,type = "l",col = "red")
# lines(predhr$pred.score,predhr$pred.HR_upper,type = "l",col = "red")
# abline(h = 1)
# 
# # plot nomogram
# ## the reference values of x
# varnames <- paste0("X",1:length(fit$w))
# marks_x <- list(marks_xvalues = rbind(cenx,rep(0.1,length(fit$w)), rep(-0.5,length(fit$w))),
#                 marks_xname = c("ref0","ref1","ref2"),
#                 marks_xname_at_hr = c(0.6,2,4),
#                 marks_xcols = c("blue","red","orange"))
# ## the number of event at cenx was set as 2000
# attri_num_plot <- list(event_num_at_cenx = 2000,attri_num_labels = seq(-2000,7000,by = 1000),
#                        attri_num_ticks = seq(-2000,8000,by = 200),attri_num_yname = "Attributable No.")
# total_score_lim = c(0,252)
# ## plot and get RIPs
# RIPs <- nnomogram(object = fit,filename = "nomogram1.tiff",cenx = cenx,
#                        marks_x = marks_x,attri_num_plot = attri_num_plot,varnames = varnames,
#                        total_score_lim = total_score_lim)
# ## plot using the default settings and get RIPS
# RIPs <- nnomogram.snwsm(object = fit,filename = "nomogram1.tiff",cenx = cenx)




