# NOTE: this a modified version of the original simulation file. To obtain
# the results of our simulations, we utilized a computer cluster. They are
# too expensive to run on 1 machine. This file is representative of what was
# done in our simulations. The script as written should run in a timely manner.

library(glmnet)
library(MASS)
source('simulation_support_funcs.R')

run1iter <- function(B = 100, p = 10, n = 100, r = 200, SNR = 4, rho = 0,
                     Pstar = 0.95, beta_type = 1, do_cvc = F){
  
  lasso_vars<-tfmat1<-fs_vars<-ssfs_vars<-tfmat4<-tfmat5<-NA
  
  ntrue <- 5
  mysigma <- diag(p)
  
  for(j in 1:p){
    for(k in 1:p){
      mysigma[j,k]<-rho^abs(j-k)
    }
  }
  
  if(beta_type == 1) {
    coefs <- rep(1, ntrue)
    s.covar <- round(seq(1, p, len = ntrue))
  } 
  if(beta_type == 2) {
    coefs <- rep(1, ntrue)
    s.covar <- 1:ntrue
  }
  if(beta_type == 3) {
    coefs <- seq(10, 0.5, len = ntrue)
    s.covar <- 1:ntrue
  } 
  
  betas <- rep(0,p)
  betas[s.covar] <- coefs
  betas <- matrix(betas,ncol=1)
  
  # generate training data
  n1<-rmvnorm(n=n,mean=rep(0,p),sigma=mysigma)
  colnames(n1)<-paste("x",1:p,sep="")
  signal <- (n1%*%betas)
  eps <- rnorm(n,0,sqrt(var(signal)/SNR))
  y <- signal + eps
  
  # generate test data
  n1.test <- rmvnorm(n = 100000, mean = rep(0, p), sigma = mysigma)
  colnames(n1.test)<-paste("x", 1:p, sep="")
  signal.test <- (n1.test %*% betas)
  eps.test <- rnorm(100000, 0, sqrt(var(signal) / SNR))
  y.test <- signal.test + eps.test
  
  # lasso
  model<-glmnet(x = n1, y = y)
  mycond<-which(colSums(ifelse(as.matrix(coef(model)[-1,])==0,0,1))<=ntrue)
  mysnum<-names(which.max(mycond))
  red.coefs<-coef(model)[-1,which(colnames(coef(model))==mysnum)]
  myranked<-rank(-abs(red.coefs))
  myrankedOr<-myranked[order(myranked)]
  lasso.sel<-	names(myrankedOr)[1:ntrue]
  lasso_vars<-lasso.sel
  pmat <- predict(model, n1.test)
  lasso_err <- mean((pmat[,colnames(pmat) == mysnum] - y.test)^2)
  
  # lasso cv
  mycv <- cv.glmnet(x = n1, y = y)
  lasso_cved <- glmnet(x = n1, y = y, lambda = mycv$lambda.min)
  red.coefs <- coef(lasso_cved)[-1,]
  lasso.sel <- names(which(red.coefs != 0))
  lasso_vars_cv <- lasso.sel
  lasso_err_cv <- mean((predict(lasso_cved, n1.test) - y.test)^2)
  
  # forward selection
  mydata<-data.frame(n1,y)
  lower<-lm(as.formula(paste("y","~","1",sep="")),mydata)
  upper<-lm(as.formula(paste("y","~",".",sep="")),mydata)
  mystep<-step(lower,scope=list(lower=lower,upper=upper),direction="forward",steps=ntrue,trace=0,k=0)
  fs_vars<-names(coef(mystep))[-1]
  
  # stability selection w/ forward selection
  mypicks<-c()
  for(j in 1:B){
    samp<-sample(1:n,n,replace=TRUE)
    train<-mydata[samp,]
    lower<-lm(as.formula(paste("y","~","1",sep="")),train)
    upper<-lm(as.formula(paste("y","~",".",sep="")),train)
    mystep<-step(lower,scope=list(lower=lower,upper=upper),direction="forward",steps=ntrue,trace=0,k=0)
    steps<-names(coef(mystep))[-1]
    mypicks<-c(mypicks,steps)
  }
  ssfs_vars<-names(table(mypicks)[order(table(mypicks),decreasing=T)][1:ntrue])
  
  # CVC
  full_cvc_res <- mps_cvc_res <- cvc_time <- mps_cvc_time <- NULL
  if(do_cvc) {
    mymat <- combn(1:p, ntrue)
    my_ols.models <- split(mymat, rep(1:ncol(mymat), each = nrow(mymat)))
    pt <- proc.time()
    mycvc <- CVC(n1, as.matrix(y), type = "ols", ols.models = my_ols.models, n.fold = as.numeric(args[9]), B = 100)
    cvc_time <- proc.time() - pt
    full_cvc_res <- matrix(paste0("x", do.call(rbind, my_ols.models[which(mycvc$p.vals.c > 0.05)])), ncol = ntrue)
    
    pt <- proc.time()
    picked_mods_ls <- mps_cvc_func(n1, as.matrix(y), ntrue, p)
    mps_cvc_time <- proc.time() - pt
    mps_cvc_res <- matrix(paste0("x", do.call(rbind, picked_mods_ls)), ncol = ntrue)
  }
  
  # MPS
  mse <- function(y,y.pred){return(mean((y-y.pred)^2))}
  pt <- proc.time()
  full_sel <- full.select.gen(myframe = data.frame(y,n1),
                              resp.name = "y", depth = ntrue, r = r,
                              Pstar = Pstar, condense = TRUE)
  mps_time <- proc.time() - pt
  myc <- full_sel$mps
  myc_counts <- full_sel$counts
  fss_vars <- myc[1,]
  mps_vars <- myc
  
  
  tfacc1<-fs_err<-ssfs_err<-tfacc4<-tfacc5<-NA
  fss_err<-mps_err<-cvc_err<-cvc_mps_err<-NA
  



  mod<-lm(as.formula(paste("y", paste(c(na.omit(fs_vars)), collapse=" + "), sep=" ~ ")),data.frame(cbind(n1,y)))
  pred<-predict(mod,data.frame(n1.test))
  fs_err<-mean((pred-y.test)^2)

  mod<-lm(as.formula(paste("y", paste(c(na.omit(ssfs_vars)), collapse=" + "), sep=" ~ ")),data.frame(cbind(n1,y)))
  pred<-predict(mod,data.frame(n1.test))
  ssfs_err<-mean((pred-y.test)^2)

  rte_cvc<-NULL
  if(do_cvc){
    for(i in 1:dim(full_cvc_res)[1]){
      mod<-lm(as.formula(paste("y", paste(c(na.omit(full_cvc_res[i,])), collapse=" + "), sep=" ~ ")),data.frame(cbind(n1,y)))
      pred<-predict(mod,data.frame(n1.test))
      cvc_err[i]<-mean((pred-y.test)^2)
    }
  }
  
  rte_mps_cvc<-NULL
  if(do_cvc){
    for(i in 1:dim(mps_cvc_res)[1]){
      mod<-lm(as.formula(paste("y", paste(c(na.omit(mps_cvc_res[i,])), collapse=" + "), sep=" ~ ")),data.frame(cbind(n1,y)))
      pred<-predict(mod,data.frame(n1.test))
      cvc_mps_err[i]<-mean((pred-y.test)^2)
    }
  }
  
  mod<-lm(as.formula(paste("y", paste(c(na.omit(fss_vars)), collapse=" + "), sep=" ~ ")),data.frame(cbind(n1,y)))
  pred<-predict(mod,data.frame(n1.test))
  fss_err<-mean((pred-y.test)^2)
  
  rte7<-NULL
  for(i in 1:dim(mps_vars)[1]){
    mod<-lm(as.formula(paste("y", paste(c(na.omit(mps_vars[i,])), collapse=" + "), sep=" ~ ")),data.frame(cbind(n1,y)))
    pred<-predict(mod,data.frame(n1.test))
    mps_err[i]<-mean((pred-y.test)^2)
  }
  
  mod<-lm(as.formula(paste("y", paste(c(colnames(n1)[s.covar]), collapse=" + "), sep=" ~ ")),data.frame(cbind(n1,y)))
  pred<-predict(mod,data.frame(n1.test))
  tacc<-mean((pred-y.test)^2)
  
  pred_true <- n1.test %*% betas
  tacc_true <- mean((pred_true - y.test)^2)
  
  return(list(truecov = colnames(n1)[s.covar],
              lasso = lasso_vars, lasso_cv = lasso_vars_cv,
              fs = fs_vars, ssfs = ssfs_vars,
              cvc = full_cvc_res, mps_cvc = mps_cvc_res,
              fss = fss_vars, mps = mps_vars, mps.counts = myc_counts,
              true_model_err = tacc_true, oracle_err = tacc,
              lasso_err = lasso_err, lasso_err_cv = lasso_err_cv,
              fs_err = fs_err, ssfs_err = ssfs_err,
              cvc_err = cvc_err, cvc_mps_err = cvc_mps_err,
              fss_err = fss_err, mps_err = mps_err, 
              cvc_time = cvc_time, mps_time = mps_time, mps_cvc_time = mps_cvc_time))
}

results1 <- run1iter()

# to obtain RTEs as done in our paper, here is an example
get_rte <- function(x) {
  x / results1$true_model_err
}
get_rte(results1$lasso_err_cv)


