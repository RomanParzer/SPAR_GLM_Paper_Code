###########################################
### benchmarks and competitors
##########################################

### overview: wrappers for functions providing predictions and set of active predictors (where applicable)
# 1 AdLASSO glmnet, 
# 2 glmnet with alpha=3/4, 
# 3 SIS and ISIS,
# 6 SPAR wrapper (incl RP + TARP)
# 7 RF
# 8 svm

# deprecated
# 4 TARP glm adaption
# 5 RP wrapper

# # install.packages("remotes")
# remotes::install_github("RomanParzer/SPAR@v3.2") # for specific release used in these simulations and applications
# # remotes::install_github("RomanParzer/SPAR@*release") # for latest release
# # remotes::install_github("RomanParzer/SPAR") # for current REPO status

pacman::p_load(pls, glmnet, SIS, MASS,SPAR,Matrix,ROCR,metrica,e1071,randomForest)

# # 1 AdLASSO
myAdLASSO <- function(x,y,xtest,family,type.measure="default") {
  if (family$family=="binomial" & family$link=="logit") {
    fit_family <- "binomial"
  } else if (family$family=="poisson" & family$link=="log") {
    fit_family <- "poisson"
  } else if (family$family=="gaussian" & family$link=="identity") {
    fit_family <- "gaussian"
  } else {
    fit_family <- family
  }
  ridge_cv <- cv.glmnet(x = x, y = y,nfold = 10,alpha = 0,family=fit_family,type.measure=type.measure)
  best_ridge_coef <- as.numeric(coef(ridge_cv, s = ridge_cv$lambda.min))[-1]
  alasso_cv <- cv.glmnet(x=x,y=y,penalty.factor = 1 / abs(best_ridge_coef),alpha=1,nfold=10,family=fit_family,type.measure=type.measure)
  # refit with optimal lambda
  alasso <- glmnet(x,y,penalty.factor = 1 / abs(best_ridge_coef),alpha=1, lambda = alasso_cv$lambda.1se,family=fit_family)
  
  yhat <- predict(alasso,s = alasso_cv$lambda.1se,newx = xtest, type = "response")
  yhat_tr <- predict(alasso,s = alasso_cv$lambda.1se,newx = x, type = "response")
  
  return(list(yhat=yhat, yhat_tr=yhat_tr,lambda=alasso_cv$lambda.1se,beta=alasso$beta,intercept=alasso$a0, scr_coef=alasso$beta))
}

# # 2 Elastic Net with alpha
myElNet <- function(x,y,xtest,family,alpha=3/4,type.measure="default") {
  if (family$family=="binomial" & family$link=="logit") {
    fit_family <- "binomial"
  } else if (family$family=="poisson" & family$link=="log") {
    fit_family <- "poisson"
  } else if (family$family=="gaussian" & family$link=="identity") {
    fit_family <- "gaussian"
  } else {
    fit_family <- family
  }
  elnet_cv <- cv.glmnet(x,y,alpha=alpha,folds=10,family=fit_family,type.measure=type.measure) # standardize = TRUE by default
  # refit with optimal lambda
  elnet <- glmnet(x,y,alpha=alpha, lambda = elnet_cv$lambda.1se,family=fit_family)
  yhat <- predict(elnet,s = elnet_cv$lambda.1se,newx = xtest, type = "response")
  yhat_tr <- predict(elnet,s = elnet_cv$lambda.1se,newx = x, type = "response")
  
  return(list(yhat=yhat, yhat_tr=yhat_tr,lambda=elnet_cv$lambda.1se,beta=elnet$beta,intercept=elnet$a0, scr_coef=elnet$beta))
}

# # 3 SIS 
mySIS <- function(x,y,xtest,family,iter=FALSE,type.measure="deviance") {

  invisible(capture.output(
    SIS_res <- SIS(x,y,tune="cv",nfolds=10,iter=iter,family=family$family,type.measure=type.measure)
  ))
  
  yhat <- predict(SIS_res,newx=xtest,type="response")
  yhat_tr <- predict(SIS_res,newx=x,type="response")
  
  beta <- Matrix(c(0),nrow=ncol(x),ncol=1,sparse=TRUE)
  beta[SIS_res$ix,1] <- SIS_res$coef.est[-1]
  
  return(list(yhat=yhat, yhat_tr=yhat_tr,lambda=SIS_res$lambda,ind_screen=SIS_res$ix0,
              beta=beta,intercept=SIS_res$coef.est[1],scr_coef=beta))
}

# # 6 SPAR
mySPAR <- function(x,y,xtest,family,doCV=FALSE,
                   nummods=c(20),nlambda=20,
                   split_data=FALSE, nscreen = 2*nrow(x),
                   type.measure="deviance",
                   type.rpm="cwdatadriven",
                   type.screening="ridge",
                   opt_par="best",
                   avg_type="link") {
  if (!doCV) {
    spar_res <- spar(x,y,family=family,nummods=nummods,nlambda = nlambda,
                     type.measure=type.measure,type.rpm=type.rpm,type.screening=type.screening,
                     control=list(scr=list(nscreen=nscreen, split_data=split_data)))
    val_sum <- spar_res$val_res
    coef <- coef(spar_res)
  } else {
    spar_res <- spar.cv(x,y,family=family,nummods=nummods,nlambda = nlambda,
                        type.measure=type.measure,type.rpm=type.rpm,type.screening=type.screening,
                        control=list(scr=list(nscreen=nscreen, split_data=split_data)))
    val_sum <- spar_res$val_sum
    coef <- coef(spar_res,opt_par = opt_par)
  }
  yhat <- predict(spar_res,xnew = xtest,coef = coef,avg_type=avg_type)
  yhat_tr <- predict(spar_res,xnew = x,coef = coef,avg_type=avg_type)
  
  return(list(yhat=yhat, yhat_tr=yhat_tr,beta=coef$beta,intercept=coef$intercept,
              betas=spar_res$betas,val_sum=val_sum,scr_coef = spar_res$scr_coef/spar_res$xscale[spar_res$xscale>0]))
}


# # 7 RF

myRF <- function(x,y,xtest,family) {
  n <- length(y)
  p <- ncol(x)
  
  if (family$family == "binomial") {
    y <- factor(y)
  } 

  capture.output(
    tRF <- tuneRF(x,y,ntreeTry = n, doBest=TRUE, stepFactor = 1.5, improve=0.01,trace=FALSE), 
    file = nullfile()
  )
  mybeta <- as.numeric(tRF$importance)
  
  if (family$family == "binomial") {
    preds <- as.numeric(predict(tRF,newdata = rbind(x,xtest),type="prob")[,2])
  } else {
    preds <- predict(tRF,newdata = rbind(x,xtest),type="response")
  }
  
  return(list(yhat=preds[-(1:n)],yhat_tr=preds[1:n],beta=mybeta))
}


# # 8 svm

mySVM <- function(x,y,xtest,family) {
  n <- length(y)
  p <- ncol(x)
  prob <- FALSE
  
  if (family$family == "binomial") {
    y <- factor(y)
    prob <- TRUE
  } 
  best_svm <- best.tune(svm, train.x=x,train.y=y,probability=prob,
                        ranges = list(cost = c(0.5,1,2),kernel=c("radial","linear"),
                        tunecontrol = tune.control(sampling = "cross",cross=10)
                        )
  )
  pred <- predict(best_svm,newdata = rbind(x,xtest),probability = prob)
  if (prob==TRUE) {
    col1_ind <- which(colnames(attributes(pred)$probabilities)=="1")
    pred <- as.numeric(attributes(pred)$probabilities[,col1_ind])
  }
  mybeta <- as.numeric(crossprod(best_svm$coefs,best_svm$SV))
  return(list(yhat=pred[-(1:n)], yhat_tr=pred[1:n],beta = mybeta))
}


# for Evaluation of Variable Ranking 
my_AUC <- function(ytest,phat) {
  if (length(phat)==0) {
    return(NA)
  }
  stopifnot(length(ytest)==length(phat))
  # ypred <- as.numeric(phat>0.5)
  if (var(ytest)==0) {
    res <- NA
  } else {
    res <- performance(prediction(phat,ytest),measure="auc")@y.values[[1]] 
  }
  return(res)
}

my_MCC <- function(ytest,phat) {
  if (length(phat)==0) {
    return(NA)
  }
  if (length(ytest)!=length(phat)) {
    warning("MCC input lengths differ!")
    return(NA)
  }
  ypred <- as.numeric(phat>0.5)
  if (var(ytest)==0 | var(ypred)==0) {
    res <- NA
  } else {
    # can lead to integer overflows, then implement by hand (see VarSel)
    res <- mcc(data=NULL,ytest,ypred)$mcc
  }
  return(res)
}

