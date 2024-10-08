
pacman::p_load(foreach,parallel,tidyr,dplyr,ROCR,SPAR,e1071)
source("../functions/glm_data_generation.R")
source("../functions/glm_methods.R")
source("../functions/multi_assign.R")

# compare prediction performance (AUC) and link estimation (MSLE) of different random projections
simulation_settings <- tibble(n=200, p=c(2000),ntest=200, act_setting="medium", cov_setting=c("group"),
                              signal_strength=c(10,1/8,100,1000,1/4), avg_exp = c(1,10,0.5,0.7,10),
                              family=list(gaussian(identity),gaussian(log),binomial(logit),binomial(cloglog),poisson(log)))

simulation_settings <- rbind(simulation_settings,
                             tibble(n=200, p=c(2000),ntest=200, act_setting="sparse", cov_setting=c("group"),
                                    signal_strength=c(10,1/8,100,1000,1/4), avg_exp = c(1,10,0.5,0.7,10),
                                    family=list(gaussian(identity),gaussian(log),binomial(logit),binomial(cloglog),poisson(log))))

simulation_settings <- rbind(simulation_settings,
                             tibble(n=200, p=c(2000),ntest=200, act_setting="dense", cov_setting=c("group"),
                                    signal_strength=c(10,1/8,100,1000,1/4), avg_exp = c(1,10,0.5,0.7,10),
                                    family=list(gaussian(identity),gaussian(log),binomial(logit),binomial(cloglog),poisson(log))))


simulation_settings <- rbind(simulation_settings,
                             tibble(n=200, p=c(500),ntest=200, act_setting="medium", cov_setting=c("group"),
                                    signal_strength=c(10,1/8,100,1000,1/4), avg_exp = c(1,10,0.5,0.7,10),
                                    family=list(gaussian(identity),gaussian(log),binomial(logit),binomial(cloglog),poisson(log))))

simulation_settings <- rbind(simulation_settings,
                             tibble(n=200, p=c(10000),ntest=200, act_setting="medium", cov_setting=c("group"),
                                    signal_strength=c(10,1/8,100,1000,1/4), avg_exp = c(1,10,0.5,0.7,10),
                                    family=list(gaussian(identity),gaussian(log),binomial(logit),binomial(cloglog),poisson(log))))

simulation_settings <- simulation_settings %>% mutate(a = round(ifelse(act_setting=="sparse",2*log(p),
                                                                       ifelse(act_setting=="medium",n/2+2*log(p),p/4))))

nset <- nrow(simulation_settings)
methods <- c("ElNet","AdLASSO","Gaussian","Sparse","SparseCW",
             "GLM_HOLP","Ridge_cv",
             "Ridge_sqrtnp","Ridge_10","Ridge_1","Ridge_01","Ridge_001",
             "Ridge_dev08","Ridge_dev095","Ridge_dev0999",
             "marGLM","True_Beta")
measures <- c("AUC","MSPE","MSLE","pAUC","Cor_scr","pAUC_scr","lambda","dev.ratio","mtr_over_mf_abs","TrRatio_3a")
nmeas <- length(measures)
nmethods <- length(methods)
nrep <- 2
res <- array(c(0),dim=c(nrep,nmeas,nmethods,nset),
             dimnames = list(rep=1:nrep,
                             measure=measures,
                             method=methods,
                             setting=1:nset))
attributes(res)$settings <- simulation_settings

unlink("../saved_results/log.txt")
my.cluster <- parallel::makeCluster(7, type = "PSOCK", outfile = "../saved_results/log.txt")
doParallel::registerDoParallel(cl = my.cluster)
# foreach::getDoParRegistered()
clusterExport(my.cluster,c('simulation_settings','nmethods','methods'), envir = environment())
clusterEvalQ(my.cluster, {  
  pacman::p_load(ROCR,SPAR,ggplot2,e1071)
  source("../functions/glm_data_generation.R")
  source("../functions/glm_methods.R")
  source("../functions/multi_assign.R")
})

# i <- j <- 3

parres <- foreach(j = 1:nset) %:%
  foreach(i=1:nrep) %dopar% {
  resij <- matrix(c(0),nmeas,nmethods)
  c(n,p,ntest,act_setting,cov_setting,signal_strength,avg_exp,family,a) %<-% simulation_settings[j,]
  family <- family[[1]]
  set.seed((1234+i)^2)
  data <- generate_data_glm(n,p,family=family,cov_setting=cov_setting,ntest=ntest,a=a,signal_strength = signal_strength, avg_exp = avg_exp)
  
  m <- n/4
  x <- data$x
  y <- data$y
  beta <- data$beta
  alpha <- data$alpha
  xtest <- data$xtest

  # # testing different scales internally
  # x[,1:(p/2)] <- x[,1:(p/2)]*2 + 3
  # xtest[,1:(p/2)] <- xtest[,1:(p/2)]*2 + 3
  # beta[1:(p/2)] <- beta[1:(p/2)]/2
  # alpha <- alpha - sum(3*beta[1:(p/2)])
  
  xcenter <- apply(x,2,mean)
  xscale <- apply(x,2,sd)
  # xcenter <- rep(0,p)
  # xscale <- rep(1,p)
  
  z <- scale(x,center=xcenter,scale=xscale)
  ztest <- scale(xtest,center = xcenter, scale = xscale)

  coefs <- matrix(c(0),p,nmethods-5)
  colnames(coefs) <- methods[-c(1:5)]
  
  tmp_sc <- apply(x,2,function(col)sqrt(var(col)*(n-1)/n))
  z2 <- scale(x,center=colMeans(x),scale=tmp_sc)
  lam_max <- 1000 * max(abs(t(y)%*%z2[,tmp_sc>0]))/n*family$mu.eta(family$linkfun(mean(y)))/family$variance(mean(y))
  
  if (family$family=="binomial") {
    y_svm <- factor(ifelse(y>0 ,1,0))
    best_svm <- best.tune(svm, train.x=z,train.y=y_svm, probability=TRUE,kernel="linear",ranges = list(cost = 1e10),
    ) 
    pred <- predict(best_svm,probability=TRUE,newdata=z)
    if (colnames(attributes(pred)$probabilities)[1]=="0") { # if goal class is 0 switch sign
      best_svm$coefs <- -best_svm$coefs
    }
    coefs[,1] <- as.numeric(crossprod(best_svm$coefs,best_svm$SV))
  } else {
    eig <- eigen(tcrossprod(z),symmetric = TRUE)
    myinv <- tcrossprod(eig$vectors[,eig$values>1e-8]%*%diag(1/sqrt(eig$values[eig$values>1e-8])))
    if (family$family=="poisson") {
      y_corr <- ifelse(y==0,1e-4,y)
    } else {
      y_corr <- y
    }
    coefs[,1] <- as.numeric(crossprod(z,myinv%*%(family$linkfun(y_corr) - mean(family$linkfun(y_corr)))))
  }
  resij[c(7,8),6] <- c(0,1)
  
  glmnet_res <- glmnet::cv.glmnet(x=z,y=y,family = family,alpha=0)
  coefs[,2] <- coef(glmnet_res)[-1]
  resij[c(7,8),7] <- c(glmnet_res$lambda.1se,glmnet_res$glmnet.fit$dev.ratio[glmnet_res$index["1se",1]])
  
  glmnet_res <- glmnet::glmnet(x=z,y=y,family = family,alpha=0,lambda = c(sqrt(n)+sqrt(p),10,1,0.1,0.01))
  coefs[,3:7] <- as.matrix(glmnet_res$beta)
  resij[c(7,8),8:12] <- rbind(glmnet_res$lambda,glmnet_res$dev.ratio)
  
  glmnet_res <- glmnet::glmnet(x=z,y=y,family = family,alpha=0,lambda.min.ratio = 1e-4 /lam_max)
  for (l in 1:3) {
    mydevrat <- c(0.8,0.95,0.999)[l]
    lam_ind <- which.min(glmnet_res$lambda[glmnet_res$dev.ratio<=mydevrat])
    coefs[,7+l] <- coef(glmnet_res,s=glmnet_res$lambda[lam_ind])[-1]
    resij[c(7,8),12+l] <- c(glmnet_res$lambda[lam_ind],glmnet_res$dev.ratio[lam_ind])
  }
  
  coefs[,11] <- apply(z,2,function(zj){
    glm_res <- glm(y~zj,family=family,start=c(1,0))
    glm_res$coefficients[2]
  })
  coefs[,12] <- as.numeric(beta*xscale)
  
  # corMat <- apply(coefs,2,function(col1){
  #   apply(coefs,2,function(col2) cor(col1,col2))
  # })
  # corMat
  # corMat2 <- apply(coefs[data$ind,],2,function(col1){
  #   cor(col1,coefs[data$ind,"True_Beta"])
  # })
  # corMat2
  # resij[c(7,8),13:15]
  
  RP_CW <- SPAR:::generate_cw_rp(m=m,p=p,coef = sample(c(-1,1),p,replace=TRUE))
  reductions <- list("Gaussian"=SPAR::generate_gaussian_rp(m,p),
                     "Sparse"=SPAR::generate_sparse_rp(m,p,psi=1/6),
                     "SparseCW"=RP_CW
                     )
  
  resij[-c(7,8),3:5] <- sapply(reductions, function(Phi){
    zz <- tcrossprod(z,Phi)
    zztest <- tcrossprod(ztest,Phi)
    glmnet_res <- glmnet::glmnet(zz,y,family = family,alpha=0)
    mar_coef <- coef(glmnet_res,s=min(glmnet_res$lambda))
    eta_hat <- as.numeric(zztest%*%mar_coef[-1] + mar_coef[1])
    if (family$family=="binomial") {
      AUC <- my_AUC(data$ytest,family$linkinv(eta_hat))
    } else {
      AUC <- NA
    }
    c(AUC,mean((data$ytest-family$linkinv(eta_hat))^2),mean((eta_hat-xtest%*%beta - alpha)^2),NA,NA,NA,NA,NA)
  })
  
  resij[-c(7,8),-c(1:5)] <- apply(coefs,2, function(coef){
    Phi <- RP_CW
    Phi@x <- coef
    zz <- tcrossprod(z,Phi)
    zztest <- tcrossprod(ztest,Phi)
    if (max(abs(coef))<1e-12) {
      mar_coef <- c(coef(glm(y~1,family = family)),rep(0,m))
    } else {
      glmnet_res <- glmnet::glmnet(zz,y,family = family,alpha=0)
      mar_coef <- coef(glmnet_res,s=min(glmnet_res$lambda))
    }
    eta_hat <- as.numeric(zztest%*%mar_coef[-1] + mar_coef[1])
    if (family$family=="binomial") {
      AUC <- my_AUC(data$ytest,family$linkinv(eta_hat))
    } else {
      AUC <- NA
    }
    c(AUC,mean((data$ytest-family$linkinv(eta_hat))^2),mean((eta_hat-xtest%*%beta - alpha)^2),
      performance(prediction(as.numeric(abs(mar_coef[-1]/xscale)),as.numeric(beta!=0)),measure="auc",fpr.stop=n/2/(p-a))@y.values[[1]]*2*(p-a)/n,
      cor(beta[data$ind],(coef/xscale)[data$ind]),
      performance(prediction(as.numeric(abs(coef/xscale)),as.numeric(beta!=0)),measure="auc",fpr.stop=n/2/(p-a))@y.values[[1]]*2*(p-a)/n,
      mean(abs(coef[data$ind]))/mean(abs(coef[-data$ind])),
      mean(data$ind %in% order(abs(coef),decreasing = TRUE)[1:(3*a)])
    )
  })
  
  # ElNet 
  tmp_res <- myElNet(x,y,xtest,family=family)
  eta_hat <- xtest%*%tmp_res$beta + tmp_res$intercept
  if (family$family=="binomial") {
    AUC <- my_AUC(data$ytest,tmp_res$yhat)
  } else {
    AUC <- NA
  }
  resij[-c(7,8),1] <- c(AUC,mean((data$ytest-tmp_res$yhat)^2),mean((eta_hat-xtest%*%beta - alpha)^2),
                 performance(prediction(as.numeric(abs(tmp_res$beta)),as.numeric(beta!=0)),measure="auc",fpr.stop=n/2/(p-a))@y.values[[1]]*2*(p-a)/n,
                 cor(beta[data$ind],tmp_res$beta[data$ind]),
                 performance(prediction(as.numeric(abs(tmp_res$beta)),as.numeric(beta!=0)),measure="auc",fpr.stop=n/2/(p-a))@y.values[[1]]*2*(p-a)/n,
                 mean(abs(tmp_res$beta[data$ind])) / mean(abs(tmp_res$beta[-data$ind])),
                 mean(data$ind %in% order(abs(tmp_res$beta),decreasing = TRUE)[1:(3*a)])
                 )
  
  # AdLASSO
  tmp_res <- myAdLASSO(x,y,xtest,family=family)
  eta_hat <- xtest%*%tmp_res$beta + tmp_res$intercept
  if (family$family=="binomial") {
    AUC <- my_AUC(data$ytest,tmp_res$yhat)
  } else {
    AUC <- NA
  }
  resij[-c(7,8),2] <- c(AUC,mean((data$ytest-tmp_res$yhat)^2),mean((eta_hat-xtest%*%beta - alpha)^2),
                        performance(prediction(as.numeric(abs(tmp_res$beta)),as.numeric(beta!=0)),measure="auc",fpr.stop=n/2/(p-a))@y.values[[1]]*2*(p-a)/n,
                        cor(beta[data$ind],tmp_res$beta[data$ind]),
                        performance(prediction(as.numeric(abs(tmp_res$beta)),as.numeric(beta!=0)),measure="auc",fpr.stop=n/2/(p-a))@y.values[[1]]*2*(p-a)/n,
                        mean(abs(tmp_res$beta[data$ind])) / mean(abs(tmp_res$beta[-data$ind])),
                        mean(data$ind %in% order(abs(tmp_res$beta),decreasing = TRUE)[1:(3*a)])
                        )

  # res[i,,,j] <- resij
  # res[i,,,j]
  cat(sprintf('Finished rep %d / %d for setting %d / %d at %s.\n',i,nrep,j,nset,Sys.time()))
  resij
}

for (j in 1:nset) {
  for (i in 1:nrep) {
    res[i,,,j] <- parres[[j]][[i]]
  }
}


saveRDS(res,"../saved_results/Result_ScrRP_Experiment_all.rds")
