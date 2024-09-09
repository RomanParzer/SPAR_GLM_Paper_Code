
pacman::p_load(foreach,parallel,tidyr,dplyr,ROCR,SPAR)
source("../functions/glm_data_generation.R")
source("../functions/glm_methods.R")
source("../functions/multi_assign.R")


# compare prediction performance (Acc) and link estimation (MSLE) of different random projections
simulation_settings <- tibble(n=200, p=c(2000),ntest=200, act_setting=c("sparse","medium","dense"), cov_setting=c("group"),
                              signal_strength=1/8, avg_exp = 10,family=list(gaussian(log)))
simulation_settings <- simulation_settings %>% mutate(a = round(ifelse(act_setting=="sparse",2*log(p),
                                                                       ifelse(act_setting=="medium",n/2+2*log(p),p/4))))

nset <- nrow(simulation_settings)
methods <- c("ElNet","AdLASSO","Gaussian","Sparse","SparseCW",
             "GLM_HOLP_unstd","GLM_HOLP","GLM_HOLP_ycent","GLM_HOLP_sqrtnp",
             "Ridge_cv","Ridge_def_min","Ridge_small","Ridge_no_intcpt",
             "marGLM","True_Beta")

nmethods <- length(methods)
nrep <- 10^2
res <- array(c(0),dim=c(nrep,6,nmethods,nset),
             dimnames = list(rep=1:nrep,
                             measure=c("MSPE","MSLE","Cor","pAUC","Cor_scr","pAUC_scr"),
                             method=methods,
                             setting=1:nset))
attributes(res)$settings <- simulation_settings

unlink("../saved_results/log.txt")
my.cluster <- parallel::makeCluster(4, type = "PSOCK", outfile = "../saved_results/log.txt")
doParallel::registerDoParallel(cl = my.cluster)
# foreach::getDoParRegistered()
clusterExport(my.cluster,c('simulation_settings','nmethods','methods'), envir = environment())
clusterEvalQ(my.cluster, {  
  pacman::p_load(ROCR,SPAR)
  source("../functions/glm_data_generation.R")
  source("../functions/glm_methods.R")
  source("../functions/multi_assign.R")
})

# i <- j <- 2

parres <- foreach(j = 1:nset) %:%
  foreach(i=1:nrep) %dopar% {
    resij <- matrix(c(0),6,nmethods)
    c(n,p,ntest,act_setting,cov_setting,signal_strength,avg_exp,family,a) %<-% simulation_settings[j,]
    family <- family[[1]]
    set.seed((1234+i)^2)
    data <- generate_data_glm(n,p,family=family,cov_setting=cov_setting,ntest=ntest,a=a,signal_strength = signal_strength, avg_exp = avg_exp)
    summary(data$y)
    
    m <- n/4
    x <- data$x
    y <- data$y
    beta <- data$beta
    alpha <- data$alpha
    xtest <- data$xtest
    
    # testing different scales internally
    # x[,1:(p/2)] <- x[,1:(p/2)]*2 + 3
    # xtest[,1:(p/2)] <- xtest[,1:(p/2)]*2 + 3
    # beta[1:(p/2)] <- beta[1:(p/2)]/2
    # alpha <- alpha - sum(3*beta[1:(p/2)])
    
    xcenter <- apply(x,2,median)
    xscale <- apply(x,2,mad)
    
    z <- scale(x,center=xcenter,scale=xscale)
    ztest <- scale(xtest,center = xcenter, scale = xscale)
    
    coefs <- matrix(c(0),p,nmethods-5)
    colnames(coefs) <- methods[-c(1:5)]
    
    y_corr <- ifelse(y<=0,1e-4,y)
    coefs[,1] <- as.numeric(crossprod(x,solve(tcrossprod(x),family$linkfun(y_corr))))*xscale
    coefs[,2] <- as.numeric(crossprod(z,solve(tcrossprod(z),family$linkfun(y_corr))))
    coefs[,3] <- as.numeric(crossprod(z,solve(tcrossprod(z),family$linkfun(y_corr)-mean(family$linkfun(y_corr)))))
    coefs[,4] <- as.numeric(crossprod(z,solve(tcrossprod(z)+diag(n)*(sqrt(n)+sqrt(p)),family$linkfun(y_corr))))
    coefs[,5] <- coef(glmnet::cv.glmnet(x=z,y=y_corr,family = family,alpha=0,start=rep(0,p+1)))[-1]
    
    glmnet_res <- glmnet::glmnet(x=z,y=y_corr,family = family,alpha=0)
    lam_def_min <- min(glmnet_res$lambda)
    coefs[,6] <- coef(glmnet_res,s=lam_def_min)[-1]
    
    glmnet_res <- glmnet::glmnet(x=z,y=y_corr,family = family,alpha=0,lambda.min.ratio = 1e-7/n /(100*lam_def_min))
    lam <- min(glmnet_res$lambda)
    coefs[,7] <- coef(glmnet_res,s=lam)[-1]
    
    glmnet_res <- glmnet::glmnet(x=z,y=y_corr,family = family,alpha=0,intercept = FALSE,lambda=lam)
    lam <- min(glmnet_res$lambda)
    coefs[,8] <- coef(glmnet_res,s=lam)[-1]
    coefs[,9] <- apply(z,2,function(zj){
      glm_res <- glm(y_corr~zj,family=family,start=c(1,0))
      glm_res$coefficients[2]
    })
    coefs[,10] <- as.numeric(beta*xscale)
    
    # corMat <- apply(coefs,2,function(col1){
    #   apply(coefs,2,function(col2) cor(col1,col2))
    # })
    # corMat
    # corMat2 <- apply(coefs[data$ind,],2,function(col1){
    #   apply(coefs[data$ind,],2,function(col2) cor(col1,col2))
    # })
    # corMat2

    RP_CW <- SPAR:::generate_cw_rp(m=m,p=p,coef = sample(c(-1,1),p,replace=TRUE))
    reductions <- list("Gaussian"=SPAR::generate_gaussian_rp(m,p),
                       "Sparse"=SPAR::generate_sparse_rp(m,p,psi=1/6),
                       "SparseCW"=RP_CW
    )
    
    resij[,3:5] <- sapply(reductions, function(Phi){
      zz <- tcrossprod(z,Phi)
      zztest <- tcrossprod(ztest,Phi)
      
      glmnet_res <- glmnet::glmnet(zz,y_corr,family = family,alpha=0)
      mar_coef <- coef(glmnet_res,s=min(glmnet_res$lambda))
      eta_hat <- as.numeric(zztest%*%mar_coef[-1] + mar_coef[1])
      c(mean((data$ytest-family$linkinv(eta_hat))^2),mean((eta_hat-xtest%*%beta - alpha)^2),NA,NA,NA,NA)
    })
    
    resij[,-c(1:5)] <- apply(coefs,2, function(coef){
      Phi <- RP_CW
      Phi@x <- coef
      zz <- tcrossprod(z,Phi)
      zztest <- tcrossprod(ztest,Phi)
      if (max(abs(coef))<1e-12) {
        mar_coef <- c(coef(glm(y_corr~1,family = family)),rep(0,m))
      } else {
        glmnet_res <- glmnet::glmnet(zz,y_corr,family = family,alpha=0)
        mar_coef <- coef(glmnet_res,s=min(glmnet_res$lambda))
      }
      eta_hat <- as.numeric(zztest%*%mar_coef[-1] + mar_coef[1])
      c(mean((data$ytest-family$linkinv(eta_hat))^2),mean((eta_hat-xtest%*%beta - alpha)^2),
        cor(as.numeric(beta),mar_coef[-1]/xscale),
        performance(prediction(as.numeric(abs(mar_coef[-1]/xscale)),as.numeric(beta!=0)),measure="auc",fpr.stop=n/2/(p-a))@y.values[[1]]*2*(p-a)/n,
        cor(beta[data$ind],(coef/xscale)[data$ind]),
        performance(prediction(as.numeric(abs(coef/xscale)),as.numeric(beta!=0)),measure="auc",fpr.stop=n/2/(p-a))@y.values[[1]]*2*(p-a)/n)
    })
    
    # ElNet 
    tmp_res <- myElNet(x,y_corr,xtest,family=family)
    eta_hat <- xtest%*%tmp_res$beta + tmp_res$intercept
    resij[,1] <- c(mean((data$ytest-tmp_res$yhat)^2),mean((eta_hat-xtest%*%beta - alpha)^2),
                   cor(as.numeric(beta),as.numeric(tmp_res$beta)),
                   performance(prediction(as.numeric(abs(tmp_res$beta)),as.numeric(beta!=0)),measure="auc",fpr.stop=n/2/(p-a))@y.values[[1]]*2*(p-a)/n,
                   cor(beta[data$ind],tmp_res$beta[data$ind]),
                   performance(prediction(as.numeric(abs(tmp_res$beta)),as.numeric(beta!=0)),measure="auc",fpr.stop=n/2/(p-a))@y.values[[1]]*2*(p-a)/n)
    
    # AdLASSO
    tmp_res <- myAdLASSO(x,y_corr,xtest,family=family)
    eta_hat <- xtest%*%tmp_res$beta + tmp_res$intercept
    resij[,2] <- c(mean((data$ytest-tmp_res$yhat)^2),mean((eta_hat-xtest%*%beta - alpha)^2),
                   cor(as.numeric(beta),as.numeric(tmp_res$beta)),
                   performance(prediction(as.numeric(abs(tmp_res$beta)),as.numeric(beta!=0)),measure="auc",fpr.stop=n/2/(p-a))@y.values[[1]]*2*(p-a)/n,
                   cor(beta[data$ind],tmp_res$beta[data$ind]),
                   performance(prediction(as.numeric(abs(tmp_res$beta)),as.numeric(beta!=0)),measure="auc",fpr.stop=n/2/(p-a))@y.values[[1]]*2*(p-a)/n)
    
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
# saveRDS(res,"../saved_results/Result_RP_Experiment_gauLog.rds")

res <- readRDS("../saved_results/Result_RP_Experiment_gauLog.rds")

methods <- dimnames(res)[[3]]
settings <- attributes(res)$settings
settings$family <- sapply(settings$family,function(obj) paste0(obj$family,"(",obj$link,")"))

mydf_all <- data.frame(pivot_longer(data.frame(res[,1,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="MSPE"),
                       MSLE=pivot_longer(data.frame(res[,2,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="MSLE")$MSLE,
                       Cor=pivot_longer(data.frame(res[,3,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Cor")$Cor,
                       pAUC=pivot_longer(data.frame(res[,4,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="pAUC")$pAUC,
                       Cor_scr=pivot_longer(data.frame(res[,5,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Cor_scr")$Cor_scr,
                       pAUC_scr=pivot_longer(data.frame(res[,6,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="pAUC_scr")$pAUC_scr,
                       settings[1,],
                       setting=1)

for (k in 2:nrow(settings)) {
  mydf_all <- rbind(mydf_all,
                    data.frame(pivot_longer(data.frame(res[,1,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="MSPE"),
                               MSLE=pivot_longer(data.frame(res[,2,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="MSLE")$MSLE,
                               Cor=pivot_longer(data.frame(res[,3,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Cor")$Cor,
                               pAUC=pivot_longer(data.frame(res[,4,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="pAUC")$pAUC,
                               Cor_scr=pivot_longer(data.frame(res[,5,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Cor_scr")$Cor_scr,
                               pAUC_scr=pivot_longer(data.frame(res[,6,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="pAUC_scr")$pAUC_scr,
                               settings[k,],
                               setting=k)
  )
}

mydf_all$Method <- factor(mydf_all$Method,levels = methods)
mydf_all$act_setting <- factor(mydf_all$act_setting,levels = c("sparse","medium","dense"))


tmp_df <- mydf_all %>% filter(MSPE<30,MSLE<0.3) %>%
  pivot_longer(c(MSPE,MSLE),names_to = "Measure",values_to = "value")
tmp_df$Measure <- factor(tmp_df$Measure,levels = c("MSLE","MSPE"))
tmp_df %>%  filter(Method %in% methods, act_setting=="medium") %>%
  ggplot(aes(x=Method,y=value,fill=Method)) +
  geom_boxplot() +
  theme(legend.position = "none") +
  facet_grid(Measure~.,scales = "free_y")+
  labs(y=" ") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# ggsave("../plots/RPexperiment_gauLog_all.pdf", height = 5, width = 8)

tmp_df <- mydf_all %>% pivot_longer(c(Cor_scr,pAUC_scr),names_to = "Measure",values_to = "value")
tmp_df$Measure <- factor(tmp_df$Measure,levels = c("Cor_scr","pAUC_scr"))
tmp_df %>%  filter(Method %in% methods[-c(3,4,5,15)]) %>%
  ggplot(aes(x=Method,y=value,fill=Method)) +
  geom_boxplot() +
  theme(legend.position = "none") +
  ggh4x::facet_grid2(Measure~act_setting, scales = "free_y",independent = "y") +
  # facet_grid(Measure~act_setting,scales = "free_y")+
  labs(y=" ") +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1))
# ggsave("../plots/ScrExperiment_gauLog_all.pdf", height = 5, width = 8)

#
# tmp_df <- mydf_all %>% pivot_longer(c(Cor,pAUC),names_to = "Measure",values_to = "value")
# tmp_df$Measure <- factor(tmp_df$Measure,levels = c("Cor","pAUC"))
# tmp_df %>%  filter(Method %in% methods[-c(3,4,5)]) %>%
#   ggplot(aes(x=Method,y=value,fill=Method)) +
#   geom_boxplot() +
#   theme(legend.position = "none") +
#   ggh4x::facet_grid2(Measure~act_setting, scales = "free_y",independent = "y") +
#   # facet_grid(Measure~act_setting,scales = "free_y")+
#   labs(y=" ") +
#   theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1))
#