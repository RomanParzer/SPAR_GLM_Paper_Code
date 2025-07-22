
pacman::p_load(foreach,parallel,tidyr,dplyr,ROCR,SPAR,e1071)
source("../functions/glm_data_generation.R")
source("../functions/glm_methods.R")
source("../functions/multi_assign.R")

logit_link <- mybinom_link <- make.link("logit")
mybinom_link$linkfun <- function(mu) logit_link$linkfun((mu-0.5)*0.99998 + 0.5)
mybinom_link$linkinv <- function(eta) (logit_link$linkinv(eta) - 0.5) / 0.99998 + 0.5
mybinom_link$mu.eta <- function(eta) logit_link$mu.eta(eta)/ 0.99998
mybinom_link$valideta <- function(eta) all(abs(eta)<=logit_link$linkfun(0.99999))
mybinom_link$name <- "logit_th" 
mybinom_fam <- binomial(link = mybinom_link)

# compare prediction performance (AUC) and link estimation (MSLE) of different random projections
simulation_settings <- tibble(n=200, p=c(500), act_setting="medium", cov_setting=c("group"),
                              signal_strength=c(10,1/8,100,1000,1/4,1,5), avg_exp = c(1,10,0.5,0.7,10,0.5,0.5),
                              family=list(gaussian(identity),gaussian(log),binomial(logit),binomial(cloglog),poisson(log),binomial(logit),mybinom_fam),
                              weight_binom = c(1,1,1,1,1,100,1),a = round(n/2 + 2*log(p)))

nrep <- 10
methods <- c("L2_dev08","L2_dev095","L2_dev0999",
             "GLM_HOLP","GLM_HOLP_yuncent","GLM_HOLP_Xunstd","GLM_HOLP_Xyunstd",
             "marGLM","true",
             "L2_nointcpt_dev08","L2_nointcpt_dev095","L2_nointcpt_dev0999",
             "manual","manual_noint")
nmethods <- length(methods)
nset <- nrow(simulation_settings)


mseMats <- corMats <- array(c(0),dim=c(nrep,nmethods,nmethods,nset),dimnames = list(reps=NULL,method1=methods,
                                                                                 method2=methods,settings=paste("Scenario",1:nset)))
cor_tr_act <- array(c(0),dim=c(nrep,nmethods,nset),dimnames = list(reps=NULL,method=methods,
                                                                            settings=paste("Scenario",1:nset)))

unlink("../saved_results/log.txt")
my.cluster <- parallel::makeCluster(7, type = "PSOCK", outfile = "../saved_results/log.txt")
doParallel::registerDoParallel(cl = my.cluster)
# foreach::getDoParRegistered()
clusterExport(my.cluster,c('simulation_settings','nmethods','methods','logit_link'), envir = environment())
clusterEvalQ(my.cluster, {  
  source("../functions/glm_data_generation.R")
  source("../functions/glm_methods.R")
  source("../functions/multi_assign.R")
})

# select random seed and setting
# i <- 1
# j <- 6
foreach(j = 1:nset) %:%
  foreach(i=1:nrep) %dopar% {
    c(n,p,act_setting,cov_setting,signal_strength,avg_exp,family2,weight_binom,a) %<-% simulation_settings[j,]
    family2 <- family2[[1]]
    set.seed((1234+i)^2)
    data <- generate_data_glm(n,p,family=family2,cov_setting=cov_setting,weight_binom=weight_binom,
                              a=a,signal_strength = signal_strength, avg_exp = avg_exp)
    
    x <- data$x
    y <- data$y
    beta <- data$beta
    alpha <- data$alpha
    
    if (family2$family=="binomial" & weight_binom>1) {
      y_fit <- cbind(y*weight_binom,(1-y)*weight_binom)
    } else {
      y_fit <- y
    }
    
    check_gy <- !is.nan(sum(family2$linkfun(y))*0)
    
    coefs <- matrix(c(0),p,nmethods)
    colnames(coefs) <- methods 
    
    # testing different scales internally
    x[,1:(p/2)] <- x[,1:(p/2)]*2 + 3
    beta[1:(p/2)] <- beta[1:(p/2)]/2
    alpha <- alpha - sum(3*beta[1:(p/2)])
    
    xcenter <- apply(x,2,mean)
    xscale <- apply(x,2,sd)
    # xcenter <- rep(0,p)
    # xscale <- rep(1,p)
    z <- scale(x,center=xcenter,scale=xscale)
    
    # xcenterrob <- apply(x,2,median)
    # xscalerob <- apply(x,2,mad)
    # zrob <- scale(x,center=xcenterrob,scale=xscalerob)
    
    eig <- eigen(tcrossprod(z),symmetric = TRUE)
    myinv <- tcrossprod(eig$vectors[,eig$values>1e-8]%*%diag(1/sqrt(eig$values[eig$values>1e-8]),nrow = sum(eig$values>1e-8)))
    if (family2$family=="poisson") {
      y_corr <- ifelse(y==0,1e-4,y)
    } else if (family2$family=="binomial") {
      y_corr <- ifelse(y==0,1e-4,ifelse(y==1,1-1e-4,y))
    } else {
      y_corr <- y
    }
    # # tmp diff link function
    family <- family2
    # family$linkfun <- function(mu){
    #   family2$linkfun(mu) - mean(family2$linkfun(y_corr))
    # }
    # family$linkinv <- function(mu){
    #   family2$linkinv(mu) + mean(family2$linkfun(y_corr))
    # }
    
    coefs[,4] <- as.numeric(crossprod(z,myinv%*%(family$linkfun(y_corr) - mean(family$linkfun(y_corr))))) / xscale
    coefs[,5] <- as.numeric(crossprod(z,myinv%*%(family$linkfun(y_corr)))) / xscale
    
    coefs[,6] <- as.numeric(crossprod(x,solve(tcrossprod(x),family$linkfun(y_corr)- mean(family$linkfun(y_corr))))) / 1
    coefs[,7] <- as.numeric(crossprod(x,solve(tcrossprod(x),family$linkfun(y_corr)))) / 1
    
    tmp_sc <- apply(x,2,function(col)sqrt(var(col)*(n-1)/n))
    z2 <- scale(x,center=colMeans(x),scale=tmp_sc)
    lam_max <- 1000 * max(abs(t(y)%*%z2[,tmp_sc>0]))/n*family$mu.eta(family$linkfun(mean(y)))/family$variance(mean(y))
    glmnet_res <- glmnet::glmnet(x=x,y=y_fit,family = family,alpha=0,lambda.min.ratio = 1e-4 /lam_max)
    for (l in 1:3) {
      mydevrat <- c(0.8,0.95,0.999)[l]
      lam_ind <- which.min(glmnet_res$lambda[glmnet_res$dev.ratio<=mydevrat])
      coefs[,l] <- coef(glmnet_res,s=glmnet_res$lambda[lam_ind])[-1]
      # # does lam * beta_lam go to zero? yes it seems
      # print(paste("lam*betalam:",glmnet_res$lambda[lam_ind] * mean(coefs[,l]^2)))
      # # is beta_lam in span of X'? not due to scaling
      # Pbetlam <- crossprod(x,solve(tcrossprod(x),x%*%coefs[,l]))
      # print(paste("proj:",mean((coefs[,l] - Pbetlam )^2)/mean(coefs[,l]^2)))
    }
    
    coefs[,8] <- apply(x,2,function(zj){
      glm_res <- glm(y_fit~zj,family=family,start=c(1,0))
      glm_res$coefficients[2]
    })
    coefs[,9] <- as.numeric(beta)
    
    glmnet_res <- glmnet::glmnet(x=x,y=y_fit,family = family,alpha=0,lambda.min.ratio = 1e-4 /lam_max,intercept = FALSE)
    for (l in 1:3) {
      mydevrat <- c(0.8,0.95,0.999)[l]
      lam_ind <- which.min(glmnet_res$lambda[glmnet_res$dev.ratio<=mydevrat])
      coefs[,9+l] <- coef(glmnet_res,s=glmnet_res$lambda[lam_ind])[-1]
      # print(paste("lam*betalam:",glmnet_res$lambda[lam_ind] * mean(coefs[,9+l]^2)))
      # Pbetlam <- crossprod(x,solve(tcrossprod(x),x%*%coefs[,9+l]))
      # print(paste("proj:",mean((coefs[,9+l] - Pbetlam )^2)/mean(coefs[,9+l]^2)))
    }
    
    myminObj <- function(mybeta,lam = 1e-3) {
      sum(family$dev.resids(y = y,mu = family$linkinv(mybeta[1] + z%*%mybeta[-1]),wt=1)) / 2 + lam/2 * sum(mybeta[-1]^2) 
    }
    
    lam <- 1e-3
    opt_res <- optim(c(family$linkfun(mean(y)),as.numeric(rep(0,p))),myminObj,lam=lam,method = "BFGS",
                     control = list(maxit=700))
    # str(opt_res)
    coefs[,13] <- opt_res$par[-1] / xscale
    # print(paste("lam*betalam:",lam * mean(coefs[,13]^2)))
    # Pbetlam <- crossprod(x,solve(tcrossprod(x),x%*%coefs[,13]))
    # print(paste("proj:",mean((coefs[,13] - Pbetlam )^2)/mean(coefs[,13]^2))) # not for z
    
    myminObj_noint <- function(mybeta,lam = 1e-3) {
      sum(family$dev.resids(y = y,mu = family$linkinv(x%*%mybeta),wt=1)) / 2 + lam/2 * sum(mybeta^2) 
    }
    lam <- 1e-4
    opt_noint <- optim(as.numeric(rep(0,p)),myminObj_noint,lam=lam,method = "BFGS",
                       control = list(maxit=700))
    # str(opt_noint)
    coefs[,14] <- opt_noint$par 
    # print(paste("lam*betalam:",lam * mean(coefs[,14]^2)))
    # Pbetlam <- crossprod(x,solve(tcrossprod(x),x%*%coefs[,14]))
    # print(paste("proj:",mean((coefs[,14] - Pbetlam )^2)/mean(coefs[,14]^2))) # yes
  
    corMat <- apply(coefs,2,function(col1){
      apply(coefs,2,function(col2) cor(col1,col2))
    })
    # corMat
    mseMat <- apply(coefs,2,function(col1){
      apply(coefs,2,function(col2) mean((col1-col2)^2)/mean(col1^2))
    })
    # mseMat
    corMat2 <- apply(coefs[data$ind,],2,function(col1){
      cor(col1,coefs[data$ind,"true"])
    })
    # corMat2
    resij <- list(cor = corMat, mse = mseMat, cor_tr = corMat2)
    cat(sprintf('Finished rep %d / %d for setting %d / %d at %s.\n',i,nrep,j,nset,Sys.time()))
    warnings()
    saveRDS(resij,paste0("../tmp_saved_files/TestGLMHOLPi",i,"j",j,".rds"))
    resij
  }

for (j in 1:nset) {
  for (i in 1:nrep) {
    resij <- readRDS(paste0("../tmp_saved_files/TestGLMHOLPi",i,"j",j,".rds"))
    mseMats[i,,,j] <- resij$mse
    corMats[i,,,j] <- resij$cor
    cor_tr_act[i,,j] <- resij$cor_tr
  }
}

saveRDS(list(mseMats=mseMats,
             corMats=corMats,
             cor_tr_act = cor_tr_act), 
        file = sprintf("../saved_results/TestGLM_HOLP_nset%d_reps%d_nmeth%d.rds",nset,nrep,nmethods))
parallel::stopCluster(cl = my.cluster)

warnings()


# # # # # # # # # # # # look at results here
# 
res <- readRDS("../saved_results/TestGLM_HOLP_nset7_reps10_nmeth14.rds")
# important part: manual close to GLM_HOLP, manual_noint close to GLM_HOLP_Xyunstd
# for canonical well defined cases (gy_i in R) j = 1,5,6,7 (j=5 poisson log numerically dificult for manual_oint, need to check convergence)
j <- 7
simulation_settings$family[[j]]
apply(res$mseMats[,,c(4,7,13,14),j],c(2,3),mean)
apply(res$corMats[,,c(4,7,13,14),j],c(2,3),mean)
# apply(res$corMats[,,,j],c(2,3),sd)
apply(res$cor_tr_act[,,j],2,mean)
