
pacman::p_load(foreach,parallel,tidyr,dplyr,ROCR,SPAR,e1071)
source("../functions/glm_data_generation.R")
source("../functions/glm_methods.R")
source("../functions/multi_assign.R")


# compare prediction performance (AUC) and link estimation (MSLE) of different random projections
simulation_settings <- tibble(n=200, p=c(2000),ntest=1000, act_setting="medium", cov_setting=c("group"),
                              signal_strength=c(5,1/16,50,500,1/8), avg_exp = c(1,10,0.5,0.7,10),
                              family=list(gaussian(identity),gaussian(log),binomial(logit),binomial(cloglog),poisson(log)),
                              weight_binom = c(1,1,1,1,1))

simulation_settings <- rbind(simulation_settings,
                             tibble(n=200, p=c(2000),ntest=1000, act_setting="sparse", cov_setting=c("group"),
                              signal_strength=c(5,1/16,50,500,1/8), avg_exp = c(1,10,0.5,0.7,10),
                              family=list(gaussian(identity),gaussian(log),binomial(logit),binomial(cloglog),poisson(log)),
                              weight_binom = c(1,1,1,1,1)))

simulation_settings <- rbind(simulation_settings,
                             tibble(n=200, p=c(2000),ntest=1000, act_setting="dense", cov_setting=c("group"),
                              signal_strength=c(5,1/16,50,500,1/8), avg_exp = c(1,10,0.5,0.7,10),
                              family=list(gaussian(identity),gaussian(log),binomial(logit),binomial(cloglog),poisson(log)),
                              weight_binom = c(1,1,1,1,1)))


simulation_settings <- simulation_settings %>% mutate(a = round(ifelse(act_setting=="sparse",2*log(p),
                                                                       ifelse(act_setting=="medium",n/2+2*log(p),p/4))))

nset <- nrow(simulation_settings)
methods <- list("Ensemble_mlogp"=function(x,y,xtest,family){mySPAR(x,y,xtest,family=family,nummods = c(20), nlambda = 1,
                                                                   nscreen = nrow(x)*2,type.rpm="cwdatadriven",
                                                                   mslow = round(log(ncol(x))), msup = 1+round(log(ncol(x))))},
                "Ensemble_m3logp"=function(x,y,xtest,family){mySPAR(x,y,xtest,family=family,nummods = c(20), nlambda = 1,
                                                                     nscreen = nrow(x)*2,type.rpm="cwdatadriven",
                                                                     mslow = round(3*log(ncol(x))), msup = 1+round(3*log(ncol(x))))},
                "Ensemble_mn4"=function(x,y,xtest,family){mySPAR(x,y,xtest,family=family,nummods = c(20), nlambda = 1,
                                                                     nscreen = nrow(x)*2,type.rpm="cwdatadriven",
                                                                     mslow = round(nrow(x)/4), msup = 1+round(nrow(x)/4))},
                "Ensemble_mn2"=function(x,y,xtest,family){mySPAR(x,y,xtest,family=family,nummods = c(20), nlambda = 1,
                                                                     nscreen = nrow(x)*2,type.rpm="cwdatadriven",
                                                                     mslow = round(nrow(x)/2), msup = 1+round(nrow(x)/2))},
                "Ensemble_m3n4"=function(x,y,xtest,family){mySPAR(x,y,xtest,family=family,nummods = c(20), nlambda = 1,
                                                                     nscreen = nrow(x)*2,type.rpm="cwdatadriven",
                                                                     mslow = round(3*nrow(x)/4), msup = 1+round(3*nrow(x)/4))},
                "Ensemble_mStoch"=function(x,y,xtest,family){mySPAR(x,y,xtest,family=family,nummods = c(20), nlambda = 1,
                                                                     nscreen = nrow(x)*2,type.rpm="cwdatadriven",
                                                                     mslow = round(log(ncol(x))), msup = round(nrow(x)/2))},
                "Ensemble_mStochHigh"=function(x,y,xtest,family){mySPAR(x,y,xtest,family=family,nummods = c(20), nlambda = 1,
                                                                     nscreen = nrow(x)*2,type.rpm="cwdatadriven",
                                                                     mslow = round(2*log(ncol(x))), msup = round(3*nrow(x)/4))}
                )

measures <- c("AUC","bAcc","Acc","rMSPE","rDev","rDev_tr","rMSLE","pAUC","Precision","Recall","Sign_ratio_Scr","Cor_Scr","NumAct","Time")

nrep <- 100

nmethods <- length(methods)
nset <- nrow(simulation_settings)
nmeas <- length(measures)
res <- array(c(0),dim=c(nrep,nmeas,nmethods,nset),dimnames = list(reps=NULL,measures=measures,
                                                                  method=names(methods),settings=paste("Scenario",1:nset)))
attributes(res)$settings <- simulation_settings


unlink("../saved_results/log.txt")
my.cluster <- parallel::makeCluster(7, type = "PSOCK", outfile = "../saved_results/log.txt")
doParallel::registerDoParallel(cl = my.cluster)
# foreach::getDoParRegistered()
clusterExport(my.cluster,c('simulation_settings','nmeas','nmethods','methods'), envir = environment())
clusterEvalQ(my.cluster, {  
  source("../functions/glm_data_generation.R")
  source("../functions/glm_methods.R")
  source("../functions/multi_assign.R")
})

# i <- j <- 3
foreach(j = 1:nset) %:%
  foreach(i=1:nrep) %dopar% {
    parresi <- matrix(c(0),nmeas,nmethods)
    c(n,p,ntest,act_setting,cov_setting,signal_strength,avg_exp,family,a) %<-% simulation_settings[j,]
    family <- family[[1]]
    set.seed((1234+i)^2)
    data <- generate_data_glm(n,p,cov_setting,ntest,signal_strength=signal_strength,avg_exp = avg_exp,a=a,family=family)
    x <- data$x
    y <- data$y
    xtest <- data$xtest
    ytest <- data$ytest
    
    rDev_const <- sum(family$dev.resids(ytest,rep(mean(y),ntest),1))
    rDev_const_tr <- sum(family$dev.resids(y,rep(mean(y),n),1))
    
    rMSPE_const <- mean((ytest-mean(y))^2)
    rMSLE_const <- mean((xtest%*%data$beta + data$alpha - family$linkfun(mean(y)))^2)
    
    for (k in 1:nmethods) {
      set.seed((1234+i)^2 + k)
      tstamp <- Sys.time()
      callres <- tryCatch( methods[[k]](x,y,xtest,family),
                           error=function(error_message) {
                             message(paste0("Error in ",names(methods)[k],": ",error_message))
                             return(NULL)
                           })
      tmp_time <- as.numeric(Sys.time() - tstamp,units="secs")
      
      if (is.null(callres)) {
        parresi[,k] <- c(rep(NA,nmeas-1), # numAct
                         tmp_time) 
      } else {
        if (!is.null(callres$scr_coef)) {
          scr_sign <- mean(sign(data$beta[data$ind])==sign(callres$scr_coef[data$ind]))
          scr_cor <- cor(as.numeric(data$beta[data$ind]),callres$scr_coef[data$ind])
        } else {
          scr_sign <- scr_cor <- NA
        }
        if ("dgCMatrix" %in% class(callres$beta)) {
          tmp_ind <- callres$beta@i+1
        } else {
          tmp_ind <- which(callres$beta!=0)
        }
        tr_ind <- est_ind <- numeric(p)
        tr_ind[data$ind] <- 1
        est_ind[tmp_ind] <- 1
        if (is.null(callres$beta) | is.null(callres$intercept) ) {
          rMSLE <- NA
        } else {
          rMSLE <- mean((xtest%*%(callres$beta-data$beta) - data$alpha + callres$intercept)^2) / rMSLE_const
        }
        if (family$family=="binomial") {
          tmpAUC <- my_AUC(ytest,callres$yhat)
          tmpbAcc <- (mean(callres$yhat[ytest==1]>0.5) + mean(callres$yhat[ytest==0]<0.5))/2
        } else {
          tmpAUC <- tmpbAcc <- NA
        }
        parresi[,k] <- c(tmpAUC,
                         tmpbAcc,
                         mean(ytest==round(callres$yhat)),
                         mean((callres$yhat-ytest)^2)/rMSPE_const, # rMSPE
                         sum(family$dev.resids(ytest,callres$yhat,1)) / rDev_const, 
                         sum(family$dev.resids(y,callres$yhat_tr,1)) /rDev_const_tr,
                         rMSLE,
                         performance(prediction(as.numeric(abs(callres$beta)),as.numeric(data$beta!=0)),measure="auc",fpr.stop=n/2/(p-a))@y.values[[1]]*2*(p-a)/n,
                         ifelse(length(tmp_ind)==0,0,mean(tmp_ind %in% data$ind)), # precision
                         mean(data$ind %in% tmp_ind), # recall
                         scr_sign,
                         scr_cor,
                         length(tmp_ind), # numAct
                         tmp_time)
        
      }
      # res[i,,k,j] <- parresi[,k]
      # res[i,,k,j]
      # message(sprintf("Finished method %d / %d!",k,nmethods))
    }
    # res[i,,,j] <- parresi
    # res[i,,,j]
    cat(sprintf('Finished rep %d / %d for setting %d / %d at %s.\n',i,nrep,j,nset,Sys.time()))
    warnings()
    saveRDS(parresi,paste0("../tmp_saved_files/GLM_Sims_M_i",i,"j",j,".rds"))
    parresi
  }

for (j in 1:nset) {
  for (i in 1:nrep) {
    res[i,,,j] <- readRDS(paste0("../tmp_saved_files/GLM_Sims_M_i",i,"j",j,".rds"))
  }
}

saveRDS(res, 
        file = sprintf("../saved_results/SPARglm_M_Ensemble_screen_nset%d_reps%d_nmeth%d.rds",nset,nrep,nmethods))
parallel::stopCluster(cl = my.cluster)

warnings()
