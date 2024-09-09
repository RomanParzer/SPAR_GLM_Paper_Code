#####################################################
### Testing sparse projected averaged regression SPAR
#####################################################

pacman::p_load(foreach,parallel,tidyr,dplyr)
source("../functions/glm_data_generation.R")
source("../functions/glm_methods.R")
source("../functions/multi_assign.R")

################ adapt #########################################################################

nrep <- 100

methods <- list("Ridge"=function(x,y,xtest,family){myElNet(x,y,xtest,alpha=0,family=family)},
                "LASSO"=function(x,y,xtest,family){myElNet(x,y,xtest,alpha=1,family=family)},
                "AdLASSO"=myAdLASSO,
                "ElNet"=myElNet,
                "SIS"=function(x,y,xtest,family){mySIS(x,y,xtest,family=family,iter=FALSE)},
                "SVM"=mySVM,
                "RF"=myRF,
                "RP_CW"=function(x,y,xtest,family){mySPAR(x,y,xtest,family=family,nummods = c(1), nlambda = 1,nscreen = ncol(x),type.rpm="cw")},
                "RP_CW_Ensemble"=function(x,y,xtest,family){mySPAR(x,y,xtest,family=family,nummods = c(50), nlambda = 1,nscreen = ncol(x),type.rpm="cw")},
                "ourRP"=function(x,y,xtest,family){mySPAR(x,y,xtest,family=family,nummods = c(1), nlambda = 1,nscreen = ncol(x),
                                                          type.rpm="cwdatadriven",type.screening = "ridge")},
                "ourRP_Ensemble"=function(x,y,xtest,family){mySPAR(x,y,xtest,family=family,nummods = c(50), nlambda = 1,nscreen = ncol(x),
                                                                   type.rpm="cwdatadriven",type.screening = "ridge")},
                "TARP"=function(x,y,xtest,family){mySPAR(x,y,xtest,family=family,nummods = c(20), nlambda = 1,
                                                         type.rpm="sparse",type.screening = "marglik")},
                "TARP+"=function(x,y,xtest,family){mySPAR(x,y,xtest,family=family,nummods = c(20), nlambda = 1,
                                                          type.rpm="cwdatadriven",type.screening = "marglik")},
                "SPARconv"=function(x,y,xtest,family){mySPAR(x,y,xtest,family=family,nummods = c(20), nlambda = 1,
                                                             type.rpm="sparse",type.screening = "ridge")},
                "SPAR"=function(x,y,xtest,family){mySPAR(x,y,xtest,family=family,nummods = c(20), nlambda = 20)},
                "SPAR split"=function(x,y,xtest,family){mySPAR(x,y,xtest,family=family,nummods = c(20), nlambda = 20, split_data=TRUE)},
                "SPAR res avg"=function(x,y,xtest,family){mySPAR(x,y,xtest,family=family,nummods = c(20), nlambda = 20, avg_type="response")},
                "SPAR CV"=function(x,y,xtest,family){mySPAR(x,y,xtest,family=family,nummods=c(10,20,30,50),doCV=TRUE)},
                "SPAR CV 1se"=function(x,y,xtest,family){mySPAR(x,y,xtest,family=family,nummods=c(10,20,30,50),opt_par = "1se",doCV=TRUE)})

measures <- c("AUC","bAcc","Acc","rMSPE","rDev","rDev_tr","rMSLE","pAUC","Precision","Recall","Sign_ratio_Scr","Cor_Scr","NumAct","Time")

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
                             tibble(n=200, p=c(2000),ntest=200, act_setting="medium", cov_setting=c("ar1"),
                                    signal_strength=c(10,1/8,100,1000,1/4), avg_exp = c(1,10,0.5,0.7,10),
                                    family=list(gaussian(identity),gaussian(log),binomial(logit),binomial(cloglog),poisson(log))))

simulation_settings <- rbind(simulation_settings,
                             tibble(n=200, p=c(2000),ntest=200, act_setting="medium", cov_setting=c("comsym"),
                                    signal_strength=c(10,1/8,100,1000,1/4), avg_exp = c(1,10,0.5,0.7,10),
                                    family=list(gaussian(identity),gaussian(log),binomial(logit),binomial(cloglog),poisson(log))))

simulation_settings <- rbind(simulation_settings,
                             tibble(n=200, p=c(2000),ntest=200, act_setting="medium", cov_setting=c("ind"),
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

################## simulations ##############################################################################

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


# i <- j <- 1

parres <- foreach(j = 1:nset) %:%
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
    
    rDev_const <- sum(family$dev.resids(ytest,rep(family$linkinv(coef(glm(y~1,family = family,start=1))),ntest),1))
    rDev_const_tr <- sum(family$dev.resids(y,rep(family$linkinv(coef(glm(y~1,family = family,start=1))),n),1))
    
    rMSPE_const <- mean((ytest-mean(y))^2)
    rMSLE_const <- mean((xtest%*%data$beta + data$alpha)^2)
    
    # k <- 15
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
    parresi
  }

for (j in 1:nset) {
  for (i in 1:nrep) {
    res[i,,,j] <- parres[[j]][[i]]
  }
}

saveRDS(res, 
        file = sprintf("../saved_results/SPARglm_sims_nset%d_reps%d_nmeth%d.rds",nset,nrep,nmethods))
parallel::stopCluster(cl = my.cluster)

warnings()
