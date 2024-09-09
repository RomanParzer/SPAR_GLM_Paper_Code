#####################################################
### Testing sparse projected averaged regression SPAR
#####################################################

pacman::p_load(foreach,parallel,tidyr,dplyr,readxl)
source("../functions/glm_methods.R")
source("../functions/multi_assign.R")

# add different versions of data?

################ adapt #########################################################################

nrep <- 100

methods <- list("Ridge"=function(x,y,xtest,family){myElNet(x,y,xtest,alpha=0,family=family)},
                "Ridge mse"=function(x,y,xtest,family){myElNet(x,y,xtest,alpha=0,family=family,type.measure = "mse")},
                "LASSO"=function(x,y,xtest,family){myElNet(x,y,xtest,alpha=1,family=family)},
                "AdLASSO"=myAdLASSO,
                "AdLASSO mse"=function(x,y,xtest,family){myAdLASSO(x,y,xtest,family,type.measure = "mse")},
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
                "SPAR"=function(x,y,xtest,family){mySPAR(x,y,xtest,family=family,nummods = c(20), nlambda = 20)},
                "SPAR split"=function(x,y,xtest,family){mySPAR(x,y,xtest,family=family,nummods = c(20), nlambda = 20, split_data=TRUE)},
                "SPAR res avg"=function(x,y,xtest,family){mySPAR(x,y,xtest,family=family,nummods = c(20), nlambda = 20, avg_type="response")},
                "SPAR CV"=function(x,y,xtest,family){mySPAR(x,y,xtest,family=family,nummods=c(10,20,30,50),doCV=TRUE)},
                "SPAR CV 1se"=function(x,y,xtest,family){mySPAR(x,y,xtest,family=family,nummods=c(10,20,30,50),opt_par = "1se",doCV=TRUE)},
                "SPAR CV mse"=function(x,y,xtest,family){mySPAR(x,y,xtest,family=family,nummods=c(10,20,30,50),type.measure="mse",doCV=TRUE)}
)

measures <- c("rMSPE","rDev","rDev_tr","NumAct","Time")


# duration_spectra <- list(x=as.matrix(rbind(cbind(scale="large",readRDS("../data/TribDataPia/X_large_scale_diff.rds")),
#                                              cbind(scale="small",readRDS("../data/TribDataPia/X_small_scale_diff.rds")))),
#                          y=round(c(unlist(readRDS("../data/TribDataPia/y_large_scale.rds")),unlist(readRDS("../data/TribDataPia/y_small_scale.rds")))))
# rownames(duration_spectra$x) <- 1:34
# names(duration_spectra$y) <- 1:34
# write.csv(duration_spectra$x,"../data/tribology/diff_spectra.csv",row.names = TRUE)
# write.csv(duration_spectra$y,"../data/tribology/duration.csv",row.names = TRUE)

diff_spectra <- read.csv("../data/tribology/diff_spectra.csv",header=TRUE,row.names = 1)
diff_spectra[,1] <- ifelse(diff_spectra[,1]=="large",1,0)
colnames(diff_spectra)[-1] <- sapply(colnames(diff_spectra)[-1],function(name) substring(name,2))
duration_spectra <- list(x=as.matrix(diff_spectra),
                         y=read.csv("../data/tribology/duration.csv")$x)

datasets <- list(duration_spectra=duration_spectra)
dataset_sizes <- tibble(dataset=names(datasets),
                        n = sapply(datasets,function(arg) nrow(arg$x) ),
                        p = sapply(datasets,function(arg) ncol(arg$x) ))

################## simulations #############################################################################

nmethods <- length(methods)
nset <- length(datasets)
nmeas <- length(measures)
res <- array(c(0),dim=c(nrep,nmeas,nmethods,3,nset),dimnames = list(reps=NULL,measures=measures,
                                                                  method=names(methods),
                                                                  flink=c("poisson(log)","gaussian(log)","quasipoisson(log)"),
                                                                  datasets=names(datasets)))

unlink("../saved_results/log.txt")
my.cluster <- parallel::makeCluster(7, type = "PSOCK", outfile = "../saved_results/log.txt")
doParallel::registerDoParallel(cl = my.cluster)
# foreach::getDoParRegistered()
clusterExport(my.cluster,c('datasets','nmeas','nmethods','methods'), envir = environment())
clusterEvalQ(my.cluster, {  
  source("../functions/multi_assign.R")
  source("../functions/glm_methods.R")
})
# j <- 1
# i <- 1
parres <- foreach(j = 1:(nset)) %:%
foreach(i=1:nrep) %dopar% {
  parresi <- array(c(0),dim=c(nmeas,nmethods,3))
  
    mydata <- datasets[[j]]
    ntot <- nrow(mydata$x)
    ntest <- ntot%/%4
    n <- ntot - ntest
    p <- ncol(mydata$x)

    set.seed((1234+i)^2)
    testind <- sample(1:ntot,ntest,replace=FALSE)
    x <- as.matrix(mydata$x[-testind,])
    xtest <- as.matrix(mydata$x[testind,])
    y <- mydata$y[-testind]
    ytest <- mydata$y[testind]
    
    rMSPE_const <- mean((ytest-mean(y))^2)
    
    family <- poisson(log)
    rDev_const <- sum(family$dev.resids(ytest,rep(family$linkinv(coef(glm(y~1,family = family,start=1))),ntest),1))
    rDev_const_tr <- sum(family$dev.resids(y,rep(family$linkinv(coef(glm(y~1,family = family,start=1))),n),1))
    
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
        parresi[,k,1] <- c(rep(NA,nmeas-1), # numAct
                         tmp_time) 
      } else {
        ypred <- round(callres$yhat)
        parresi[,k,1] <- c(mean((callres$yhat-ytest)^2)/rMSPE_const, # rMSPE
                         sum(family$dev.resids(ytest,callres$yhat,1)) / rDev_const, 
                         sum(family$dev.resids(y,callres$yhat_tr,1)) /rDev_const_tr,
                         sum(callres$beta!=0), # numAct
                         tmp_time)
        
      }
      # res[i,,k,j] <- parresi[,k]
      # res[i,,k,j]
      # message(sprintf("Finished method %d / %d!",k,nmethods))
    }
    family <- gaussian(log)
    rDev_const <- sum(family$dev.resids(ytest,rep(family$linkinv(coef(glm(y~1,family = family,start=1))),ntest),1))
    rDev_const_tr <- sum(family$dev.resids(y,rep(family$linkinv(coef(glm(y~1,family = family,start=1))),n),1))
    for (k in which(!names(methods)%in%c("SVM","RF"))) {
      set.seed((1234+i)^2 + k)
      tstamp <- Sys.time()
      callres <- tryCatch( methods[[k]](x,y,xtest,family),
                           error=function(error_message) {
                             message(paste0("Error in ",names(methods)[k],": ",error_message))
                             return(NULL)
                           })
      tmp_time <- as.numeric(Sys.time() - tstamp,units="secs")
      
      if (is.null(callres)) {
        parresi[,k,2] <- c(rep(NA,nmeas-1), # numAct
                           tmp_time) 
      } else {
        ypred <- round(callres$yhat)
        parresi[,k,2] <- c(mean((callres$yhat-ytest)^2)/rMSPE_const, # rMSPE
                           sum(family$dev.resids(ytest,callres$yhat,1)) / rDev_const, 
                           sum(family$dev.resids(y,callres$yhat_tr,1)) /rDev_const_tr,
                           sum(callres$beta!=0), # numAct
                           tmp_time)
        
      }
      # message(sprintf("Finished method %d / %d!",k,nmethods))
    }
    family <- quasipoisson(log)
    rDev_const <- sum(family$dev.resids(ytest,rep(family$linkinv(coef(glm(y~1,family = family,start=1))),ntest),1))
    rDev_const_tr <- sum(family$dev.resids(y,rep(family$linkinv(coef(glm(y~1,family = family,start=1))),n),1))
    for (k in which(!names(methods)%in%c("SIS","SVM","RF"))) {
      set.seed((1234+i)^2 + k)
      tstamp <- Sys.time()
      callres <- tryCatch( methods[[k]](x,y,xtest,family),
                           error=function(error_message) {
                             message(paste0("Error in ",names(methods)[k],": ",error_message))
                             return(NULL)
                           })
      tmp_time <- as.numeric(Sys.time() - tstamp,units="secs")
      
      if (is.null(callres)) {
        parresi[,k,3] <- c(rep(NA,nmeas-1), # numAct
                           tmp_time) 
      } else {
        ypred <- round(callres$yhat)
        parresi[,k,3] <- c(mean((callres$yhat-ytest)^2)/rMSPE_const, # rMSPE
                           sum(family$dev.resids(ytest,callres$yhat,1)) / rDev_const, 
                           sum(family$dev.resids(y,callres$yhat_tr,1)) /rDev_const_tr,
                           sum(callres$beta!=0), # numAct
                           tmp_time)
        
      }
      # message(sprintf("Finished method %d / %d!",k,nmethods))
    }
    # res[i,,,,j] <- parresi
    # res[i,,,,j]
    cat(sprintf('Finished rep %d / %d for setting %d / %d at %s.\n',i,nrep,j,nset,Sys.time()))
    warnings()
    parresi
}

for (j in 1:nset) {
  for (i in 1:nrep) {
    res[i,,,,j] <- parres[[j]][[i]]
  }
}

saveRDS(list(res=res,dataset_sizes=dataset_sizes), 
        file = sprintf("../saved_results/SPARglm_data_trib_nset%d_reps%d_nmeth%d.rds",nset,nrep,nmethods))

parallel::stopCluster(cl = my.cluster)

warnings()