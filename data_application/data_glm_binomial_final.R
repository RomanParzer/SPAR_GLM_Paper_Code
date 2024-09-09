#####################################################
### Testing sparse projected averaged regression SPAR
#####################################################

pacman::p_load(foreach,parallel,tidyr,dplyr)
source("../functions/glm_methods.R")
source("../functions/multi_assign.R")

################ adapt #########################################################################

nrep <- 100

methods <- list("Ridge"=function(x,y,xtest,family){myElNet(x,y,xtest,alpha=0,family=family)},
                "Ridge auc"=function(x,y,xtest,family){myElNet(x,y,xtest,alpha=0,family=family,type.measure = "auc")},
                "LASSO"=function(x,y,xtest,family){myElNet(x,y,xtest,alpha=1,family=family)},
                "AdLASSO"=myAdLASSO,
                "AdLASSO auc"=function(x,y,xtest,family){myAdLASSO(x,y,xtest,family,type.measure = "auc")},
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
                "SPAR CV auc"=function(x,y,xtest,family){mySPAR(x,y,xtest,family=family,nummods=c(10,20,30,50),type.measure="1-auc",doCV=TRUE)}
                )

measures <- c("AUC","Acc","Sens","Spec","bAcc","rMSPE","rDev","rDev_tr","NumAct","Time")


DLBCL <- read.csv("../data/lymphoma/DLBCL",skip=5475,header=FALSE)
# str(DLBCL)
# DLBCL$V5470
lymphoma <- list(x=as.matrix(DLBCL[,-5470]),y=DLBCL[,5470])
# mean(lymphoma$y)

set.seed(1234)
lymphoma_big <- list(x=cbind(lymphoma$x,lymphoma$x[sample(length(lymphoma$y)),],lymphoma$x[sample(length(lymphoma$y)),]),y=lymphoma$y)

darwin_tmp <- read.csv("../data/darwin/data.csv",stringsAsFactors = TRUE)
darwin_orig <- list(x=as.matrix(darwin_tmp[,-c(1,452)]),y=as.numeric(darwin_tmp$class)-1)
# implausible values (time, pressure 0), impute outlying values by DDC Roussew
# for (i in 1:18) {
#   print(colnames(darwin_orig$x)[i])
#   print(summary(c(darwin$x[,c(i + (0:24)*18)])))
# }
tmp <- cellWise::DDC(darwin_orig$x,list(returnBigXimp=TRUE,tolProb=0.999,silent=TRUE))
# summary(tmp$remX[tmp$indcells])
darwin <- list(x=tmp$Ximp,y=darwin_orig$y)
set.seed(1234)
darwin_big <- list(x=cbind(darwin$x,darwin$x[sample(nrow(darwin$x)),],darwin$x[sample(nrow(darwin$x)),]),y=darwin$y)

datasets <- list(lymphoma=lymphoma,
                 lymphoma_big=lymphoma_big,
                 darwin=darwin,
                 darwin_big=darwin_big)


dataset_sizes <- tibble(dataset=names(datasets),
                        n = sapply(datasets,function(arg) nrow(arg$x) ),
                        p = sapply(datasets,function(arg) ncol(arg$x) ))


################## simulations #############################################################################

nmethods <- length(methods)
nset <- length(datasets)
nmeas <- length(measures)
res <- array(c(0),dim=c(nrep,nmeas,nmethods,2,nset),dimnames = list(reps=NULL,measures=measures,
                                                                  method=names(methods),link=c("logit","cloglog"),datasets=names(datasets)))

unlink("../saved_results/log.txt")
my.cluster <- parallel::makeCluster(7, type = "PSOCK", outfile = "../saved_results/log.txt")
doParallel::registerDoParallel(cl = my.cluster)
# foreach::getDoParRegistered()
clusterExport(my.cluster,c('datasets','nmeas','nmethods','methods'), envir = environment())
clusterEvalQ(my.cluster, {  
  source("../functions/multi_assign.R")
  source("../functions/glm_methods.R")
})
# i <- j <- 1
parres <- foreach(j = 1:(nset)) %:%
foreach(i=1:nrep) %dopar% {
  parresi <- array(c(0),dim=c(nmeas,nmethods,2))
  
    mydata <- datasets[[j]]
    dataname <- names(datasets)[j]
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
    
    family <- binomial(logit)
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
        parresi[,k,1] <- c(my_AUC(ytest,callres$yhat),
                           mean(ytest==ypred),
                           mean(ypred[ytest==1]==1),
                           mean(ypred[ytest==0]==0),
                           (mean(ypred[ytest==1]==1) + mean(ypred[ytest==0]==0))/2,
                         mean((callres$yhat-ytest)^2)/rMSPE_const, # rMSPE
                         sum(family$dev.resids(ytest,callres$yhat,1)) / rDev_const, 
                         sum(family$dev.resids(y,callres$yhat_tr,1)) /rDev_const_tr,
                         sum(callres$beta!=0), # numAct
                         tmp_time)
        
      }
    }
    
    if (dataname%in% c("lymphoma","lymphoma_big")) {
      family <- binomial(cloglog)
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
          parresi[,k,2] <- c(rep(NA,nmeas-1), # numAct
                             tmp_time) 
        } else {
          ypred <- round(callres$yhat)
          parresi[,k,2] <- c(my_AUC(ytest,callres$yhat),
                             mean(ytest==ypred),
                             mean(ypred[ytest==1]==1),
                             mean(ypred[ytest==0]==0),
                             (mean(ypred[ytest==1]==1) + mean(ypred[ytest==0]==0))/2,
                             mean((callres$yhat-ytest)^2)/rMSPE_const, # rMSPE
                             sum(family$dev.resids(ytest,callres$yhat,1)) / rDev_const, 
                             sum(family$dev.resids(y,callres$yhat_tr,1)) /rDev_const_tr,
                             sum(callres$beta!=0), # numAct
                             tmp_time)
          
        }
      }
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
        file = sprintf("../saved_results/SPARglm_data_binom_nset%d_reps%d_nmeth%d.rds",nset,nrep,nmethods))

parallel::stopCluster(cl = my.cluster)

warnings()