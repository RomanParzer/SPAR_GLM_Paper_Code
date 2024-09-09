
pacman::p_load(SPAR,ggplot2,dplyr)

# 1st section binomial, lymphoma cancer and darwin alzeimer
# 2nd poisson tribology
# 3rd simulation group block and sparse ar1

# # # # # # binomial

# # # # darwin
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
# darwin_big <- list(x=cbind(darwin$x,darwin$x[sample(nrow(darwin$x)),],darwin$x[sample(nrow(darwin$x)),]),y=darwin$y)

feat_names <- sapply(colnames(darwin$x)[1:18],function(name)substr(name, 1, nchar(name)-1))
set.seed(1234)
spar_res_darwin <- spar.cv(darwin$x,darwin$y,family=binomial(logit),nummods=c(10,20,30,50))
plot(spar_res_darwin)
# ggsave(paste0("../plots/spar_lambda_darwin.pdf"), height = 6, width = 10)

# summary(as.numeric(abs(spar_res_darwin$betas)))
reorder_ind <- c(NULL)
for (i in 1:18) {
  reorder_ind <- c(reorder_ind,i + (0:24)*18)
}
plot(spar_res_darwin,"coefs",coef_order = reorder_ind) + 
  geom_vline(xintercept = 0.5 + 1:17*25,alpha=0.2,linetype=2) +
  annotate("text",x=0:17*25 + 12,y=40,label=feat_names,angle=90,size=3) +
  labs(fill="coef value")
# ggsave(paste0("../plots/spar_coef_darwin.pdf"), height = 4, width = 8)

# colnames(darwin$x)[reorder_ind]

# compare individual variable importance to RF and ElNet
coef_spar <- coef(spar_res_darwin)
darwin_importance <-matrix(c(0),18,3)
colnames(darwin_importance) <- c("SPAR","RF","ElNet")
rownames(darwin_importance) <- feat_names

# standardized scale
avg_coef_th <- colMeans(matrix((coef_spar$beta*spar_res_darwin$xscale)[reorder_ind],25,18))
names(avg_coef_th) <- feat_names
# round(avg_coef_th*100,2)

# closer look at 'number of pendowns'
avg_coef <- colMeans(matrix(rowMeans(as.matrix(spar_res_darwin$betas[reorder_ind,])),25,18))
names(avg_coef) <- feat_names
# round(avg_coef*100,2)
darwin_importance[,1] <- avg_coef*100
plot(spar_res_darwin,"coefs",coef_order = reorder_ind,prange = c(325,350)) +
  annotate("text",label=1:25,x=326:350,y=18)
344%%25
339%%25
# tasks 19 (copy postal code): , 11-14,16,17 slightly pos
# here, high number of pendowns is an indicator of healthy poeple, while on all tasks 1-10 + 15,18 it is rather associated with AD

# RF
RF_res <- randomForest::tuneRF(darwin$x,factor(darwin$y),ntreeTry = nrow(darwin$x), importance=TRUE,
                               doBest=TRUE, stepFactor = 1.5, improve=0.01,trace=FALSE)
# str(RF_res)

# randomForest::importance(RF_res)
avg_impM <- apply(randomForest::importance(RF_res)[reorder_ind,],2,function(tmp){
  avg_imp <- colMeans(matrix(tmp,25,18))
  round(avg_imp,2)
})
row.names(avg_impM) <- feat_names
# avg_impM
darwin_importance[,2] <- avg_impM[,3]
# RF_res$importance
# randomForest::importance(RF_res)[reorder_ind[c(330:350)],]

# can not tell the direction, importance similar

# glmnet

elnet_res <- glmnet::cv.glmnet(darwin$x,darwin$y,alpha=3/4,family="binomial")
elnet_res <- glmnet::glmnet(darwin$x,darwin$y,alpha=3/4,family="binomial",lambda=elnet_res$lambda.min)
# str(elnet_res)
elnet_coef <- colMeans(matrix((elnet_res$beta*spar_res_darwin$xscale)[reorder_ind],25,18))
names(elnet_coef) <- feat_names
# round(elnet_coef*100,2)
darwin_importance[,3] <- elnet_coef*100

round(darwin_importance,2)
# similar, ElNet sparser (eg no papertime), only different sign: mean_jerk_on_paper & max_x__extension (both rather low abs value)
# RF only gives positive value, also similar

# # # #  lymphoma

DLBCL <- read.csv("../data/lymphoma/DLBCL",skip=5475,header=FALSE)
# str(DLBCL)
# DLBCL$V5470
lymphoma <- list(x=as.matrix(DLBCL[,-5470]),y=DLBCL[,5470])
# mean(lymphoma$y)

set.seed(1234)
lymphoma_big <- list(x=cbind(lymphoma$x,
                             lymphoma$x[sample(length(lymphoma$y)),],
                             lymphoma$x[sample(length(lymphoma$y)),]),
                     y=lymphoma$y)



set.seed(1234)
spar_res_lymphoma_big <- spar.cv(lymphoma_big$x,lymphoma_big$y,family=binomial(logit),
                                 nummods=c(10,20,30,50))

# may need nscreen more than 2n for nicer plot of coefs

plot(spar_res_lymphoma_big)
plot(spar_res_lymphoma_big,"coefs") +
  geom_vline(xintercept=ncol(lymphoma$x)+0.5,alpha=0.2,linetype=2) +
  labs(fill="coef value") +
  scale_fill_gradient2(transform = scales::transform_modulus(p=1/10))
# ggsave(paste0("../plots/spar_coef_lymphoma_big_mod.pdf"), height = 4, width = 8)

plot(spar_res_lymphoma_big,"coefs") +
  geom_vline(xintercept=ncol(lymphoma$x)+0.5,alpha=0.2,linetype=2) +
  labs(fill="coef value") 
# ggsave(paste0("../plots/spar_coef_lymphoma_big.pdf"), height = 4, width = 8)

data.frame(predictor=1:ncol(lymphoma_big$x),sum_abs_coef_value=apply(abs(spar_res_lymphoma_big$betas),1,sum)) %>%
  ggplot(aes(x=predictor,y=sum_abs_coef_value)) +
  geom_point() +
  geom_vline(xintercept=ncol(lymphoma$x)+0.5,alpha=0.2,linetype=2) +
  labs(y="sum of absolute coefficients")
# ggsave(paste0("../plots/spar_coef_lymphoma_big_abs.pdf"), height = 4, width = 8)

# # # # # # poisson / gaussian log
# # # # tribology

diff_spectra <- read.csv("../data/tribology/diff_spectra.csv",header=TRUE,row.names = 1)
diff_spectra[,1] <- ifelse(diff_spectra[,1]=="large",1,0)
colnames(diff_spectra)[-1] <- sapply(colnames(diff_spectra)[-1],function(name) substring(name,2))
duration_spectra <- list(x=as.matrix(diff_spectra),
                         y=read.csv("../data/tribology/duration.csv")$x)
set.seed(1234)
spar_res_trib <- spar.cv(duration_spectra$x,duration_spectra$y,family = gaussian(log),
                         nummods = c(10,20,30,50),control=list(scr=list(nscreen = 2*length(duration_spectra$y))))
plot(spar_res_trib)
plot(spar_res_trib,"coef") +
  labs(fill="coef value")+
  scale_y_continuous(breaks=c(2,4,6,8,10,12))
# ggsave("../plots/spar_coef_tribology.pdf",height = 4,width = 8)

# Pfeiffer et al (2022) had problems to fit high durations, look at y vs yfit here (just training data)
plot(spar_res_trib,"res",xfit=duration_spectra$x,yfit=duration_spectra$y)

data.frame(y=duration_spectra$y,yfit = predict(spar_res_trib,duration_spectra$x)) %>%
  ggplot(aes(x=yfit,y=y)) +
  geom_point() + 
  geom_abline(slope=1,intercept = 0)

break_labels <- c(1000,2000,3000,4000)
break_inds <- sapply(break_labels,function(brk)which.min(abs(as.numeric(names(duration_spectra$x[1,-1]))-brk) ))

# # 3030 - 2770 cm-1 and 1480-1430 
# names(duration_spectra$x[1,-1])[c(503,636,1306,1331)]
# str(spar_res_trib$betas)
plot(spar_res_trib,"coef",prange = c(2,1815)) +
  labs(x="wavenumber") +
  scale_x_continuous(breaks=break_inds,labels=break_labels) +
  geom_vline(xintercept = c(503,636,1306,1331),linetype=2,alpha=0.2) +
  labs(fill="coef value") +
  scale_y_continuous(breaks=c(2,4,6,8,10,12))
# ggsave("../plots/spar_coef_tribology_wavenumber.pdf",height = 3,width = 8)

# first full absorption area
names(duration_spectra$x[1,-1])[503:636] # exact indizes in range
plot(spar_res_trib,"coef",prange=c(500,638)) +
  labs(x="wavenumber") +
  scale_x_continuous(breaks=seq(500,638,length.out=10),labels=names(duration_spectra$x[1,-1])[seq(500,638,length.out=10)])

# second full absorption area
names(duration_spectra$x[1,-1])[1306:1331] # exact indizes in range
plot(spar_res_trib,"coef",prange=c(1300,1335)) +
  labs(x="wavenumber") +
  scale_x_continuous(breaks=seq(1300,1335,length.out=10),labels=names(duration_spectra$x[1,-1])[seq(1300,1335,length.out=10)])


duration_spectra$y[1]
data.frame(`difference absorbance`=duration_spectra$x[1,-1],
           wavenumber = as.numeric(names(duration_spectra$x[1,-1]))) %>%
  ggplot(aes(x=wavenumber,y=difference.absorbance)) +
  geom_line() +
  scale_x_reverse() +
  geom_vline(xintercept = c(3030,2770,1480,1430),linetype=2,alpha=0.2)
# ggsave("../plots/x_observation_tribology.pdf",height = 3,width = 8)



# # # # # # simulations
source("../functions/glm_data_generation.R")
# sparse group
act_inds <- c(40:50,75,125,140:150,175,215:225,250)
beta <- numeric(300)
beta[c(40:50,140:150,215:225)] <- 1
beta[c(75,125,175,250)] <- -3
set.seed(1234)
data <- generate_data_glm(200,300,cov_setting = "group",family=poisson(log),ind = act_inds,beta=beta,
                          signal_strength = 1/2, avg_exp = 5)
# c(mean(data$y),var(data$y)) # model correct, not overdispersed

spar_res <- spar.cv(data$x,data$y,family=poisson(log),nummods=c(10,20,30,50))
plot(spar_res)
plot(spar_res,"coefs",prange = c(30,260)) +
  geom_vline(xintercept = act_inds[c(1,11,12,13,14,24,25,26,36,37)],alpha=0.5,linetype=2) +
  annotate("text",x=c(45,145,220),y=40,angle=c(90,90,90),label="active") +
  labs(fill="coef value")
# ggsave(paste0("../plots/spar_coef_sparse_group.pdf"), height = 4, width = 8)
