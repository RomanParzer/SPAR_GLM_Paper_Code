## read in file and plot  -----------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------

pacman::p_load(dplyr, ggplot2, tidyr, ggrepel,knitr,kableExtra)

resobj <- readRDS("../saved_results/SPARglm_M_Ensemble_screen_nset15_reps100_nmeth7.rds")
res <- resobj

methods <- dimnames(res)[[3]]
settings <- attributes(res)$settings
settings$family <- sapply(settings$family,function(obj) paste0(obj$family,"(",obj$link,")"))


mydf_all <- data.frame(pivot_longer(data.frame(res[,1,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="AUC"),
                       bAcc=pivot_longer(data.frame(res[,2,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="bAcc")$bAcc,
                       Acc=pivot_longer(data.frame(res[,3,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Acc")$Acc,
                       rMSPE=pivot_longer(data.frame(res[,4,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="rMSPE")$rMSPE,
                       rDev=pivot_longer(data.frame(res[,5,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="rDev")$rDev,
                       rDev_tr=pivot_longer(data.frame(res[,6,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="rDev_tr")$rDev_tr,
                       rMSLE=pivot_longer(data.frame(res[,7,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="rMSLE")$rMSLE,
                       pAUC=pivot_longer(data.frame(res[,8,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="pAUC")$pAUC,
                       Precision=pivot_longer(data.frame(res[,9,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Precision")$Precision,
                       Recall=pivot_longer(data.frame(res[,10,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Recall")$Recall,
                       Sign_ratio_Scr=pivot_longer(data.frame(res[,11,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Sign_ratio_Scr")$Sign_ratio_Scr,
                       Cor_Scr=pivot_longer(data.frame(res[,12,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Cor_Scr")$Cor_Scr,
                       NumAct=pivot_longer(data.frame(res[,13,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="NumAct")$NumAct,
                       Time=pivot_longer(data.frame(res[,14,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Time")$Time,
                       settings[1,],
                       setting=1)
for (k in 2:nrow(settings)) {
  mydf_all <- rbind(mydf_all,
                    data.frame(pivot_longer(data.frame(res[,1,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="AUC"),
                               bAcc=pivot_longer(data.frame(res[,2,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="bAcc")$bAcc,
                               Acc=pivot_longer(data.frame(res[,3,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Acc")$Acc,
                               rMSPE=pivot_longer(data.frame(res[,4,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="rMSPE")$rMSPE,
                               rDev=pivot_longer(data.frame(res[,5,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="rDev")$rDev,
                               rDev_tr=pivot_longer(data.frame(res[,6,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="rDev_tr")$rDev_tr,
                               rMSLE=pivot_longer(data.frame(res[,7,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="rMSLE")$rMSLE,
                               pAUC=pivot_longer(data.frame(res[,8,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="pAUC")$pAUC,
                               Precision=pivot_longer(data.frame(res[,9,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Precision")$Precision,
                               Recall=pivot_longer(data.frame(res[,10,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Recall")$Recall,
                               Sign_ratio_Scr=pivot_longer(data.frame(res[,11,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Sign_ratio_Scr")$Sign_ratio_Scr,
                               Cor_Scr=pivot_longer(data.frame(res[,12,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Cor_Scr")$Cor_Scr,
                               NumAct=pivot_longer(data.frame(res[,13,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="NumAct")$NumAct,
                               Time=pivot_longer(data.frame(res[,14,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Time")$Time,
                               settings[k,],
                               setting=k)
  )
}

# mydf_all$Method <- stringr::str_replace_all(mydf_all$Method,"\\.","")

mydf_all$Method <- factor(mydf_all$Method,levels=methods)

# plot binomial class measures
mydf_all %>% filter(Method %in% methods,family=="binomial(logit)") %>%
  pivot_longer(c(Acc,bAcc,AUC,rMSPE,rDev),names_to = "Measure",values_to = "Value") %>% 
  ggplot(aes(x=Method,y=Value,fill=Method)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  # ggh4x::facet_grid2(cov_setting~act_setting, scales = "free_y",independent = "y") +
  facet_grid(Measure~act_setting, scales = "free_y") +
  theme(legend.position = "none") #+

# select AUC
mydf_all <- mutate(mydf_all, 
                   pred_error = ifelse(family %in% c("binomial(logit)","binomial(cloglog)"),1-AUC,rMSPE))



short_fam_names1 <- labeller(
  family=c(`binomial(cloglog)` = "bin(cll)", `binomial(logit)` = "bin(logit)",`gaussian(identity)` = "gau(id)", 
           `gaussian(log)` = "gau(log)",`poisson(log)` = "poi(log)"),
  act_setting=c(`sparse`="sparse",`medium`="medium",`dense`="dense"))


# plot families pred
mydf_all %>% filter(Method %in% methods) %>%
  ggplot(aes(x=Method,y=rMSPE,fill=Method)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  facet_grid(family~act_setting, scales = "free_y",labeller =short_fam_names1) +
  # coord_cartesian(ylim=c(0,1.0)) +
  theme(legend.position = "none") +
  labs(y="prediction error")
# ggsave(paste0("../plots/glm_pred_error_m_ensemble.pdf"), height = 6, width = 8)

myqnorm <- qnorm(0.975)
# plot mean ranks per facet
rank_df <- mydf_all %>% group_by(rep,setting) %>%
  mutate(rank=rank(pred_error))
rank_df %>% group_by(Method,act_setting,family) %>%
  summarize(meanRank = mean(rank),sdRank=sd(rank)) %>%
  ggplot(aes(x=Method,y=meanRank,col=Method)) +
  geom_point() +
  geom_errorbar(aes(ymin = meanRank - myqnorm*sdRank/sqrt(100),ymax = meanRank + myqnorm*sdRank/sqrt(100))) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  facet_grid(family~act_setting, scales = "free_y",labeller =short_fam_names1) +
  # coord_cartesian(ylim=c(0,1.0)) +
  theme(legend.position = "none") +
  labs(y="Mean rank for prediction")
# ggsave(paste0("../plots/glm_pred_rank_m_ensemble.pdf"), height = 6, width = 8)
