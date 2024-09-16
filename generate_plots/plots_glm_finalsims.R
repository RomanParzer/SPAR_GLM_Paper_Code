## read in file and plot  -----------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------

pacman::p_load(dplyr, ggplot2, tidyr, ggrepel,knitr,kableExtra)

resobj <- readRDS("../saved_results/SPARglm_sims_nset40_reps100_nmeth19.rds")
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


# plot binomail class measures
mydf_all %>% filter(Method %in% methods,family=="binomial(logit)",
                    p==2000,cov_setting=="group") %>%
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
# str(mydf_all)
mydf_all$act_setting <- factor(mydf_all$act_setting,levels = c("sparse","medium","dense"))

mydf_all$Method <- stringr::str_replace_all(mydf_all$Method,"TARP.","TARP+")
mydf_all$Method <- stringr::str_replace_all(mydf_all$Method,"\\."," ")
mydf_all$Method <- factor(mydf_all$Method,levels = methods)

# mydf_all <- mutate(mydf_all,"isSparse"=ifelse(Method%in%c("AdLASSO","ElNet","SIS"),TRUE,FALSE))
# ,linetype=isSparse


# plot families pred
mydf_all %>% filter(Method %in% methods,p==2000,cov_setting=="group") %>%
  ggplot(aes(x=Method,y=pred_error,fill=Method)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  # ggh4x::facet_grid2(cov_setting~act_setting, scales = "free_y",independent = "y") +
  facet_grid(family~act_setting, scales = "free_y") +
  # coord_cartesian(ylim=c(0,1.0)) +
  theme(legend.position = "none") +
  labs(y="prediction error")

# select methods to show
methods
show_methods <- methods[c(1,3,4,6,7,9,12,15,18)]
scr_coef_meth <- methods[c(1,2,3,4,5,12,15,16)]
rank_methods <- methods[c(1,2,3,4,5,6,7,9,12,15,18)]
show_families <- c("binomial(logit)","gaussian(identity)","poisson(log)")


# plot families pred
mydf_all %>% filter(Method %in% show_methods,p==2000,cov_setting=="group") %>%
                    #family%in%show_families) %>%
  ggplot(aes(x=Method,y=pred_error,fill=Method)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  # ggh4x::facet_grid2(cov_setting~act_setting, scales = "free_y",independent = "y") +
  facet_grid(family~act_setting, scales = "free_y") +
  # coord_cartesian(ylim=c(0,1.0)) +
  theme(legend.position = "none") +
  labs(y="prediction error")
# ggsave(paste0("../plots/glm_pred_error_families.pdf"), height = 6, width = 8)


short_fam_names1 <- labeller(
  family=c(`binomial(cloglog)` = "bin(cll)", `binomial(logit)` = "bin(logit)",`gaussian(identity)` = "gau(id)", 
           `gaussian(log)` = "gau(log)",`poisson(log)` = "poi(log)"),
  cov_setting=c(`ar1`="ar1",`comsym`="comsym",`group`="group",`ind`="ind"))

mydf_all %>% filter(Method %in% show_methods,p==2000,act_setting=="medium",
                    #family%in%show_families,
                    pred_error<1.5) %>%
  ggplot(aes(x=Method,y=pred_error,fill=Method)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  # ggh4x::facet_grid2(cov_setting~act_setting, scales = "free_y",independent = "y") +
  facet_grid(family~cov_setting, scales = "free_y",labeller = short_fam_names1) +
  # coord_cartesian(ylim=c(0,1.0)) +
  theme(legend.position = "none") +
  labs(y="prediction error")
# ggsave(paste0("../plots/glm_pred_error_cov_settings_med.pdf"), height = 5, width = 8)

short_fam_names2 <- labeller(
  family=c(`binomial(cloglog)` = "bin(cll)", `binomial(logit)` = "bin(logit)",`gaussian(identity)` = "gau(id)", 
    `gaussian(log)` = "gau(log)",`poisson(log)` = "poi(log)"),
  p=c(`500`="500",`2000`="2000",`10000`="10000"))

mydf_all %>% filter(Method %in% show_methods,cov_setting=="group",act_setting=="medium",
                    #family%in%show_families,
                    pred_error < 1.5) %>%
  ggplot(aes(x=Method,y=pred_error,fill=Method)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  # ggh4x::facet_grid2(cov_setting~act_setting, scales = "free_y",independent = "y") +
  facet_grid(family~p, scales = "free_y",labeller = short_fam_names2) +
  # coord_cartesian(ylim=c(0,1.0)) +
  theme(legend.position = "none") +
  labs(y="prediction error")
# ggsave(paste0("../plots/glm_pred_error_p.pdf"), height = 5, width = 8)

mydf_all %>% filter(cov_setting=="group",act_setting=="medium",
                    family=="poisson(log)",p==500,
                    pred_error<2) %>%
  ggplot(aes(x=Method,y=pred_error,fill=Method)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  # ggh4x::facet_grid2(cov_setting~act_setting, scales = "free_y",independent = "y") +
  # facet_grid(family~p, scales = "free_y",labeller = short_fam_names2) +
  coord_cartesian(ylim=c(0,1.75)) +
  theme(legend.position = "none") +
  labs(y="prediction error")


# plot families rDev
mydf_all %>% filter(Method %in% show_methods,p==2000,cov_setting=="group",
                    family%in%show_families) %>%
  ggplot(aes(x=Method,y=rDev,fill=Method)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  # ggh4x::facet_grid2(cov_setting~act_setting, scales = "free_y",independent = "y") +
  facet_grid(family~act_setting, scales = "free_y") +
  # coord_cartesian(ylim=c(0,1.0)) +
  theme(legend.position = "none")

# plot families rDev_tr
mydf_all %>% filter(Method %in% show_methods,p==2000,cov_setting=="group",
                    family%in%show_families) %>%
  ggplot(aes(x=Method,y=rDev_tr,fill=Method)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  # ggh4x::facet_grid2(cov_setting~act_setting, scales = "free_y",independent = "y") +
  facet_grid(family~act_setting, scales = "free_y") +
  coord_cartesian(ylim=c(0,1.1)) +
  theme(legend.position = "none") #+


# plot families rMSLE
mydf_all %>% filter(Method %in% show_methods,Method!="RF",Method!="SVM",
                    #family%in%show_families,
                    p==2000,cov_setting=="group") %>%
  ggplot(aes(x=Method,y=rMSLE,fill=Method)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  # ggh4x::facet_grid2(family~act_setting, scales = "free_y",independent = "y") +
  facet_grid(family~act_setting, scales = "free_y") +
  theme(legend.position = "none") #+
# ggsave(paste0("../plots/glm_rMSLE_families.pdf"), height = 6, width = 8)

# plot families pAUC
mydf_all %>% filter(Method %in% show_methods,p==2000,cov_setting=="group",
                    family%in%show_families) %>%
  ggplot(aes(x=Method,y=pAUC,fill=Method)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggh4x::facet_grid2(family~act_setting, scales = "free_y",independent = "y") +
  # facet_grid(family~act_setting+p, scales = "free_y") +
  theme(legend.position = "none") #+
# ggsave(paste0("../plots/glm_pAUC_families.pdf"), height = 5, width = 8)


# plot families Prec
mydf_all %>% filter(Method %in% show_methods,p==2000,cov_setting=="group",
                    family%in%show_families) %>%
  ggplot(aes(x=Method,y=Precision,fill=Method)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  # ggh4x::facet_grid2(cov_setting~act_setting, scales = "free_y",independent = "y") +
  facet_grid(family~act_setting, scales = "free_y") +
  theme(legend.position = "none") #+

# plot families Rec
mydf_all %>% filter(Method %in% show_methods,p==2000,cov_setting=="group",
                    family%in%show_families) %>%
  ggplot(aes(x=Method,y=Recall,fill=Method)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  # ggh4x::facet_grid2(cov_setting~act_setting, scales = "free_y",independent = "y") +
  facet_grid(family~act_setting, scales = "free_y") +
  theme(legend.position = "none") #

# plot families Sign_ratio_Scr
Scr_df <- mydf_all %>% filter(Method %in% scr_coef_meth,p==2000,cov_setting=="group",family%in%show_families)
Scr_df$Method <- factor(Scr_df$Method,levels = scr_coef_meth[c(2:8,1)])

Scr_df %>%  
  ggplot(aes(x=Method,y=Sign_ratio_Scr,fill=Method)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  # ggh4x::facet_grid2(cov_setting~act_setting, scales = "free_y",independent = "y") +
  facet_grid(family~act_setting, scales = "free_y") +
  theme(legend.position = "none") #+´

# plot families Cor_Scr
Scr_df %>% 
  ggplot(aes(x=Method,y=Cor_Scr,fill=Method)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  facet_grid(family~act_setting, scales = "free_y") +
  theme(legend.position = "none") +
  labs(y= "Cor") +
  scale_x_discrete(labels=rep(c("LASSO","AdLASSO","ElNet","SIS","marGLM","Ridge0","Rid0_sm_n","Ridge"),1))
# already inspected earlier in extrea simulations

# plot families NumAct
mydf_all %>% filter(Method %in% show_methods,p==2000,cov_setting=="group",family%in%show_families) %>%
  ggplot(aes(x=Method,y=NumAct,fill=Method)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  # ggh4x::facet_grid2(cov_setting~act_setting, scales = "free_y",independent = "y") +
  facet_grid(family~act_setting, scales = "free_y") +
  theme(legend.position = "none") #+
# ggsave(paste0("../plots/glm_NumAct_families.pdf"), height = 6, width = 8)

# plot families Time
mydf_all %>% filter(Method %in% show_methods,p==2000,cov_setting=="group",family%in%show_families) %>%
  ggplot(aes(x=Method,y=Time,color=Method)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  # ggh4x::facet_grid2(cov_setting~act_setting, scales = "free_y",independent = "y") +
  facet_grid(.~family, scales = "free_y") +
  theme(legend.position = "none") +
  scale_y_log10()


myp <- c(500,2000,10^4)
mydf_time_sum_fam <- mydf_all %>% filter(Method %in% rank_methods,act_setting=="medium",cov_setting=="group") %>% 
  group_by(Method,p,family) %>%
  summarise(Time=mean(Time)) %>% mutate(is_ref=FALSE)
mydf_time_sum_fam <- rbind(mydf_time_sum_fam,
                       data.frame(p=rep(myp,3),
                                  Time=c(log(myp)*myp/30,myp/30,log(myp)/10),
                                  Method=c(rep("O(plog(p))",length(myp)),rep("O(p)",length(myp)),rep("O(log(p))",length(myp))),
                                  is_ref=TRUE,family="poisson(log)"),
                       data.frame(p=rep(myp,3),
                                  Time=c(log(myp)*myp/30,myp/30,log(myp)/10),
                                  Method=c(rep("O(plog(p))",length(myp)),rep("O(p)",length(myp)),rep("O(log(p))",length(myp))),
                                  is_ref=TRUE,family="binomial(logit)"),
                       data.frame(p=rep(myp,3),
                                  Time=c(log(myp)*myp/30,myp/30,log(myp)/10),
                                  Method=c(rep("O(plog(p))",length(myp)),rep("O(p)",length(myp)),rep("O(log(p))",length(myp))),
                                  is_ref=TRUE,family="binomial(cloglog)"),
                       data.frame(p=rep(myp,3),
                                  Time=c(log(myp)*myp/30,myp/30,log(myp)/10),
                                  Method=c(rep("O(plog(p))",length(myp)),rep("O(p)",length(myp)),rep("O(log(p))",length(myp))),
                                  is_ref=TRUE,family="gaussian(log)"),
                       data.frame(p=rep(myp,3),
                                  Time=c(log(myp)*myp/30,myp/30,log(myp)/10),
                                  Method=c(rep("O(plog(p))",length(myp)),rep("O(p)",length(myp)),rep("O(log(p))",length(myp))),
                                  is_ref=TRUE,family="gaussian(identity)")
)
mydf_time_sum_fam$Method <- factor(mydf_time_sum_fam$Method,levels=c(rank_methods,"O(plog(p))","O(p)","O(log(p))"))
mydf_time_sum_fam %>% 
  ggplot(aes(x=p,y=Time,col=Method,linetype=is_ref)) +
  geom_line() +
  geom_point() +
  # facet_grid(family~., scales = "free_y") +
  scale_color_manual(values=c(scales::hue_pal()(12)[1:11],"darkgrey","darkgrey","darkgrey"))+ # careful about length!
  scale_x_log10() +
  scale_y_log10() +
  geom_text_repel(data = filter(mydf_time_sum_fam,p==10000),
                  aes(x=p,y=Time,label=Method),show.legend = FALSE) +
  theme(legend.position="none") +
  facet_grid(.~family) +
  labs(y="Time in s")
# ggsave(paste0("../plots/glm_Time_p_families.pdf"), height = 5, width = 10)

mydf_time_sum <- mydf_all %>% filter(Method %in% rank_methods,act_setting=="medium",
                                     # family%in%show_families,
                                     cov_setting=="group") %>% 
  group_by(Method,p) %>%
  summarise(Time=mean(Time)) %>% mutate(is_ref=FALSE)
mydf_time_sum <- rbind(mydf_time_sum,
                       data.frame(p=rep(myp,3),
                                  Time=c(log(myp)*myp/30,myp/20,log(myp)/3),
                                  Method=c(rep("O(plog(p))",length(myp)),rep("O(p)",length(myp)),rep("O(log(p))",length(myp))),
                                  is_ref=TRUE)
)
mydf_time_sum$Method <- factor(mydf_time_sum$Method,levels=c(rank_methods,"O(plog(p))","O(p)","O(log(p))"))
mydf_time_sum %>% 
  ggplot(aes(x=p,y=Time,col=Method,linetype=is_ref)) +
  geom_line() +
  geom_point() +
  # facet_grid(family~., scales = "free_y") +
  scale_color_manual(values=c(scales::hue_pal()(12)[1:11],"darkgrey","darkgrey","darkgrey"))+ # careful about length!
  scale_x_log10() +
  scale_y_log10() +
  geom_text_repel(data = filter(mydf_time_sum,p==10000),
                  aes(x=p,y=Time,label=Method),show.legend = FALSE) +
  theme(legend.position="none") +
  labs(y="Time in s")
# ggsave(paste0("../plots/glm_Time_p.pdf"), height = 4.6, width = 8)


# # rank tables

# # pred / LE / pAUC
# each act_setting

rank_methods <- methods
n_showm <- length(rank_methods)
mydf_all_rank <- mydf_all %>% filter(Method %in% rank_methods) %>% group_by(rep,setting) %>% mutate(rank_pred=rank(pred_error),rank_LE = rank(rMSLE), rank_pAUC = n_showm + 1 - rank(pAUC))
rank_tab <- mydf_all_rank %>% group_by(Method,family,act_setting) %>% summarise(mean_rank_pred = mean(rank_pred), se_rank_pred = sd(rank_pred)/sqrt(100),
                                                                                         mean_rank_LE = mean(rank_LE), se_rank_LE = sd(rank_LE)/sqrt(100),
                                                                                         mean_rank_pAUC = mean(rank_pAUC), se_rank_pAUC = sd(rank_pAUC)/sqrt(100)) 

filter(rank_tab,family=="binomial(logit)",act_setting=="medium")
filter(rank_tab,family=="poisson(log)",act_setting=="dense")
# together



OVtab <- mydf_all_rank %>% group_by(Method) %>% summarise(mean_rank_pred = mean(rank_pred), se_rank_pred = sd(rank_pred)/sqrt(100*length(unique(mydf_all_rank$setting))),
                                                 mean_rank_LE = mean(rank_LE), se_rank_LE = sd(rank_LE)/sqrt(100*length(unique(mydf_all_rank$setting))),
                                                 mean_rank_pAUC = mean(rank_pAUC), se_rank_pAUC = sd(rank_pAUC)/sqrt(100*length(unique(mydf_all_rank$setting)))) 

OVtab[,-1] <- round(OVtab[,-1],3)
OVtab[,c(3,5,7)] <- apply(OVtab[,c(3,5,7)],2,function(col)paste0("(",col,")"))
OVtab
kable(OVtab,format = "latex",booktabs=TRUE)


# special plots

mydf_all %>% filter(Method %in% show_methods,p==2000,cov_setting=="group",
                    family %in%show_families) %>%
  pivot_longer(c(pred_error,rMSLE),names_to = "Measure",values_to = "Value") %>% 
  ggplot(aes(x=Method,y=Value,fill=Method)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  # ggh4x::facet_grid2(family+Measure~act_setting, scales = "free_y",independent = "y") +
  facet_grid(family+Measure~act_setting, scales = "free_y") +
  theme(legend.position = "none") +
  labs(y=" ")
# ggsave(paste0("../plots/glm_ov_pred_families.pdf"), height = 8, width = 8)

