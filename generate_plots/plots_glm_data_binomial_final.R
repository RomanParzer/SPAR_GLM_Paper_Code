## read in file and plot  -----------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------

pacman::p_load(dplyr, ggplot2, tidyr, ggrepel,knitr,kableExtra)

resobj <- readRDS("../saved_results/SPARglm_data_binom_nset4_reps100_nmeth20.rds")
res <- resobj$res
methods <- dimnames(res)[[3]]
resobj$dataset_sizes


mydf_all <- data.frame(pivot_longer(data.frame(res[,1,,1,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="AUC"),
                       Acc=pivot_longer(data.frame(res[,2,,1,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Acc")$Acc,
                       Sens=pivot_longer(data.frame(res[,3,,1,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Sens")$Sens,
                       Spec=pivot_longer(data.frame(res[,4,,1,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Spec")$Spec,
                       bAcc=pivot_longer(data.frame(res[,5,,1,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="bAcc")$bAcc,
                       rMSPE=pivot_longer(data.frame(res[,6,,1,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="rMSPE")$rMSPE,
                       rDev=pivot_longer(data.frame(res[,7,,1,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="rDev")$rDev,
                       rDev_tr=pivot_longer(data.frame(res[,8,,1,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="rDev_tr")$rDev_tr,
                       NumAct=pivot_longer(data.frame(res[,9,,1,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="NumAct")$NumAct,
                       Time=pivot_longer(data.frame(res[,10,,1,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Time")$Time,
                       resobj$dataset_sizes[1,],
                       link="logit")

mydf_all <- rbind(mydf_all,
      data.frame(pivot_longer(data.frame(res[,1,,2,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="AUC"),
                 Acc=pivot_longer(data.frame(res[,2,,2,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Acc")$Acc,
                 Sens=pivot_longer(data.frame(res[,3,,2,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Sens")$Sens,
                 Spec=pivot_longer(data.frame(res[,4,,2,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Spec")$Spec,
                 bAcc=pivot_longer(data.frame(res[,5,,2,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="bAcc")$bAcc,
                 rMSPE=pivot_longer(data.frame(res[,6,,2,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="rMSPE")$rMSPE,
                 rDev=pivot_longer(data.frame(res[,7,,2,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="rDev")$rDev,
                 rDev_tr=pivot_longer(data.frame(res[,8,,2,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="rDev_tr")$rDev_tr,
                 NumAct=pivot_longer(data.frame(res[,9,,2,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="NumAct")$NumAct,
                 Time=pivot_longer(data.frame(res[,10,,2,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Time")$Time,
                 resobj$dataset_sizes[1,],
                 link="cloglog"))
      
      
for (k in 2:nrow(resobj$dataset_sizes)) {
  mydf_all <- rbind(mydf_all,
                    data.frame(pivot_longer(data.frame(res[,1,,1,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="AUC"),
                               Acc=pivot_longer(data.frame(res[,2,,1,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Acc")$Acc,
                               Sens=pivot_longer(data.frame(res[,3,,1,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Sens")$Sens,
                               Spec=pivot_longer(data.frame(res[,4,,1,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Spec")$Spec,
                               bAcc=pivot_longer(data.frame(res[,5,,1,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="bAcc")$bAcc,
                               rMSPE=pivot_longer(data.frame(res[,6,,1,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="rMSPE")$rMSPE,
                               rDev=pivot_longer(data.frame(res[,7,,1,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="rDev")$rDev,
                               rDev_tr=pivot_longer(data.frame(res[,8,,1,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="rDev_tr")$rDev_tr,
                               NumAct=pivot_longer(data.frame(res[,9,,1,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="NumAct")$NumAct,
                               Time=pivot_longer(data.frame(res[,10,,1,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Time")$Time,
                               resobj$dataset_sizes[k,],
                               link="logit"))
  if (k<=2) {
    mydf_all <- rbind(mydf_all,
                      data.frame(pivot_longer(data.frame(res[,1,,2,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="AUC"),
                                 Acc=pivot_longer(data.frame(res[,2,,2,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Acc")$Acc,
                                 Sens=pivot_longer(data.frame(res[,3,,2,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Sens")$Sens,
                                 Spec=pivot_longer(data.frame(res[,4,,2,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Spec")$Spec,
                                 bAcc=pivot_longer(data.frame(res[,5,,2,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="bAcc")$bAcc,
                                 rMSPE=pivot_longer(data.frame(res[,6,,2,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="rMSPE")$rMSPE,
                                 rDev=pivot_longer(data.frame(res[,7,,2,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="rDev")$rDev,
                                 rDev_tr=pivot_longer(data.frame(res[,8,,2,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="rDev_tr")$rDev_tr,
                                 NumAct=pivot_longer(data.frame(res[,9,,2,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="NumAct")$NumAct,
                                 Time=pivot_longer(data.frame(res[,10,,2,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Time")$Time,
                                 resobj$dataset_sizes[k,],
                                 link="cloglog"))
  }
  
}

mydf_all$Method <- stringr::str_replace_all(mydf_all$Method,"\\."," ")
mydf_all$Method <- factor(mydf_all$Method,levels = methods)
mydf_all$dataset <- factor(mydf_all$dataset,levels = resobj$dataset_sizes$dataset)
# mydf_all <- mutate(mydf_all,"isSparse"=ifelse(Method%in%c("AdLASSO","ElNet","SIS"),TRUE,FALSE))
# ,linetype=isSparse
methods
# show_methods <- methods
# split, AdLASSO always worse here
show_methods <- methods[c(1,3,6,7,8,9,11,14,15,18)]
rank_methods <- methods[c(1,3,4,6,7,8,9,11,14,15,18)]

mydf_all %>% filter(Method %in% show_methods,dataset%in%c("lymphoma","lymphoma_big")) %>%
  ggplot(aes(x=link,y=rMSPE,fill=Method)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  # ggh4x::facet_grid2(cov_setting~act_setting, scales = "free_y",independent = "y") +
  facet_grid(dataset~Method, scales = "free_y") +
  theme(legend.position = "none")
# ggsave(paste0("../plots/glm_databinom_lymphoma_link_rMSPE.pdf"), height = 6, width = 16)
# logit always better, only for TARP on lymphoma_big, the cloglog link is better than logit
# so always use logit


mydf_all %>% filter(Method %in% show_methods,
                    dataset!="darwin_big",
                    link=="logit") %>%
  ggplot(aes(x=Method,y=AUC,fill=Method)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  # ggh4x::facet_grid2(cov_setting~act_setting, scales = "free_y",independent = "y") +
  facet_grid(.~dataset, scales = "free_y") +
  theme(legend.position = "none") +
  coord_cartesian(ylim=c(0.7,1.0))
# ggsave(paste0("../plots/glm_databinom_AUC.pdf"), height = 6, width = 8)

mydf_all %>% filter(Method %in% show_methods,
                    dataset!="darwin_big",
                    link=="logit") %>%
  ggplot(aes(x=Method,y=rMSPE,fill=Method)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  # ggh4x::facet_grid2(cov_setting~act_setting, scales = "free_y",independent = "y") +
  facet_grid(.~dataset, scales = "free_y") +
  theme(legend.position = "none") +
  coord_cartesian(ylim=c(0,1.0))


mydf_all %>% filter(Method %in% show_methods,
                    dataset!="darwin_big",
                    link=="logit") %>%
  ggplot(aes(x=Method,y=rDev,fill=Method)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  # ggh4x::facet_grid2(cov_setting~act_setting, scales = "free_y",independent = "y") +
  facet_grid(.~dataset, scales = "free_y") +
  theme(legend.position = "none") +
  coord_cartesian(ylim=c(0,1.0))

mydf_all %>% filter(Method %in% show_methods,
                    dataset!="darwin_big",
                    link=="logit") %>%
  ggplot(aes(x=Method,y=rDev_tr,fill=Method)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  # ggh4x::facet_grid2(cov_setting~act_setting, scales = "free_y",independent = "y") +
  facet_grid(.~dataset, scales = "free_y") +
  theme(legend.position = "none") +
  coord_cartesian(ylim=c(0,1.0))


mydf_all %>% filter(Method %in% show_methods,
                    dataset!="darwin_big",
                    link=="logit") %>%
  ggplot(aes(x=Method,y=bAcc,fill=Method)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  # ggh4x::facet_grid2(cov_setting~act_setting, scales = "free_y",independent = "y") +
  facet_grid(.~dataset, scales = "free_y") +
  theme(legend.position = "none") +
  coord_cartesian(ylim=c(0.4,1.0))

mydf_all %>% filter(Method %in% show_methods,
                    dataset!="darwin_big",
                    link=="logit") %>%
  ggplot(aes(x=Method,y=NumAct,fill=Method)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  # ggh4x::facet_grid2(cov_setting~act_setting, scales = "free_y",independent = "y") +
  facet_grid(dataset~., scales = "free_y") +
  theme(legend.position = "none")

mydf_all %>% filter(Method %in% show_methods,
                    dataset!="darwin_big",
                    link=="logit") %>%
  ggplot(aes(x=Method,y=Time,colour=Method)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  # ggh4x::facet_grid2(cov_setting~act_setting, scales = "free_y",independent = "y") +
  facet_grid(dataset~., scales = "free_y") +
  theme(legend.position = "none") +
  scale_y_log10()


# # rank tables, combine with tribology

# each dataset
n_showm <- length(rank_methods)
myrankdf <- mydf_all %>% filter(Method %in% rank_methods,
                                dataset!="darwin_big",
                                link=="logit") %>% 
  group_by(rep,dataset) %>% mutate(rank_AUC=n_showm + 1 - rank(AUC) ,rank_rMSPE=rank(rMSPE), rank_bAcc = n_showm + 1 - rank(bAcc))

myranksumdf <- myrankdf %>% group_by(Method,dataset) %>% summarise(mean_rank_AUC = mean(rank_AUC), se_rank_AUC = sd(rank_AUC)/sqrt(100),
                                                                   mean_rank_rMSPE = mean(rank_rMSPE), se_rank_rMSPE = sd(rank_rMSPE)/sqrt(100),
                                                                   mean_rank_bAcc = mean(rank_bAcc), se_rank_bAcc = sd(rank_bAcc)/sqrt(100)
                                                                   ) 
print(filter(myranksumdf,dataset=="darwin"),n=n_showm)
# together over all 3 data sets
print(myrankdf %>% group_by(Method) %>% summarise(mean_rank_AUC = mean(rank_AUC), se_rank_AUC = sd(rank_AUC)/sqrt(100*3),
                                                  mean_rank_rMSPE = mean(rank_rMSPE), se_rank_rMSPE = sd(rank_rMSPE)/sqrt(100*3),
                                                  mean_rank_bAcc = mean(rank_bAcc), se_rank_bAcc = sd(rank_bAcc)/sqrt(100*3)), 
      n=n_showm)

RankTab <- myranksumdf %>% pivot_wider(names_from = dataset,values_from = c(3:8),names_vary = "slowest")
RankTab[,-1] <- round(RankTab[,-1],3)
RankTab[,1+1:9*2] <- apply(RankTab[,1+1:9*2],2,function(col)paste0("(",col,")"))
colnames(RankTab) <- c("Method",rep(c("mean","se"),9))
RankTab
kable(RankTab,format = "latex",booktabs=TRUE) %>% 
  add_header_above(c(" "=1, "AUC"=2,"rMSPE"=2,"bAcc"=2, "AUC"=2,"rMSPE"=2,"bAcc"=2, "AUC"=2,"rMSPE"=2,"bAcc"=2)) %>%
  add_header_above(c(" "=1, "lymphoma"=6,"lymphoma_big"=6,"darwin"=6))

# # # summary table, combine with tribology, prefer this for datasets

n_showm <- length(rank_methods)
mysumdf <- mydf_all %>% filter(Method %in% rank_methods,
                               dataset!="darwin_big",
                               link==ifelse(Method=="TARP"&dataset=="lymphoma_big","cloglog","logit")) %>% 
  group_by(Method,dataset) %>% summarise(mean_AUC=mean(AUC,na.rm=TRUE),se_AUC=sd(AUC,na.rm=TRUE)/sqrt(100),
                                   mean_rMSPE=mean(rMSPE,na.rm=TRUE),se_rMSPE=sd(rMSPE,na.rm=TRUE)/sqrt(100),
                                   mean_bAcc=mean(bAcc,na.rm=TRUE),se_bAcc=sd(bAcc,na.rm=TRUE)/sqrt(100))
SumTab <- mysumdf %>% pivot_wider(names_from = dataset,values_from = c(3:8),names_vary = "slowest")
SumTab[,-1] <- round(SumTab[,-1],3)
SumTab[,1+1:9*2] <- apply(SumTab[,1+1:9*2],2,function(col)paste0("(",col,")"))
SumTab <- cbind(SumTab,readRDS("../saved_results/table_tribology_rMSPE.rds")[,-1])
colnames(SumTab) <- c("Method",rep(c("mean","se"),10))
SumTab
kable(SumTab,format = "latex",booktabs=TRUE) %>% 
  add_header_above(c(" "=1, "AUC"=2,"rMSPE"=2,"bAcc"=2, "AUC"=2,"rMSPE"=2,"bAcc"=2, "AUC"=2,"rMSPE"=2,"bAcc"=2,"rMSPE"=2)) %>%
  add_header_above(c(" "=1, "lymphoma"=6,"lymphoma_big"=6,"darwin"=6,"tribology"=2))


# smaller version with AUC + rMSPE
kable(SumTab[,c(1,20,21,14:17,2:5,8:11)],format = "latex",booktabs=TRUE) %>% 
  add_header_above(c(" "=1, "rMSPE"=2,"AUC"=2,"rMSPE"=2, "AUC"=2,"rMSPE"=2, "AUC"=2,"rMSPE"=2)) %>%
  add_header_above(c(" "=1, "FTIR spectra"=2,"Darwin"=4,"DLBCL"=4,"DLBCL_extended"=4))

# copy output latex code to latex file


kable(t(SumTab[c(4,6,7,10,11),c(1,20,21)]),format = "latex",booktabs=TRUE)
