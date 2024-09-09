## read in file and plot  -----------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------

pacman::p_load(dplyr, ggplot2, tidyr, ggrepel,knitr,kableExtra)

resobj <- readRDS("../saved_results/SPARglm_data_trib_nset1_reps100_nmeth20.rds")

res <- resobj$res
methods <- dimnames(res)[[3]]
resobj$dataset_sizes

# RF and SVM independent of flink, so copy to all
res[,,c(8,9),2,1] <- res[,,c(8,9),3,1] <- res[,,c(8,9),1,1]

mydf_all <- data.frame(pivot_longer(data.frame(res[,1,,1,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="rMSPE"),
                       rDev=pivot_longer(data.frame(res[,2,,1,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="rDev")$rDev,
                       rDev_tr=pivot_longer(data.frame(res[,3,,1,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="rDev_tr")$rDev_tr,
                       NumAct=pivot_longer(data.frame(res[,4,,1,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="NumAct")$NumAct,
                       Time=pivot_longer(data.frame(res[,5,,1,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Time")$Time,
                       resobj$dataset_sizes[1,],
                       flink="poisson(log)")
mydf_all <- rbind(mydf_all,
                  data.frame(pivot_longer(data.frame(res[,1,,2,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="rMSPE"),
                             rDev=pivot_longer(data.frame(res[,2,,2,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="rDev")$rDev,
                             rDev_tr=pivot_longer(data.frame(res[,3,,2,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="rDev_tr")$rDev_tr,
                             NumAct=pivot_longer(data.frame(res[,4,,2,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="NumAct")$NumAct,
                             Time=pivot_longer(data.frame(res[,5,,2,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Time")$Time,
                             resobj$dataset_sizes[1,],
                             flink="gaussian(log)"))

mydf_all <- rbind(mydf_all,
                  data.frame(pivot_longer(data.frame(res[,1,,3,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="rMSPE"),
                             rDev=pivot_longer(data.frame(res[,2,,3,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="rDev")$rDev,
                             rDev_tr=pivot_longer(data.frame(res[,3,,3,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="rDev_tr")$rDev_tr,
                             NumAct=pivot_longer(data.frame(res[,4,,3,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="NumAct")$NumAct,
                             Time=pivot_longer(data.frame(res[,5,,3,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Time")$Time,
                             resobj$dataset_sizes[1,],
                             flink="quasipoisson(log)"))

      
# # here we only have one data set      
# for (k in 2:nrow(resobj$dataset_sizes)) {
#   mydf_all <- rbind(mydf_all,
#                     data.frame(pivot_longer(data.frame(res[,1,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="rMSPE"),
#                                rDev=pivot_longer(data.frame(res[,2,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="rDev")$rDev,
#                                rDev_tr=pivot_longer(data.frame(res[,3,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="rDev_tr")$rDev_tr,
#                                NumAct=pivot_longer(data.frame(res[,4,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="NumAct")$NumAct,
#                                Time=pivot_longer(data.frame(res[,5,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Time")$Time,
#                                resobj$dataset_sizes[k,]))
# 
# }

mydf_all$Method <- stringr::str_replace_all(mydf_all$Method,"\\."," ")
mydf_all$Method <- factor(mydf_all$Method,levels = methods)
mydf_all$flink <- factor(mydf_all$flink,levels = c("poisson(log)","gaussian(log)","quasipoisson(log)"))

# mydf_all <- mutate(mydf_all,"isSparse"=ifelse(Method%in%c("AdLASSO","ElNet","SIS"),TRUE,FALSE))
# ,linetype=isSparse
methods
# show_methods <- methods
show_methods <- methods[c(1,3,6,7,8,9,11,14,15,18)]
rank_methods <- methods[c(1,3,4,6,7,8,9,11,14,15,18)]


mydf_all %>% filter(Method %in% show_methods,
                    !(Method=="SIS" & flink=="quasipoisson(log)")) %>%
  ggplot(aes(x=flink,y=rMSPE,fill=Method)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  # ggh4x::facet_grid2(cov_setting~act_setting, scales = "free_y",independent = "y") +
  facet_grid(.~Method, scales = "free_y") +
  theme(legend.position = "none") +
  coord_cartesian(ylim=c(0,0.9))
# ggsave(paste0("../plots/glm_data_trib_flink_rMSPE.pdf"), height = 6, width = 16)
# not much difference, gaussian best for all

mydf_all %>% filter(Method %in% show_methods,
                    !(Method=="SIS" & flink=="quasipoisson(log)")) %>%
  ggplot(aes(x=flink,y=Time,fill=Method)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  # ggh4x::facet_grid2(cov_setting~act_setting, scales = "free_y",independent = "y") +
  facet_grid(.~Method, scales = "free_y") +
  theme(legend.position = "none") +
  scale_y_log10()
# poisson(log) would be faster

mydf_all %>% filter(Method %in% show_methods,
                    flink=="gaussian(log)") %>%
  ggplot(aes(x=Method,y=rMSPE,fill=Method)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  # ggh4x::facet_grid2(cov_setting~act_setting, scales = "free_y",independent = "y") +
  theme(legend.position = "none") +
  coord_cartesian(ylim=c(0,1))
# ggsave(paste0("../plots/glm_data_trib_rMSPE.pdf"), height = 6, width = 10)

mydf_all %>% filter(Method %in% show_methods,
                    flink=="gaussian(log)") %>%
  ggplot(aes(x=Method,y=rDev,fill=Method)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  # ggh4x::facet_grid2(cov_setting~act_setting, scales = "free_y",independent = "y") +
  theme(legend.position = "none") +
  coord_cartesian(ylim=c(0,1))


mydf_all %>% filter(Method %in% show_methods,
                    flink=="gaussian(log)") %>%
  ggplot(aes(x=Method,y=rDev_tr,fill=Method)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  # ggh4x::facet_grid2(cov_setting~act_setting, scales = "free_y",independent = "y") +
  theme(legend.position = "none") +
  coord_cartesian(ylim=c(0,1))


mydf_all %>% filter(Method %in% show_methods,
                    flink=="gaussian(log)") %>%
  ggplot(aes(x=Method,y=NumAct,fill=Method)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  # ggh4x::facet_grid2(cov_setting~act_setting, scales = "free_y",independent = "y") +
  theme(legend.position = "none")


mydf_all %>% filter(Method %in% show_methods,
                    flink=="gaussian(log)") %>%
  ggplot(aes(x=Method,y=Time,fill=Method)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  # ggh4x::facet_grid2(cov_setting~act_setting, scales = "free_y",independent = "y") +
  theme(legend.position = "none") +
  scale_y_log10()


# # rank tables

# each dataset
n_showm <- length(rank_methods)
myrankdf <- mydf_all %>% filter(Method %in% rank_methods,
                                flink=="gaussian(log)") %>% group_by(rep,dataset) %>% 
  mutate(rank_rDev=rank(rDev) ,rank_rMSPE=rank(rMSPE))
myranksumdf <- myrankdf %>% group_by(Method,dataset) %>% summarise(mean_rank_rMSPE = mean(rank_rMSPE), se_rank_rMSPE = sd(rank_rMSPE)/sqrt(100),
                                                                                 mean_rank_rDev = mean(rank_rDev), se_rank_rDev = sd(rank_rDev)/sqrt(100)) 
myranksumdf


# # summarize rMSPE
n_showm <- length(rank_methods)
mysumdf <- mydf_all %>% filter(Method %in% rank_methods,
                                flink=="gaussian(log)") %>% group_by(Method,dataset) %>% 
  summarise(mean_rMSPE=mean(rMSPE,na.rm=TRUE),se_rMSPE=sd(rMSPE,na.rm=TRUE)/sqrt(100))
SumTab <- mysumdf[,-2]
SumTab[,-1] <- round(SumTab[,-1],3)
SumTab[,3] <- apply(SumTab[,3],2,function(col)paste0("(",col,")"))
SumTab
# saveRDS(SumTab,"../saved_results/table_tribology_rMSPE.rds")
