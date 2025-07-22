pacman::p_load(tidyr,dplyr,ggplot2)
res <- readRDS("../saved_results/Result_ScrRP_Experiment_25.rds")

methods <- dimnames(res)[[3]]

nmethods <- length(methods)
settings <- attributes(res)$settings
settings$family <- sapply(settings$family,function(obj) paste0(obj$family,"(",obj$link,")"))
settings$family[settings$weight_binom>1] <- sapply(settings$family[settings$weight_binom>1],function(tmp) paste0(tmp,"n100"))


mydf_all <- data.frame(pivot_longer(data.frame(res[,1,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="AUC"),
                       MSPE=pivot_longer(data.frame(res[,2,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="MSPE")$MSPE,
                       MSLE=pivot_longer(data.frame(res[,3,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="MSLE")$MSLE,
                       pAUC=pivot_longer(data.frame(res[,4,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="pAUC")$pAUC,
                       Cor_scr=pivot_longer(data.frame(res[,5,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Cor_scr")$Cor_scr,
                       pAUC_scr=pivot_longer(data.frame(res[,6,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="pAUC_scr")$pAUC_scr,
                       lambda=pivot_longer(data.frame(res[,7,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="lambda")$lambda,
                       dev.ratio=pivot_longer(data.frame(res[,8,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="dev.ratio")$dev.ratio,
                       mtr_over_mf_abs=pivot_longer(data.frame(res[,9,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="mtr_over_mf_abs")$mtr_over_mf_abs,
                       TrRatio_3a=pivot_longer(data.frame(res[,10,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="TrRatio_3a")$TrRatio_3a,
                       check_gy=pivot_longer(data.frame(res[,11,,1],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="check_gy")$check_gy,
                       settings[1,],
                       setting=1)
for (k in 2:nrow(settings)) {
  mydf_all <- rbind(mydf_all,
                    data.frame(pivot_longer(data.frame(res[,1,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="AUC"),
                               MSPE=pivot_longer(data.frame(res[,2,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="MSPE")$MSPE,
                               MSLE=pivot_longer(data.frame(res[,3,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="MSLE")$MSLE,
                               pAUC=pivot_longer(data.frame(res[,4,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="pAUC")$pAUC,
                               Cor_scr=pivot_longer(data.frame(res[,5,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="Cor_scr")$Cor_scr,
                               pAUC_scr=pivot_longer(data.frame(res[,6,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="pAUC_scr")$pAUC_scr,
                               lambda=pivot_longer(data.frame(res[,7,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="lambda")$lambda,
                               dev.ratio=pivot_longer(data.frame(res[,8,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="dev.ratio")$dev.ratio,
                               mtr_over_mf_abs=pivot_longer(data.frame(res[,9,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="mtr_over_mf_abs")$mtr_over_mf_abs,
                               TrRatio_3a=pivot_longer(data.frame(res[,10,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="TrRatio_3a")$TrRatio_3a,
                               check_gy=pivot_longer(data.frame(res[,11,,k],rep=1:100),1:(dim(res)[3]),names_to="Method",values_to="check_gy")$check_gy,
                               settings[k,],
                               setting=k)
  )
}


mydf_all$Method <- factor(mydf_all$Method,levels = methods[c(1:5,7:12,6,13:nmethods)])

levels(mydf_all$Method)[12] <- "L2_limit0"
methods[6] <-  "L2_limit0"

mydf_all$act_setting <- factor(mydf_all$act_setting,levels = c("sparse","medium","dense"))

mydf_all <- mutate(mydf_all,signal_level = ifelse(family=="poisson(log)",
                                                  ifelse(signal_strength==1/4,"high","low"),
                                                  ifelse(family=="gaussian(identity)",
                                                         ifelse(signal_strength==10,"high","low"),
                                                         ifelse(signal_strength %in% c(1/8,100,1,5,1000),"high","low"))),
                   pred_error=ifelse(family=="binomial(logit)",1-AUC,MSPE))

tmp_df <- mydf_all %>% filter(MSPE<150) %>%
  pivot_longer(c(MSPE,MSLE),names_to = "Measure",values_to = "value")
tmp_df$Measure <- factor(tmp_df$Measure,levels = c("MSLE","MSPE"))
tmp_df %>%  filter(Method %in% methods,
                   p==2000, 
                   m==50,
                   signal_level=="low",
                   act_setting=="medium") %>%
  ggplot(aes(x=Method,y=value,fill=Method)) +
  geom_boxplot() +
  theme(legend.position = "none") +
  ggh4x::facet_grid2(Measure~family, scales = "free_y",independent = "y") +
  # facet_grid(Measure~family,scales = "free_y")+
  labs(y=" ") +
  theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust=1))
# ggsave("../plots/RPexperiment_families_all.pdf", height = 5, width = 8)

tmp_df <- mydf_all %>% 
  filter(MSPE<105) %>%
  filter(Method %in% methods,
                   p==2000, 
         m==50,
         signal_level=="high",
                   act_setting=="medium") 
tmp_df %>%
  ggplot(aes(x=Method,y=Cor_scr,fill=Method)) +
  geom_boxplot() +
  theme(legend.position = "none") +
  ggh4x::facet_grid2(check_gy~family, scales = "free_y",independent = "y") +
  # labs(y=" ") +
  theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust=1))
# ggsave("../plots/RPexperiment_families_gy.pdf", height = 10, width = 20)

tmp_df <- mydf_all %>% 
  filter(MSPE<105) %>%
  filter(Method %in% methods,
         p==2000, 
         m==50,
         signal_level=="high",
         act_setting=="sparse") 
tmp_df %>%
  ggplot(aes(x=Method,y=Cor_scr,fill=Method)) +
  geom_boxplot() +
  theme(legend.position = "none") +
  ggh4x::facet_grid2(check_gy~family, scales = "free_y",independent = "y") +
  labs(y=" ") +
  theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust=1))
# ggsave("../plots/RPexperiment_families_cor_sparse.pdf", height = 10, width = 20)

tmp_df %>% filter(startsWith(family,"binom")) %>%
  ggplot(aes(x=Method,y=AUC,fill=Method)) +
  geom_boxplot() +
  theme(legend.position = "none") +
  ggh4x::facet_grid2(check_gy~family, scales = "free_y",independent = "y") +
  labs(y=" ") +
  theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust=1))
# ggsave("../plots/RPexperiment_families_AUC.pdf", height = 10, width = 20)

tmp_df <- mydf_all %>% 
  filter(Method %in% methods,
         MSPE<105,
         m==50,
         p==2000, 
         act_setting=="medium")
tmp_df %>% 
  ggplot(aes(x=Method,y=MSPE,fill=Method)) +
  geom_boxplot() +
  theme(legend.position = "none") +
  ggh4x::facet_grid2(signal_level~family, scales = "free_y",independent = "y") +
  labs(y=" ") +
  theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust=1))
# ggsave("../plots/RPexperiment_families_signal.pdf", height = 10, width = 20)


tmp_df <- mydf_all %>% 
  filter(Method %in% methods,
         p==2000, 
         m==50,
         act_setting=="medium")
tmp_df %>% 
  ggplot(aes(x=Method,y=Cor_scr,fill=Method)) +
  geom_boxplot() +
  theme(legend.position = "none") +
  facet_grid(signal_level~family, scales = "free_y") +
  labs(y=" ") +
  theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust=1))
# ggsave("../plots/RPexperiment_families_signal_cor.pdf", height = 10, width = 20)



tmp_df <- mydf_all %>% 
  pivot_longer(c(Cor_scr,pAUC_scr,mtr_over_mf_abs,TrRatio_3a),names_to = "Measure",values_to = "value")
tmp_df$Measure <- factor(tmp_df$Measure,levels = c("Cor_scr","pAUC_scr","mtr_over_mf_abs","TrRatio_3a"))
tmp_df %>%  filter(Method %in% methods[-c(1:5,nmethods)],
                   p==2000, 
                   act_setting=="medium",
                   m==50,
                   signal_level=="high") %>%
  ggplot(aes(x=Method,y=value,fill=Method)) +
  geom_boxplot() +
  theme(legend.position = "none") +
  ggh4x::facet_grid2(Measure~family, scales = "free_y",independent = "y") +
  # facet_grid(Measure~family,scales = "free_y")+
  labs(y=" ") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# ggsave("../plots/Screxperiment_families_all.pdf", height = 5, width = 8)

tmp_df <- mydf_all %>% 
  pivot_longer(c(Cor_scr,pAUC_scr,mtr_over_mf_abs,TrRatio_3a),names_to = "Measure",values_to = "value")
tmp_df$Measure <- factor(tmp_df$Measure,levels = c("Cor_scr","pAUC_scr","mtr_over_mf_abs","TrRatio_3a"))
tmp_df %>%  filter(Method %in% methods[-c(1:5,nmethods)],
                   p==2000, 
                   signal_level=="high",
                   m==50,
                   act_setting=="sparse") %>%
  ggplot(aes(x=Method,y=value,fill=Method)) +
  geom_boxplot() +
  theme(legend.position = "none") +
  ggh4x::facet_grid2(Measure~family, scales = "free_y",independent = "y") +
  # facet_grid(Measure~family,scales = "free_y")+
  labs(y=" ") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


tmp_df <- mydf_all %>% 
  pivot_longer(c(Cor_scr,pAUC_scr,mtr_over_mf_abs,TrRatio_3a),names_to = "Measure",values_to = "value")
tmp_df$Measure <- factor(tmp_df$Measure,levels = c("Cor_scr","pAUC_scr","mtr_over_mf_abs","TrRatio_3a"))
tmp_df %>%  filter(Method %in% methods[-c(1:5,nmethods)],
                   p==500, 
                   m==50,
                   signal_level=="high",
                   act_setting=="medium") %>%
  ggplot(aes(x=Method,y=value,fill=Method)) +
  geom_boxplot() +
  theme(legend.position = "none") +
  ggh4x::facet_grid2(Measure~family, scales = "free_y",independent = "y") +
  # facet_grid(Measure~family,scales = "free_y")+
  labs(y=" ") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


# table
pacman::p_load(knitr,kableExtra)

corTab <- mydf_all %>%  filter(Method %in% methods[-c(1:5,8,nmethods)],
                               act_setting=="medium",
                               signal_level=="high",
                               family!="binomial(cloglog)",
                               m==50,
                               p==2000) %>%
  group_by(Method,family) %>% summarize(mCor=mean(Cor_scr,na.rm=TRUE),seCor = sd(Cor_scr,na.rm=TRUE)/sqrt(100))
corTab <- corTab %>% pivot_wider(names_from = family,values_from = c(mCor,seCor),names_vary = "slowest")
corTab[,-1] <- round(corTab[,-1],3)
corTab[,1+1:4*2] <- apply(corTab[,1+1:4*2],2,function(col)paste0("(",col,")"))
kable(corTab,format = "latex",booktabs=TRUE) %>%
  add_header_above(c("Method"=1, "binomial(logit)"=2,"gaussian(identity)"=2,"gaussian(log)"=2, "poisson(log)"=2))


mydf_all %>%  filter(Method %in% methods[-c(1:5,8:12,length(methods))],
                     act_setting=="medium",
                     family!="binomial(cloglog)",
                     family!="binomial(logit_th)",
                     family!="binomial(logit)n100",
                     signal_level=="high",
                     m==50,
                     p==2000) %>%
  ggplot(aes(x=Method,y=Cor_scr,fill=Method)) +
  geom_boxplot() +
  theme(legend.position = "none") +
  facet_grid(.~family,scales = "free_y")+
  labs(y=" ") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# ggsave("../plots/ScrCor_families_p2000med.pdf", height = 3, width = 8)

# # average lambdas
lam_tab <- round(t(rbind(apply(res[,c(7),c(7,1+12:15),3],2,function(x) c(mean(x),sd(x)/sqrt(100))),
                 apply(res[,c(7),c(7,1+12:15),1],2,function(x) c(mean(x),sd(x)/sqrt(100))),
                 apply(res[,c(7),c(7,1+12:15),2],2,function(x) c(mean(x),sd(x)/sqrt(100))),
                 apply(res[,c(7),c(7,1+12:15),7],2,function(x) c(mean(x),sd(x)/sqrt(100))))),3)
lam_tab[,1:4*2] <- apply(lam_tab[,1:4*2],2,function(col)paste0("(",col,")"))
rownames(lam_tab) <- c("L2_cv","L2_dev06","L2_dev08","L2_dev095","L2_dev0999")
# settings$family[c(3,1,2,7)]
kable(lam_tab,format = "latex",booktabs=TRUE) %>%
  add_header_above(c("Method"=1, "binomial(logit)"=2,"gaussian(identity)"=2,"gaussian(log)"=2, "poisson(log)"=2))

mydf_all %>%  filter(Method %in% methods[-c(4,8:12)],
                               act_setting=="medium",
                               p==2000,
                     m==50,
                     MSPE<75,
                     signal_level=="high",
                     family!="binomial(cloglog)",
                     family!="binomial(logit_th)",
                     family!="binomial(logit)n100") %>% 
  pivot_longer(c(MSLE,pred_error),names_to = "Measure",values_to = "value") %>%
  ggplot(aes(x=Method,y=value,fill=Method)) +
  geom_boxplot() +
  theme(legend.position = "none") +
  ggh4x::facet_grid2(Measure~family, scales = "free_y",independent = "y") +
  labs(y=" ") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
# ggsave("../plots/RP_families_p2000med.pdf", height = 4, width = 10)




# # # # effect of m
tmp_df <- mydf_all %>% 
  filter(Method %in% methods[c(7,13:17)],
         family %in% c("binomial(logit)","poisson(log)","gaussian(log)","gaussian(identity)"),
         p==2000,
         act_setting=="medium",
         signal_level=="low") %>%
  group_by(Method,family,m) %>% 
  summarise(meanPredError = mean(pred_error),
            sePredError = sd(pred_error)/sqrt(100))


tmp_df %>% 
  ggplot(aes(x=m,y=meanPredError,col=Method,shape=Method)) +
  geom_line() +
  geom_point() +
  # geom_ribbon(aes(ymin=meanMSPE-qnorm(0.975)*seMSPE,ymax=meanMSPE+qnorm(0.975)*seMSPE),
              # alpha=0.2,linetype=2,linewidth=0.1) +
  facet_grid(family~., scales = "free_y") +
  labs(y="prediction error") +
  scale_x_continuous(breaks = c(8,23,50,100,150),
                     labels = c("log(p)","3log(p)","n/4","n/2","3n/4"))
# ggsave("../plots/RPexperiment_effect_m.pdf", height = 5, width = 8)

# # # # new proposed changes:
# remove specific lambda, show MSPE/MSLE for sparse med dense, high low signal

# first cor
tmp_df <- mydf_all %>% 
  filter(Method %in% methods[c(6,7,13,14,15,16,17)],
         family!="binomial(cloglog)",
         family!="binomial(logit_th)",
         family!="binomial(logit)n100",
         m==50,
         p==2000,
         signal_level=="high"
         )

tmp_df %>% 
  ggplot(aes(x=Method,y=Cor_scr,fill=Method)) +
  geom_boxplot() +
  theme(legend.position = "none") +
  ggh4x::facet_grid2(family~act_setting, scales = "free_y",independent = "y") +
  labs(y="Correlation") +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))
# ggsave("../plots/ScrRPExp_act_cor.pdf", height = 6, width = 10)

# now MSLE and MSPE
tmp_df <- mydf_all %>% 
  filter(Method %in% methods[c(1,2,3,5,6,7,13,14,15,16,17)],
         family!="binomial(cloglog)",
         family!="binomial(logit_th)",
         family!="binomial(logit)n100",
         m==50,
         p==2000,
         signal_level=="high",
         MSPE<100
  )

tmp_df %>% 
  ggplot(aes(x=Method,y=MSLE,fill=Method)) +
  geom_boxplot() +
  theme(legend.position = "none") +
  ggh4x::facet_grid2(family~act_setting, scales = "free_y",independent = "y") +
  # labs(y=" ") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# ggsave("../plots/ScrRPExp_act_MSLE.pdf", height = 6, width = 10)

tmp_df %>% 
  ggplot(aes(x=Method,y=pred_error,fill=Method)) +
  geom_boxplot() +
  theme(legend.position = "none") +
  ggh4x::facet_grid2(family~act_setting, scales = "free_y",independent = "y") +
  labs(y="prediction error") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# ggsave("../plots/ScrRPExp_act_pred_error.pdf", height = 6, width = 10)


# now signal strength
tmp_df <- mydf_all %>% 
  filter(Method %in% methods[c(1,2,3,5,6,7,13,14,15,16,17)],
         family!="binomial(cloglog)",
         family!="binomial(logit_th)",
         family!="binomial(logit)n100",
         m==50,
         p==2000,
         act_setting=="medium",
         MSPE < 100
  )
tmp_df %>% 
  pivot_longer(c(Cor_scr,pred_error,MSLE),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=Method,y=value,fill=Method)) +
  geom_boxplot() +
  theme(legend.position = "none") +
  ggh4x::facet_grid2(family~type+signal_level, scales = "free_y",independent = "y") +
  labs(y=" ") +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
# ggsave("../plots/ScrRPExp_signal_all.pdf", height = 7, width = 12)
