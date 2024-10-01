pacman::p_load(tidyr,dplyr,ggplot2)
res <- readRDS("../saved_results/Result_ScrRP_Experiment_all.rds")

methods <- dimnames(res)[[3]]

nmethods <- length(methods)
settings <- attributes(res)$settings
settings$family <- sapply(settings$family,function(obj) paste0(obj$family,"(",obj$link,")"))

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
                               settings[k,],
                               setting=k)
  )
}


mydf_all$Method <- factor(mydf_all$Method,levels = methods[c(1:5,7:12,6,13:nmethods)])

levels(mydf_all$Method)[12] <- "Ridge_limit0"
methods[6] <-  "Ridge_limit0"

mydf_all$act_setting <- factor(mydf_all$act_setting,levels = c("sparse","medium","dense"))


tmp_df <- mydf_all %>% filter(MSPE<150) %>%
  pivot_longer(c(MSPE,MSLE),names_to = "Measure",values_to = "value")
tmp_df$Measure <- factor(tmp_df$Measure,levels = c("MSLE","MSPE"))
tmp_df %>%  filter(Method %in% methods,
                   p==2000, 
                   act_setting=="medium") %>%
  ggplot(aes(x=Method,y=value,fill=Method)) +
  geom_boxplot() +
  theme(legend.position = "none") +
  ggh4x::facet_grid2(Measure~family, scales = "free_y",independent = "y") +
  # facet_grid(Measure~family,scales = "free_y")+
  labs(y=" ") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# ggsave("../plots/RPexperiment_families_all.pdf", height = 5, width = 8)


tmp_df %>%  filter(Method %in% methods,
                   p==500, 
                   act_setting=="medium") %>%
  ggplot(aes(x=Method,y=value,fill=Method)) +
  geom_boxplot() +
  theme(legend.position = "none") +
  ggh4x::facet_grid2(Measure~family, scales = "free_y",independent = "y") +
  # facet_grid(Measure~family,scales = "free_y")+
  labs(y=" ") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


tmp_df %>%  filter(Method %in% methods,
                   p==2000, 
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
                   p==2000, 
                   act_setting=="medium") %>%
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
                               family!="binomial(cloglog)",
                               p==2000) %>%
  group_by(Method,family) %>% summarize(mCor=mean(Cor_scr,na.rm=TRUE),seCor = sd(Cor_scr,na.rm=TRUE)/sqrt(100))
corTab <- corTab %>% pivot_wider(names_from = family,values_from = c(mCor,seCor),names_vary = "slowest")
corTab[,-1] <- round(corTab[,-1],3)
corTab[,1+1:4*2] <- apply(corTab[,1+1:4*2],2,function(col)paste0("(",col,")"))
kable(corTab,format = "latex",booktabs=TRUE) %>%
  add_header_above(c("Method"=1, "binomial(logit)"=2,"gaussian(identity)"=2,"gaussian(log)"=2, "poisson(log)"=2))




mydf_all %>%  filter(Method %in% methods[-c(1:5,8,nmethods)],
                     act_setting=="medium",
                     family!="binomial(cloglog)",
                     p==2000) %>%
  ggplot(aes(x=Method,y=Cor_scr,fill=Method)) +
  geom_boxplot() +
  theme(legend.position = "none") +
  facet_grid(.~family,scales = "free_y")+
  labs(y=" ") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_x_discrete(labels= c("L2_cv","L2_10","L2_1","L2_01","L2_001","L2_limit0","L2_dev08","L2_dev095","L2_dev0999","marGLM"))
# ggsave("../plots/ScrCor_families_p2000med.pdf", height = 3, width = 8)

# # average lambdas
lam_tab <- round(t(rbind(apply(res[,c(7),c(7,13:15),3],2,function(x) c(mean(x),sd(x)/sqrt(100))),
                 apply(res[,c(7),c(7,13:15),1],2,function(x) c(mean(x),sd(x)/sqrt(100))),
                 apply(res[,c(7),c(7,13:15),2],2,function(x) c(mean(x),sd(x)/sqrt(100))),
                 apply(res[,c(7),c(7,13:15),5],2,function(x) c(mean(x),sd(x)/sqrt(100))))),3)
lam_tab[,1:4*2] <- apply(lam_tab[,1:4*2],2,function(col)paste0("(",col,")"))
rownames(lam_tab) <- c("L2_cv","L2_dev08","L2_dev095","L2_dev0999")
kable(lam_tab,format = "latex",booktabs=TRUE) %>%
  add_header_above(c("Method"=1, "binomial(logit)"=2,"gaussian(identity)"=2,"gaussian(log)"=2, "poisson(log)"=2))

# 1 - AUC for binomial
mydf_all %>%  filter(Method %in% methods[-c(4,8:12)],
                               act_setting=="medium",
                               p==2000,MSPE<75,
                               family!="binomial(cloglog)") %>% 
  mutate(pred_error=ifelse(family=="binomial(logit)",1-AUC,MSPE))  %>%
  pivot_longer(c(MSLE,pred_error),names_to = "Measure",values_to = "value") %>%
  ggplot(aes(x=Method,y=value,fill=Method)) +
  geom_boxplot() +
  theme(legend.position = "none") +
  ggh4x::facet_grid2(Measure~family, scales = "free_y",independent = "y") +
  labs(y=" ") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_x_discrete(labels= c("ElNet","AdLASSO","Gaussian","SparseCW","L2_cv","L2_limit0","L2_dev08","L2_dev095","L2_dev0999","marGLM","True_Beta"))
# ggsave("../plots/RP_families_p2000med.pdf", height = 4, width = 10)

