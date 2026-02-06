##-----------------------------------------------------##
##--------------------Fig.S6&S7---------------------------##
##-----------------------------------------------------##

library("bayesplot")
library("rstanarm")
library("ggplot2")
library(dplyr)
library(ggplot2)
library(tidyr)
require(boot)
library(Rfast2)
library(emmeans)
library(tidyverse)
library(cowplot)
source("forest.plot.R")
library(gridExtra)
theme_set(theme_bw(base_size = 10)+ theme(legend.position = "right"))


# treatment color palette
trts = c('Saline',"Flu shot",
         # 'CDN nanovax', 'CpG nanovax', 
         'MPLA nanovax',  'R848 nanovax')
trt.colorPal =c(rgb(220,105,125,maxColorValue = 255),
                rgb(51,117,56, maxColorValue = 255), 
                # rgb(47, 37, 217,maxColorValue = 255),
                # rgb(112, 112, 112,maxColorValue = 255),
                rgb(230,159,0,maxColorValue = 255),
                rgb(86,180,233,maxColorValue = 255))
trt.shapePal = c(1, 2, 
                 # 15, 17, 
                 4, 18)


# Data read-in
dat0 = read.csv('./data/Aged HET3 mice data -All  January 2024- Ana.csv')
head(dat0)
dim(dat0)
## convert the Mass of dead mouse to NA
dat0$MassDay6 = as.numeric(dat0$MassDay6) 
dat0$MassDay5 = as.numeric(dat0$MassDay5)
dat0 = dat0 %>%
  # dplyr::filter(INCLUDE == 1) %>% # only include the reasonable mice
  mutate(wtLossThroughD4 = (MassDay4-MassDay0)/MassDay0)%>%
  mutate(wtLossThroughD5 = (MassDay5-MassDay0)/MassDay0)%>%
  mutate(wtLossThroughD6 = (MassDay6-MassDay0)/MassDay0)

data = dat0 %>%
  # dplyr::select(Cage:Sex, ViralLoadmL) %>%
  # mutate(logViralLoad = log2(ViralLoadmL + 1)) %>%
  mutate(nonzero = ViralLoadmL>0)
ytemp = data$ViralLoadmL
y = ytemp
y[ytemp>0] = log10(y[ytemp>0])

data$logViralLoad = y


## set Saline to be base level
unique(data$Treatment)
data$Treatment = factor(data$Treatment, levels = trts)


##--------------------FigS7(A)---------------------------##
library(combinat)
library(caTools)
#Permutation test
exact.permutation.test <- function(y1,y2, nperm = NULL, exact = FALSE){
  distribution=c()
  np1 = length(y1)
  np2 = length(y2)
  np = np1 + np2
  obs_diff = mean(y1) - mean(y2)
  y_all = c(y1,y2)
  if(exact){
    # exact test
    ind1 = caTools::combs(1:np, np1)
    n = nrow(ind1)
    for(ii in 1:n){
      mu1 = mean(y_all[ind1[ii,]])
      mu2 = mean(y_all[-ind1[ii,]])
      distribution[ii] = mu1-mu2
    }
  } else {
    n = nperm
    # not exact test, MCMC
    for(ii in 1:n){
      ind1 = sample(1:np, np1,replace = FALSE)
      mu1 = mean(y_all[ind1])
      mu2 = mean(y_all[-ind1])
      distribution[ii] = mu1-mu2
      # distribution[i]=diff(by(outcome, sample(treatment, length(treatment), FALSE), mean))
    }
  }
  
  if(obs_diff>0){
    result.doubletail=min(sum(distribution >= abs(obs_diff))/n *2,1)
  } else{
    result.doubletail=min(sum(distribution <=-abs(obs_diff))/n *2,1)
  }
  result.twosided=min(mean(abs(distribution) >= abs(obs_diff)),1)
  
  return(list(p.value.doubletail = result.doubletail,p.value.twosided = result.twosided, obs_diff = obs_diff, np1 = np1, np2 = np2, n = n, distribution = distribution))
}
trt.list = levels(data$Treatment)
aa = "Aged"


res_all_all = c()
res_all = data.frame()
for(i in 1:(length(trt.list)-1)){
  for(j in (i+1):(length(trt.list))){
    # i = 1; j = 4
    y1 = log10(unlist(data %>% dplyr::filter(Treatment == trt.list[i]) %>%
                        dplyr::select(ViralLoadmL))+1) 
    y2 = log10(unlist(data %>% dplyr::filter(Treatment == trt.list[j]) %>%
                        dplyr::select(ViralLoadmL))+1)
    test1 = exact.permutation.test(y1, y2, exact = TRUE)
    test1a = perm::permTS(y1, y2, alternative = "two.sided", exact = "exact.ce")
    
    test1$p.value.twosided
    test1a$p.value
    res_all = rbind(res_all, c(paste0(aa, collapse = ''), trt.list[i],  trt.list[j], 
                               test1$obs_diff,
                               test1$np1, test1$np2,
                               test1$n, 
                               test1$p.value.doubletail, test1$p.value.twosided,
                               test1a$p.value)
    )
  }
}

colnames(res_all) = c("Age", "Trt1", "Trt2", "Trt1-Trt2", "n1", "n2", "n_perm","p_value_DT", "p_value_TS", "p_value_TS2")
res_all = data.frame(res_all)
res_all = res_all %>% dplyr::mutate(p_adjust_DT = p.adjust(res_all$p_value_DT, method = "fdr")) %>%
  dplyr::mutate(p_adjust_TS = p.adjust(res_all$p_value_TS, method = "fdr")) %>%
  dplyr::mutate(p_adjust_TS2 = p.adjust(res_all$p_value_TS2, method = "fdr"))
write.csv(res_all, file = paste0("res/perm_test_transl_", paste0(aa, collapse = ''),".csv"),row.names = F)
res_all_all = rbind(res_all_all, res_all)



res_all_all$Trt1.Trt2 = round(as.numeric(res_all_all$Trt1.Trt2), digits = 3)
res_all_all$p_adjust_DT = round(as.numeric(res_all_all$p_adjust_DT), digits = 4)
res_all_all$p_adjust_TS = round(as.numeric(res_all_all$p_adjust_TS), digits = 4)
res_all_all$p_adjust_TS2 = round(as.numeric(res_all_all$p_adjust_TS2), digits = 4)
res_all_all$p_value_DT = round(as.numeric(res_all_all$p_value_DT), digits = 4)
res_all_all$p_value_TS = round(as.numeric(res_all_all$p_value_TS), digits = 4)
res_all_all$p_value_TS2 = round(as.numeric(res_all_all$p_value_TS2), digits = 4)


res_all_all


##--------------------Fig.S.6(F to K)---------------------------##

data1 = data %>% mutate(WeightLossDay1 = (MassDay1-MassDay0)/MassDay0, 
                        WeightLossDay2 = (MassDay2-MassDay0)/MassDay0, 
                        WeightLossDay3 = (MassDay3-MassDay0)/MassDay0, 
                        WeightLossDay4 = (MassDay4-MassDay0)/MassDay0, 
                        WeightLossDay5 = (MassDay5-MassDay0)/MassDay0, 
                        WeightLossDay6 = (MassDay6-MassDay0)/MassDay0)

dd = 2
plot_res_all_both = c()

for (dd in 1:6){
  # res2 = lm(paste0('WeightLossDay',dd,' ~  Treatment + MassDay0'), data = data1)
  res2 = lm(paste0('WeightLossDay',dd,'~  Treatment'), data = data1)
  # res2 = lm(paste0('log2(MassDay',dd,')-log2(MassDay0) ~  Treatment'), data = data1)
  
  summary(res2)
  anova(res2)
  res.lsm = lsmeans(res2, pairwise ~ Treatment) # pairwise testing averaged over age & days
  res_means = res.lsm$lsmeans
  plot_res = as.data.frame(res_means)
  plot_res$Day = dd
  plot_res_all_both = rbind(plot_res_all_both, plot_res)
  # p1 = forest.plot(res.lsm, -10,5,title =  paste0("log2 Weight Change through Day", dd,  ": Aged"))
  p1 = forest.plot(res.lsm$contrasts, -0.2,0.1,nudge_x = 0.04,title =  paste0("Percent Weight loss through Day", dd, ": Aged"))
  # png(paste0('res/Weight',dd,'.png'),res = 200,width = 1800, height = 800)
  # print(p1)
  # dev.off()
  print(p1)
  res_means = res.lsm$lsmeans
  plot_res = as.data.frame(res_means)
  plot_res$Day = dd
  plot_res_all_both = rbind(plot_res_all_both, plot_res)
}

##--------------------Fig.S.6(P to U)---------------------------##
data.long = data  %>%
  tidyr::pivot_longer(MassDay0:MassDay6,names_to = "Days", values_to = "Mass") %>%
  mutate(Days = as.integer(gsub("MassDay","",Days))) 
data.long$Days = factor(data.long$Days)
data.long$ID = factor(data.long$ID)
data.long$Treatment = factor(data.long$Treatment)

data.long = data  %>%
  tidyr::pivot_longer(PenHDay0:PenHDay6,names_to = "Days", values_to = "PenH") %>%
  mutate(Days = as.integer(gsub("PenHDay","",Days))) 
data.long$Days = factor(data.long$Days)
data.long$ID = factor(data.long$ID)
data.long$Treatment = factor(data.long$Treatment)
plot_res_all = c()

for(dd in c(1,2,3,4,5,6)){
  # dd = 3
  data.long2 = data.long %>% dplyr::filter(Days == dd)
  # data.long %>% dplyr::filter(!is.na(PenH)) %>% group_by(Days) %>% summarise(n = n())
  # o.un = lm(PenH ~ Age + Treatment + Age:Treatment + MassDay0, data = data.long2)
  o.un = lm(PenH ~ Treatment + MassDay0, data = data.long2)
  
  library(emmeans)
  EMM <- emmeans(o.un,  ~Treatment)   # or emmeans(model, ~ treatment)
  CON <- pairs(EMM, by = "Treatment")
  # test(pairs(EMM, simple = "Treatment"), by = NULL, adjust = "Tukey")    # compare treats for each dose -- "simple effects"
  
  # summary(o.un)
  # getVarCov(o.un)
  # heatmap(getVarCov(o.un), Rowv = NA, Colv = NA, revC =T,symm = T)
  # anova(o.un)
  par(mfrow = c(1,3))
  res.lsm = lsmeans(o.un, pairwise ~ Treatment) # pairwise testing averaged over age & days
  res_means = res.lsm$lsmeans
  plot_res = as.data.frame(res_means)
  plot_res$Day = dd
  plot_res_all_both = rbind(plot_res_all_both, plot_res)
  p1 = forest.plot(res.lsm$contrasts, -1,5,title =  paste0("PenHDay",dd, ": Aged"))
  res.lsm = lsmeans(o.un, pairwise ~ Treatment) # pairwise testing for each age group averaged over days
  # png(paste0('res/PenH',dd,'.png'),res = 200,width = 1800, height = 800)
  print(p1)
  # dev.off()
  # just plot the histograms of objective functions for each age
  res_means = res.lsm$lsmeans
  plot_res = as.data.frame(res_means)
  p2 = ggplot(plot_res, aes(color =Treatment), group = interaction(Treatment)) + 
    geom_point( aes(lsmean,Treatment)) +
    geom_errorbar(aes(y = Treatment, xmin = lower.CL, xmax = upper.CL), width =0.5) +
    theme_bw() +  xlab(paste0('Estimated PenH Day', dd))+
    scale_color_manual(name = "", values = trt.colorPal) + 
    theme(legend.position = "none") 
  # print(p2)
  plot_res = as.data.frame(res_means)
  plot_res$Day = dd
  plot_res_all = rbind(plot_res_all, plot_res)
  
}
