##------------------------------------------------------##
##----------------------Fig.S5--------------------------##
##------------------------------------------------------##

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
trts = c("Saline", "CDN nanovax", "CpG nanovax", "MPLA nanovax", "R848 nanovax")
trt.shapePal = c(1, 15, 17, 4, 18)
# Saline- Red 220, green 105, blue 125
# CDN nanovax-Red 47, green 37, blue 217
# CpG nanovax- Red 112, green 112, blue 112
# MPLA- nanovax Red 230, green 159, blue 0
# R848 nanovaz- Red 86, Green 180, blue, 233
# Flu shot- Red 51, green 117, blue 56.
trt.colorPal =c(rgb(220,105,125,maxColorValue = 255), 
                rgb(47, 37, 217,maxColorValue = 255),
                rgb(112, 112, 112,maxColorValue = 255),
                rgb(230,159,0,maxColorValue = 255),
                rgb(86,180,233,maxColorValue = 255))
# rgb(51,117,56, maxColorValue = 255)

#### Data read-in
dat0 = read.csv('./data/1st invivo study young vs aged all data tables.csv')
head(dat0)
dim(dat0)
## convert the Mass of dead mouse to NA
dat0$MassDay6 = as.numeric(dat0$MassDay6) 
dat0 = dat0 %>%
  # dplyr::filter(INCLUDE == 1) %>% # only include the reasonable mice
  mutate(wtLossThroughD4 = (MassDay4-MassDay0)/MassDay0)%>%
  mutate(wtLossThroughD5 = (MassDay5-MassDay0)/MassDay0)%>%
  mutate(wtLossThroughD6 = (MassDay6-MassDay0)/MassDay0)

data = dat0 %>%
  # dplyr::select(Cage:Sex, ViralLoadmL) %>%
  mutate(logViralLoad = log2(ViralLoadmL + 1)) %>%
  mutate(nonzero = ViralLoadmL>0)

## set Saline to be base level
unique(data$Treatment)
data$Treatment = factor(data$Treatment, levels = trts)
data_noncontrol = data %>% dplyr::filter(!Treatment %in% "Saline")
##----------------------Fig.S5(C,D,E,F,G,H)--------------------------##
data1 = data %>% mutate(WeightLossDay1 = (MassDay1-MassDay0)/MassDay0, 
                        WeightLossDay2 = (MassDay2-MassDay0)/MassDay0, 
                        WeightLossDay3 = (MassDay3-MassDay0)/MassDay0, 
                        WeightLossDay4 = (MassDay4-MassDay0)/MassDay0, 
                        WeightLossDay5 = (MassDay5-MassDay0)/MassDay0, 
                        WeightLossDay6 = (MassDay6-MassDay0)/MassDay0)

plot_res_all = c()
plot_res_all_both = c()


for (dd in 1:6){
  res2 = lm(paste0('PenHDay',dd,' ~ Age + Treatment + Age:Treatment +  MassDay0'), data = data1)
  # summary(res2)
  # anova(res2)
  res.lsm = lsmeans(res2, pairwise ~ Treatment) # pairwise testing averaged over age & days
  p1 = forest.plot(res.lsm$contrasts, min_val = -4,max_val = 12,nudge_x = 2,title =  paste0("PenHDay", dd, ": Aged & Young"))
  res_means = res.lsm$lsmeans
  plot_res = as.data.frame(res_means)
  plot_res$Day = dd
  plot_res_all_both = rbind(plot_res_all_both, plot_res)
  res.lsm = lsmeans(res2, pairwise ~ Treatment|Age) # pairwise testing for each age group averaged over days
  # save(res.lsm, file = "res.lsm.Penh.rda")
  # load(file = "res.lsm.Penh.rda")
  typeof(res.lsm$contrasts)
  contrast = res.lsm$contrasts 
  res_aged = contrast[1:10]
  res_young = contrast[11:20]
  p2 = forest.plot(res_aged, min_val = -4,max_val = 12,nudge_x = 2,title = paste0("PenHDay", dd, ": Aged"))
  p3 = forest.plot(res_young, min_val = -4,max_val = 12,nudge_x = 3,title =  paste0("PenHDay", dd, ": Young"))
  require(gridExtra)
  # grid.arrange(p1, p2, p3,  ncol = 1)
  # png(paste0('res/PP(enH',dd,'.png'),res = 110,width = 900, height = 1200)
  # grid.arrange(p1, p2, p3,  ncol = 1)
  # dev.off()
  res_means = res.lsm$lsmeans
  plot_res = as.data.frame(res_means)
  plot_res$Day = dd
  plot_res_all = rbind(plot_res_all, plot_res)
  p4 = ggplot(plot_res, aes(color =Treatment), group = interaction(Age, Treatment)) + 
    geom_point( aes(lsmean,Treatment)) +
    geom_errorbar(aes(y = Treatment, xmin = lower.CL, xmax = upper.CL), width =0.5) +
    theme_bw(base_size=15) +  xlab(paste0('Estimated PenH Day', dd))+
    scale_color_manual(name = "", values = trt.colorPal) + 
    theme(legend.position = "none") +
    facet_wrap(~Age, ncol = 1) + 
    xlim(-5,10)
  # print(p4)
  grid.arrange( p2, p3 , ncol = 1)
}