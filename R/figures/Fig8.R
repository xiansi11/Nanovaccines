##-----------------------------------------------------##
##--------------------Fig.8---------------------------##
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
  mutate(logViral = log10(ViralLoadmL + 1)) %>%
  mutate(nonzero = ViralLoadmL>0)

## set Saline to be base level
unique(data$Treatment)
data$Treatment = factor(data$Treatment, levels = trts)

source("my_ggbiplot.R")
library(ggrepel) # to avoid overlapping text
trt.Pal = as.data.frame(cbind(trts, trt.colorPal, trt.shapePal))
arrow.length = 1.2
label.adjust = 1.5
## convert the Mass of dead mouse to NA
# Data read-in
dat0 = read.csv('./data/Aged HET3 mice data -All  January 2024- Ana.csv')
# head(dat0)
# dim(dat0)
## convert the Mass of dead mouse to NA
dat0$MassDay6 = as.numeric(dat0$MassDay6) 
dat0$MassDay5 = as.numeric(dat0$MassDay5)
data = dat0 %>%
  mutate(sum.score = Bronchiolitis.bronchitis+ 
           Eptihelial.necrosis+ 
           Interstitial.inflammation+
           BALT.Hyperplasia+
           Hemorrhage+
           Edema)%>%
  mutate(logViral = log10(ViralLoadmL + 1),
         logAbs = log10(Antibodies)) 
data$Treatment = factor(data$Treatment, levels = trts)


dat.all = data %>% dplyr::select(ID, Treatment, sum.score, logViral, logAbs, PenHDay2, PenHDay3, PenHDay4)
headers = c("ID", "Treatment")
## check if log scale is needed
# par(mfrow = c(2,4))
# apply(dat.all[,-(1:3)],2,hist)
dat.all = dat.all %>% 
  mutate(PenH234 = 1/3*(PenHDay2 + PenHDay3 + PenHDay4)) %>%
  mutate(Histopath = sum.score) %>%
  dplyr::select(-PenHDay2, -PenHDay3, -PenHDay4, -sum.score)

scale = T

if(scale == F){
  dat.all.uni = cbind(dat.all %>% dplyr::select(all_of(headers)),
                      apply(dat.all %>% dplyr::select(-all_of(headers)),2,function(x){(x-min(x))/(max(x)-min(x))})
  )
  dat.all = dat.all.uni
} 
##--------------------Fig.8(J)---------------------------##

temp = dat.all
temp.Pal = trt.Pal %>% dplyr::filter(trts %in% trts)

# Check point for raw dataset
temp %>% dplyr::select(ID)
# apply(temp,2,function(x)sum(!is.na(x)))
# table(temp$Treatment)

# Note for each mice, we only have 1 treatment, Select only the mice with required responses
temp = as.data.frame(temp[,apply(temp, 2, function(x){sum(!is.na(x))})>0]) %>%
  tidyr::drop_na()


if(scale == F){
  res.pca <- prcomp(temp %>% 
                      dplyr::select(-all_of(headers)),
                    center = T,scale. = F)
} else{
  res.pca <- prcomp(temp %>% 
                      dplyr::select(-all_of(headers)),
                    center = T,scale. = T)
}



##### PLOTTING
# rotate the loading and scores if loading vector is negative
# ind.pc1 = which(res.pca$rotation[,1] <0);ind.pc1
# if(length(ind.pc1)>0){  
#   res.pca$rotation[ind.pc1,1] = - res.pca$rotation[ind.pc1,1]
#   res.pca$x[ind.pc1] = -res.pca$x[ind.pc1]
# }
if(mean(res.pca$rotation[,1])<0){
  res.pca$rotation[,1] = - res.pca$rotation[,1]
  res.pca$x[,1] = -res.pca$x[,1]
}

if(mean(res.pca$rotation[,2])<0){
  res.pca$rotation[,2] = - res.pca$rotation[,2]
  res.pca$x[,2] = -res.pca$x[,2]
}



##### DISTANCE TABLE
temp$PC1 = res.pca$x[,1]
temp$PC2 = res.pca$x[,2]

out = temp %>% dplyr::group_by(Treatment) %>% 
  dplyr::summarize(meanPC1 = mean(PC1), sdPC1 = sd(PC1),
                   meanPC2 = mean(PC2), sdPC2 = sd(PC2),
                   n = n())
res0 = out %>% dplyr::filter(Treatment == trts[1])
out2 = out %>% 
  mutate(disPC1 = meanPC1 - mean(res0$meanPC1),
         # sddisPC1 = sqrt(sdPC1^2/n + (res0$sdPC1)^2/(res0$n)),
         disPC2 = meanPC2 - mean(res0$meanPC2)
         # sddisPC2 = sqrt(sdPC2^2/n + (res0$sdPC2)^2/(res0$n))
  ) %>%
  mutate(disEucli = sqrt(disPC1^2 + disPC2^2)) 
# %>%
# dplyr::filter(Treatment !=  trts[1])
out2 = out2[order(out2$disPC1),]

age = 'Aged'
write.csv(out2 %>% dplyr::select(Treatment,meanPC1,sdPC1,n,disPC1), file = paste0("res/PCA_combine_",paste0(age,collapse = ""),".csv"))
out3 = temp
write.csv(out3, file = paste0("res/PCA_combine_detail_",paste0(age,collapse = ""), ".csv"))

knitr::kable(out2 %>% dplyr::select(Treatment, meanPC1, sdPC1, n, disPC1))


##--------------------Fig.8(I)---------------------------##
cc = summary(res.pca)
cc$importance
## By Treatment
# biplot(res.pca)
fig1 = my_ggbiplot(res.pca, ellipse=T, 
                   # labels=temp$ID,
                   groups = factor(temp$Treatment, levels = trts),
                   # obs.scale = 2,
                   # var.scale = 2,
                   labels.size = 4, 
                   varname.size = 4,
                   # varname.adjust = 1.2,
                   varname.abbrev = F,
                   alpha=0,
                   arrow.length = arrow.length,
                   label.adjust = label.adjust
)  +
  geom_point(aes(shape=factor(temp$Treatment, levels = trts),
                 colour=factor(temp$Treatment, levels = trts))
  ) +
  scale_shape_manual(values=as.numeric(temp.Pal$trt.shapePal))+
  scale_color_manual(values=(temp.Pal$trt.colorPal))+
  theme_bw(base_size=10) + theme(legend.position = "right")+ 
  labs(color = "Adjuvants", shape = "Adjuvants", 
       x = paste0("PC1 (", round(cc$importance[2,1],3)* 100, "%)"),
       y = paste0("PC2 (", round(cc$importance[2,2],3)* 100, "%)"))+ 
  coord_equal(ratio = 0.7) +
  ggtitle(paste0(age, collapse = ''))
fig1$layers <- c(fig1$layers, fig1$layers[[2]]) # , fig1$layers[[1]]
fig1


##--------------------Fig.8(L)---------------------------##
# linear equal weight
dat0 = read.csv('./data/Aged HET3 mice data -All  January 2024- Ana.csv')
data = read.csv("res/PCA_combine_detail_Aged.csv")
data = data %>% mutate(MassDay0 = dat0$MassDay0)
data$Treatment[data$Treatment == 'Control (no treatment)'] = "Saline"
data$Treatment = factor(data$Treatment, levels = c("Saline","CDN nanovax" , "CpG nanovax", "Flu shot",
                                                   "MPLA nanovax","R848 nanovax"))
# summary(data$logViral)


obj1 = function(logViral, logAbs, PenH234, Histopath){
  1/3 * (logViral + PenH234 + Histopath) - logAbs 
}


obj2 = function(logViral, logAbs, PenH234, Histopath){
  exp(1/3 * (logViral + PenH234 + Histopath) - logAbs)
}


obj3 = function(logViral, logAbs, PenH234, Histopath){
  (1/3 * (logViral + PenH234 + Histopath) - logAbs)^3
  # (1/3 * (logViral + PenH234 + Histopath) - logAbs)^5
}

obj4 = function(logViral, logAbs, PenH234, Histopath){
  (logViral + PenH234 + Histopath) - logAbs 
}

obj5 = function(logViral, logAbs, PenH234, Histopath){
  exp( (logViral + PenH234 + Histopath) - logAbs )
}


obj6 = function(logViral, logAbs, PenH234, Histopath){
  ( (logViral + PenH234 + Histopath) - logAbs )^3
}



datanew = data %>% 
  mutate(Y1 = obj1(logViral, logAbs, PenH234, Histopath),
         Y2 = obj2(logViral, logAbs, PenH234, Histopath),
         Y3 = obj3(logViral, logAbs, PenH234, Histopath),
         Y4 = obj4(logViral, logAbs, PenH234, Histopath),
         Y5 = obj5(logViral, logAbs, PenH234, Histopath),
         Y6 = obj6(logViral, logAbs, PenH234, Histopath))


o.un = lm(Y1 ~ Treatment + MassDay0, data = datanew)

# summary(o.un)

anova(o.un)
res.lsm = lsmeans(o.un, pairwise ~ Treatment) # pairwise testing averaged over age & days
typeof(res.lsm$contrasts)
contrast = res.lsm$contrasts 
p1 = forest.plot(res.lsm$contrasts, -1,5,title = "Aged")
forest.plot(res.lsm$contrasts, -10,30,title = "Aged")
# forest.plot(res.lsm$contrasts, -4e3,1e4,title = "Aged")

##--------------------Fig.8(K)---------------------------##
res_means = res.lsm$lsmeans
plot_res = as.data.frame(res_means)
ggplot(plot_res, aes(color =Treatment), group = interaction(Treatment)) +
  geom_point( aes(Treatment,Y1), data = datanew,  alpha = 0.3) +
  geom_point(aes(Treatment,lsmean), size = 2) +
  geom_errorbar(aes(x = Treatment, ymin = lower.CL, ymax = upper.CL), width =0.4) +
  theme_bw(base_size=15) +  ylab('Y = 1/3(logViral + PenH234 + Histopath) - logAbs')+
  scale_color_manual(name = "", labels =trts, values = trt.colorPal) +
  scale_fill_manual(name = "", labels =trts, values = trt.colorPal) +
  # facet_wrap(~Age) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) 

