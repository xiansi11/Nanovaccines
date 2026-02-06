##------------------------------------------------------###
##---------------------Fig.2-----------------------------####
##------------------------------------------------------###


library(ggrepel)
library(ggplot2)
library(devtools)
# devtools::install_github("vqv/ggbiplot")
library(ggbiplot)
library(tidyr)
library(dplyr)
# install.packages("ggh4x")
library(ggh4x)
library(lme4)
library(lmerTest)
library(psych)
source("plot_pca.R")
source("my_fviz_cluster.R")

##---------------------Fig.2(C to J)-----------------------------####
map = read.csv("data/Treatment_map.csv")
ex.resp.ls = NULL
cbPalette <- c("#E69F00", "#0072B2","#CC79A7",
               "#D55E00","#999999",  "#56B4E9", "#009E73", "#F0E442")

aa = "Aged"
filename = "SplenicDC_combined_costim_March_April_June_August"
for (tt in c("cDC1", "cDC2", "CD4+cDC2","CD4-cDC2")){
  res_tb = plot_pca(filename, combine = T, age =aa,
                    resp.ex = ex.resp.ls,
                    costim.group = tt, no.double = T)
  cat(paste0(aa,tt))
  # print(res_tb)
  dat = read.csv(paste0("test.res/PCA_combine_",tt,"_", filename, aa,".csv"))
  df = dat %>% dplyr::select(disPC1, disPC2)
  rownames(df) = dat$Treatment
  res = factoextra::fviz_nbclust(df, kmeans, method = c("silhouette", "wss", "gap_stat"))
  # print(res)
  km.res <- kmeans((df), centers = which.max(res$data$y), nstart = 25)
  fig.clust = my_fviz_cluster(km.res, df, cbPalette,tt,aa)
  
   print(fig.clust)
  cat("\n")
}



aa = "Young"
filename = "SplenicDC_combined_costim_March_April_June_August"
for (tt in c("cDC1", "cDC2", "CD4+cDC2","CD4-cDC2")){
  cat("\n")
  res_tb = plot_pca(filename, combine = T, age =aa, resp.ex = ex.resp.ls, costim.group = tt, no.double = T)
  cat(paste0(aa,tt))
  # print(res_tb)
  dat = read.csv(paste0("test.res/PCA_combine_",tt,"_", filename, aa,".csv"))
  df = dat %>% dplyr::select(disPC1, disPC2)
  rownames(df) = dat$Treatment
  res = factoextra::fviz_nbclust(df, kmeans, method = c("silhouette", "wss", "gap_stat"))
  # print(res)
  km.res <- kmeans((df), centers = which.max(res$data$y), nstart = 25)
  fig.clust = my_fviz_cluster(km.res, df, cbPalette,tt,aa)
  
   print(fig.clust)
  cat("\n")
}

#