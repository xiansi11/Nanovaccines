# install.packages("factoextra")
library(factoextra)
library(dplyr)
Tcell.activ = c("IL12p70", "IFNalpha","IFNbeta", "IFNg",
                "IL6", "IL27")
inflam = c("IL1a", "IL1b", "IL6",  "TNFa", "IL12p40")

# Fixed color pallettes
Treatment = c('Unstimulated','NP','Micelles','CpG','CDN',
              'R848','FLA-BS','MDP','LPS','MPLA',
              'NP+Mi+CpG','NP+Mi+CDN','NP+Mi+R848',
              'NP+CpG','Mi+CpG','NP+CDN','Mi+CDN')
# gg_color_hue <- function(n) {
#   hues = seq(15, 375, length = n + 1)
#   hcl(h = hues, l = 65, c = 100)[1:n]
# }
# trt.colorPal = c("black", gg_color_hue(16))
trt.colorPal = c(#'#990F0F','#B22C2C','#CC5151','#E57E7E','#FFB2B2',
  '#990F0F','#CC5151','#E57E7E','#FFB2B2',
  #'#99540F',
  '#B26F2C','#CC8E51','#E5B17E',#'#FFD8B2',
  #'#6B990F',
  '#85B22C','#A3CC51','#C3E57E',#'#E5FFB2',
  #'#0F6B99',
  '#2C85B2','#51A3CC','#7EC3E5',#'#B2E5FF',
  #'#260F99',
  '#422CB2','#6551CC','#8F7EE5','#BFB2FF'
)
# colorBlindness::displayAvailablePalette()
# paste0(colorBlindness::SteppedSequential5Steps,collapse = "\',\'"
# trt.shapePal = c(19, 15, 17, 18, 0:12)
trt.shapePal = c(1, rep(c(15, 17, 4),5), 18)
trt.Pal = as.data.frame(cbind(Treatment, trt.colorPal, trt.shapePal))


## 1a. Separate PCA on Tcell and Inflam 
dat1 = read.csv("./test.res/PCA_separate_detail_August_Inflammation_BMDC_combined_cytokines_June_August.csv")
dat2 = read.csv("./test.res/PCA_separate_detail_August_T cell activation_BMDC_combined_cytokines_June_August.csv")
head(dat1)
head(dat2)
dim(dat1)
dim(dat2)
df = data.frame(dat1%>% dplyr::select(Mouse.ID, Treatment,Treatment.ID,Age),dat1$PC1,dat2$PC1)
colnames(df) = c("Mouse.ID", "Treatment","Treatment.ID", "Age", "Inflam.PC1", "Tcell.PC1")
headers = c("Mouse.ID","Treatment.ID","Treatment","Age")
temp = df %>% dplyr::filter(Age == "Aged") 
temp.Pal = trt.Pal %>% dplyr::filter(Treatment %in% unique(temp$Treatment))

ggplot(aes(x = Inflam.PC1, y = Tcell.PC1), data = df)+
  geom_point(aes(shape=factor(temp$Treatment, levels = temp.Pal$Treatment),
            colour=factor(temp$Treatment, levels = temp.Pal$Treatment))
  ) +
  scale_shape_manual(values=as.numeric(temp.Pal$trt.shapePal))+
  scale_color_manual(values=(temp.Pal$trt.colorPal))+ theme_bw()+
  labs(color = "Adjuvants", shape = "Adjuvants")+ 
  coord_equal(ratio = 1.2) +
  ggtitle(paste0(""))

df%>% group_by(Treatment) %>% dplyr::summarise(n = n())
temp1 = temp %>% dplyr::select(Inflam.PC1, Tcell.PC1)


# 
# ## 1b. simple average on Tcell and Inflam
# dat = read.csv(paste0("./test.res/PCA_combine_detail_BMDC_combined_cytokines_June_August.csv"))
# colnames(dat)
# Tcell.activ = c("IL12p70", "IFNalpha","IFNbeta", "IFNg",
#                 "IL6", "IL27")
# inflam = c("IL1a", "IL1b", "IL6",  "TNFa", "IL12p40")
# 
# ind1 = match(Tcell.activ, colnames(dat)); ind1 = ind1[which(!is.na(ind1))]
# ind2 = match(inflam, colnames(dat)); ind2 = ind2[which(!is.na(ind2))]
# colnames(dat)[ind1]
# colnames(dat)[ind2]
# df = dat %>% mutate(
#   Tcell.ave = rowMeans(dat[,ind1]),
#   Inflam.ave = rowMeans(dat[,ind2])
# )
# temp = df %>% dplyr::filter(Age == "Aged") 
# temp.Pal = trt.Pal %>% dplyr::filter(Treatment %in% unique(temp$Treatment))
# ggplot(aes(x = Inflam.ave, y = Tcell.ave), data = df)+
#   geom_point(aes(shape=factor(temp$Treatment, levels = temp.Pal$Treatment),
#                  colour=factor(temp$Treatment, levels = temp.Pal$Treatment))
#   ) +
#   scale_shape_manual(values=as.numeric(temp.Pal$trt.shapePal))+
#   scale_color_manual(values=(temp.Pal$trt.colorPal))+ theme_bw()+
#   labs(color = "Adjuvants", shape = "Adjuvants")+ 
#   coord_equal(ratio = 1.2) +
#   ggtitle(paste0(""))
# temp1 = temp %>% dplyr::select(Inflam.ave, Tcell.ave)

## 2. PCA on Tcell.PC1 and Inflam.PC1
# Select responses that were repeated in two experiments
res.pca <- prcomp(temp1,
                  center = TRUE,scale. = TRUE)
cc = summary(res.pca)
cc$importance
arrow.length = 1
label.adjust = 1.4
my_ggbiplot(res.pca, ellipse=T, 
            # labels=temp$Treatment, 
            groups = factor(temp$Treatment, levels = temp.Pal$Treatment),
            # obs.scale = 2,
            # var.scale = 2,
            labels.size = 1, 
            varname.size = 3,
            # varname.adjust = 1.2,
            varname.abbrev = F,
            alpha=0,
            arrow.length = arrow.length,
            label.adjust = label.adjust
) + theme_bw(base_size=10) + theme(legend.position = "right") +
  geom_point(aes(shape=factor(temp$Treatment, levels = temp.Pal$Treatment),
                 colour=factor(temp$Treatment, levels = temp.Pal$Treatment))
  ) +
  scale_shape_manual(values=as.numeric(temp.Pal$trt.shapePal))+
  scale_color_manual(values=(temp.Pal$trt.colorPal))+
  theme(legend.position = "right") +
  labs(color = "Adjuvants", shape = "Adjuvants", 
       x = paste0("PC1 (", round(cc$importance[2,1],3)* 100, "%)"),
       y = paste0("PC2 (", round(cc$importance[2,2],3)* 100, "%)"))+ 
  coord_equal(ratio = 0.7) +
  ggtitle(paste0(""))



