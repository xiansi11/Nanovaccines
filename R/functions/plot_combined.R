## Plot combined
# library(devtools)
# install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)
plot_combined = function(filename, 
                         p = 0.01, method = "none",age = "Aged", test.items = NULL,
                         no.double = FALSE, perc.cap = NULL, old.plotnames = F,
                         trt.pre = NULL, exclude.resp = c("GMCSF","IL27"),resp.align = NULL, trt.align = NULL,
                         cytokine.order = FALSE, costim.order = FALSE, relative_path = "data/"){
  cell_fun = function(j, i, x, y, w, h, fill) {
    if(!is.na(plot.Padj[i, j])){
      if(plot.Padj[i, j] < 0.001) {
        grid.text("****", x, y)
      } else if(plot.Padj[i, j] < 0.01) {
        grid.text("***", x, y)
      } else if(plot.Padj[i, j] < 0.05) {
        grid.text("**", x, y)
      }
    }
  }
  # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",  "fdr", "none")
  # filename = out.filename; no.double = T; exclude.resp = NULL;method = "fdr"; p = 0.01
  load(paste0("test.res/", filename,age,".res.all.rda"))
  if(no.double){res.all = res.all %>% dplyr::filter(!Treatment %in%
                            c("PN+CpG", "Mi+CpG", "PN+CDN", "Mi+CDN"))
  } else {
    res.all = res.all %>% dplyr::filter(Treatment %in%
                          c("PN+Mi+CpG","PN+CpG", "Mi+CpG", "CpG", "PN", "Micelles", "PN+Mi+CDN", "PN+CDN", "Mi+CDN", "CDN"))
  }
  
  if(!is.null(exclude.resp)){res.all = res.all %>% dplyr::filter(!Response %in% exclude.resp)}
  # head(res.all)
  res.all$Treatment = factor(res.all$Treatment, levels = unique(res.all$Treatment))
  if(old.plotnames){
    res.all$Response = factor(res.all$Response, levels = unique(res.all$Response))
  }else{
    # change response to the new plot response
    map = read.csv(paste0(relative_path,"Response_map.csv"), fileEncoding = 'UTF-8-BOM'); 
    response.new = unlist(res.all %>% dplyr::left_join(map, by ='Response.ID') %>% dplyr::select(ResponsePlot2))
    res.all$Response = factor(response.new, levels = unique(response.new))
  }
 
  
  # Prepare for plotting, calculate q value
  # Select the hypothesis testing we are interested in
  if(is.null(test.items)){test.items =paste0("Trt-Unsti|", age)}
  temp = res.all%>%dplyr::filter(Items == all_of(test.items))
  pvalues = temp$`Pr(>|t|)`; # hist(pvalues)
  p.adj = p.adjust(pvalues, method); # hist(p.adj)
  temp$p.adj = p.adj
  temp$log2FC = log(x = exp(temp$Estimate), base = 2)
  
  ### Heatmap
  mat.Estimate = temp %>% dplyr::select(c("log2FC", "Response", "Treatment")) %>%
                          tidyr::pivot_wider(values_from = "log2FC", 
                                             names_from = "Treatment")
  plot.Estimate = as.matrix(mat.Estimate[,-1])
  rownames(plot.Estimate) = unlist(mat.Estimate[,1])
  plot.Estimate
  
  
  
  mat.Padj = temp %>% dplyr::select(c("p.adj", "Response", "Treatment")) %>%
    tidyr::pivot_wider(values_from = "p.adj", 
                       names_from = "Treatment")
  plot.Padj = as.matrix(mat.Padj[,-1])
  rownames(plot.Padj) = unlist(mat.Padj[,1])
  plot.Padj
  if(max(temp$log2FC)>9){
    upper = 10
  } else{
    upper = 2
  }
  
  
  
  if(!is.null(resp.align)){
    plot.Estimate = plot.Estimate [match(resp.align, rownames(plot.Estimate)),]
    rownames(plot.Estimate) = resp.align
    plot.Padj = plot.Padj [match(resp.align, rownames(plot.Padj)),]
    rownames(plot.Padj) = resp.align
  }
  if(!is.null(trt.align)){
    plot.Estimate  =   plot.Estimate[,match(trt.align, colnames(plot.Estimate))]
    colnames(plot.Estimate) = trt.align
    plot.Padj = plot.Padj [,match(trt.align, colnames(plot.Padj))]
    colnames(plot.Padj) = trt.align
  }
  
  col_fun = colorRamp2(c(-1, 0, upper), c("blue", "white", "red"))
  ha = HeatmapAnnotation(simple_anno_size = unit(1, "cm"),
                         foo = 1:10, col = list(foo = col_fun))
  
  

  
  if(cytokine.order){
    ### Plot blocked cytokines plots
    inhibitory = c("IL10", "CCL22")
    Tcell.activ = c("IL12p70", "IFNalpha","IFNbeta", "IFNg",
                    "IL6", "IL27")
    inflam = c("IL1a", "IL1b", "IL6",  "TNFa", "IL12p40")
    chemokines = c("CXCL1", "CXCL2","CCL2", "CCL3", "CXCL9","CXCL10", "CCL5")
    list.all = c(Tcell.activ, inhibitory, inflam, chemokines, "GCSF")
    ind.all = match(list.all, rownames(plot.Estimate))
    ind.all = ind.all[which(!is.na(ind.all))]
    # Exlude "Other Category" by uncommenting the following 
    # ind.all = c(ind.all, setdiff(seq(nrow(plot.Estimate)), ind.all))
    mat = plot.Estimate [ind.all,]
    split = sapply(rownames(mat), function(x){
      if(x %in% Tcell.activ){
        return("Tactiv")
      } else if(x %in% chemokines){
        return("Chem")
      } else if(x %in% inflam){
        return("Inflam")
      } else if(x %in% inhibitory){
        return("Regu")
      } else if(x == "GCSF"){
        return("Growth")
      }else{
        return("Other")
      }
    })
    if("IL6" %in% names(split)){
      # Change 1 of IL6 to inflam
      split[which("IL6"==names(split))[1]] = "Inflam"
    }
    # Set the order of clusters
    split = factor(split, levels = c("Tactiv", "Inflam", "Regu", "Chem", "Growth", "Other"))
    # Updated q value
    plot.Padj = plot.Padj[ind.all, ]
    mat; split
    fig = Heatmap(mat, 
                  col = col_fun,
                  cell_fun = cell_fun,
                  show_row_dend = F,
                  show_column_dend = F,
                  # Reorder
                  row_split = split,
                  border = TRUE,
                  cluster_row_slices = F,
                  cluster_rows = F,
                  heatmap_legend_param = list(col_fun = col_fun, 
                                              title = test.items, 
                                              direction = "vertical")
    )
    if((!is.null(trt.align))|(!is.null(resp.align))){
      fig = Heatmap(mat, 
                    col = col_fun,
                    cell_fun = cell_fun,
                    show_row_dend = F,
                    show_column_dend = F,
                    # Reorder
                    row_split = split,
                    cluster_row_slices = FALSE,
                    border = TRUE,
                    cluster_rows = F,
                    cluster_columns = F,
                    heatmap_legend_param = list(col_fun = col_fun, 
                                                title = test.items, 
                                                direction = "vertical")
      )
    } 
    print(fig)
  } else if(costim.order){
    ### Plot blocked costims
    if(old.plotnames){
      CD40_subsets = c('CD40+CD11c+','CD40+DC','CD40+DC1','CD40+CD4+DC2','CD40+CD4-DC2','CD40+DC2','CD40+DNpDC','CD40+pDC','CD40+CD4+pDC','CD40+CD8+pDC')
      CD80_subsets = c('CD80+CD11c+','CD80+DC','CD80+DC1','CD80+CD4+DC2','CD80+CD4-DC2','CD80+DC2','CD80+DNpDC','CD80+pDC','CD80+CD4+pDC','CD80+CD8+pDC')
      CD86_subsets = c('CD86+CD11c+','CD86+DC','CD86+DC1','CD86+CD4+DC2','CD86+CD4-DC2','CD86+DC2','CD86+DNpDC','CD86+pDC','CD86+CD4+pDC','CD86+CD8+pDC')
    } else{
      CD40_subsets = c('CD40+CD11c+','CD40+CD4-cDC2','CD40+CD4+cDC2','CD40+CD4+CD11c(low)','CD40+CD8+CD11c(low)','CD40+CD11c(high)','CD40+cDC1','CD40+cDC2','CD40+DNpDC','CD40+CD11c(low)')
      CD80_subsets = c('CD80+CD11c+','CD80+CD4-cDC2','CD80+CD4+cDC2','CD80+CD4+CD11c(low)','CD80+CD8+CD11c(low)','CD80+CD11c(high)','CD80+cDC1','CD80+cDC2','CD80+DNpDC','CD80+CD11c(low)')
      CD86_subsets = c('CD86+CD11c+','CD86+CD4-cDC2','CD86+CD4+cDC2','CD86+CD4+CD11c(low)','CD86+CD8+CD11c(low)','CD86+CD11c(high)','CD86+cDC1','CD86+cDC2','CD86+DNpDC','CD86+CD11c(low)')
    }
    list.all = c(CD40_subsets, CD80_subsets, CD86_subsets)
    ind.all = match(list.all, rownames(plot.Estimate))
    ind.all = ind.all[which(!is.na(ind.all))]
    # ind.all = c(ind.all, setdiff(seq(nrow(plot.Estimate)), ind.all))
    mat = plot.Estimate [ind.all,]
    split = sapply(rownames(mat), function(x){
      if(x %in% CD40_subsets){
        return("CD40")
      } else if(x %in% CD80_subsets){
        return("CD80")
      } else if(x %in% CD86_subsets){
        return("CD86")
      } else{
        return("Other")
      }
    })
    # Set the order of clusters
    split = factor(split, levels = c("CD40", "CD80", "CD86","Other"))
    # Updated q value
    plot.Padj = plot.Padj[ind.all, ]
    mat; split
    print(Heatmap(mat, 
                  col = col_fun,
                  cell_fun = cell_fun,
                  show_row_dend = F,
                  show_column_dend = F,
                  # Reorder
                  row_split = split,
                  border = TRUE,
                  cluster_rows = F,
                  cluster_row_slices = FALSE,
                  heatmap_legend_param = list(col_fun = col_fun, 
                                              title = test.items, 
                                              direction = "vertical")
    ))
    
  } else {
    fig = Heatmap(plot.Estimate, 
                  col = col_fun,
                  cell_fun = cell_fun,
                  show_row_dend = F,
                  show_column_dend = F,
                  heatmap_legend_param = list(col_fun = col_fun, 
                                              title =  test.items, 
                                              direction = "vertical")
    )
    print(fig)
  }
  
  # fig1b= ggplot(temp) +
  #   theme_bw()+geom_point(aes(Response, Treatment, color = Test1, size = 1-p.adj))+
  #   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=0.95))+
  #   ggtitle(label = "Significance of Treatment") +
  #   scale_colour_manual(values = c("Down" = "blue",
  #                                  "Up"= "red"))+
  #   scale_size(range=c(0,3),breaks=c(1-p),labels=c(1-p),guide="legend")+ coord_flip()+
  #   theme(legend.position = "bottom") +
  #   ggtitle(label = paste0("Response of Aged mice(",filename,")"))
  # print(fig1b)
  # 
  # # Plot1b: regulation of treatments
  # fig1b= ggplot(temp) +
  #   theme_bw()+
  #   geom_point(aes(Response, Treatment, color = Test, 
  #                  shape = 'square', size = log2FC))+
  #   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=0.95))+
  #   scale_colour_manual(values = c("Down" = "blue",
  #                                  "Up"= "red",
  #                                  "Insignicant"= "grey"))+
  #   scale_size(range=c(1,5), breaks = c(-1,0,1,2,4), guide="legend")+ 
  #   guides(color = guide_legend(
  #     override.aes=list(shape = "square")),
  #     size = guide_legend(
  #       override.aes=list(shape = "square")))+ coord_flip()+
  #   scale_shape_identity()  +
  #   labs(x="Cell responses", y="Adjuvants", size="Log2FC", 
  #        col="Regulation test") +
  #   theme(legend.position = "right") +
  #   ggtitle(label = paste0("Response of Aged mice(",filename,")"))
  # print(fig1b)
  # 
  # fig1b= ggplot(temp) +
  #   theme_bw()+geom_point(aes(Response, Treatment, color = Test, 
  #                             shape = 'circle', size = -log10p.adj))+ 
  #   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=0.95))+
  #   ggtitle(label = "Significance of Treatment") +
  #   scale_colour_manual(values = c("Down" = "blue",
  #                                  "Up"= "red",
  #                                  "Insignicant"= "grey"))+
  #   scale_size(range=c(1,5),breaks=c(2,10),guide="legend")+ 
  #   scale_shape_identity() +
  #   guides(color = guide_legend(
  #     override.aes=list(shape = "circle")),
  #     size = guide_legend(
  #       override.aes=list(shape = "circle")))+ coord_flip()+
  #   theme(legend.position = "bottom") +
  #   ggtitle(label = paste0("Response of Aged mice(",filename,")"))
  # print(fig1b)
  # 
  # 
}

