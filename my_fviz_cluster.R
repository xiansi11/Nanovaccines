my_fviz_cluster =
  function(km.res, df, cbPalette, tt, aa){
    fig.clust = factoextra::fviz_cluster(km.res, data = df,
                                         palette = cbPalette, 
                                         ggtheme = theme_minimal(), labelsize = 10, geom=c("point","text"),
                                         main = paste(tt, paste0(aa,collapse = "")), xlab = "", ylab ="",show.legend = F,repel = 1) + 
      theme_bw() + 
      xlab(expression(Delta ~ PC1)) + ylab(expression(Delta ~ PC2)) +
      guides(fill=guide_legend(title="Cluster"), 
                  color = guide_legend(title="Cluster"),
                  shape = guide_legend(title="Cluster"))+
      theme(legend.position = "right",
            axis.text=element_text(size=14),
            axis.title=element_text(size=16,face="bold"),
            legend.title = element_text(size=14), #change legend title font size
            legend.text = element_text(size=14))
    fig.clust$layers[[4]]$mapping$colour = NULL
    fig.clust$layers[[4]]$aes_params$size = 5
    fig.clust
  }