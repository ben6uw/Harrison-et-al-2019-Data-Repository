pca_set <- function(input, sample_df, title){
  pca <- prcomp(t(na.omit(input)),scale. = TRUE)
  var_exp = pca$sdev^2
  pve= round((var_exp/sum(var_exp))*100,1)
  #print(pve[1:3])
  
  pca <- as.data.frame(pca$x[,1:3])
  pca <- merge(pca, sample_df, by.x ="row.names", by.y = "Sample.ID")
  pca = pca[match(sample_df$Sample.ID, pca$Row.names),]
  
  pca.3d = plot_ly(data = pca, x= ~PC1, y= ~PC2, z= ~PC3, 
                   color = ~trait, colors = c("#8B8378", "#EE9A00"),  
                   symbol = ~treatment,
                   symbols = c("circle","x"),
                   hoverinfo="text", 
                   text = ~paste0("<br>Sample:",Row.names,
                                  "<br>trait:",trait,
                                  "<br>lifespan:",lifespan,
                                  "<br>treatment:",treatment,
                                  "<br>Run.order:",runOrder,
                                  "<br>line.weight:",line.weight,
                                  "<br>sampleWeight:",sampleWeight,
                                  "<br>batch:",block),
                   marker = list(size=6))%>%
    add_markers() %>%
    layout(title = "3D PCA - genotype coded by color,treatment type coded by symbol",
           legend = list(size =10),
           scene = list(xaxis = list(title = "PC1"),
                        yaxis = list(title = "PC2"),
                        zaxis = list(title = "PC3")))
  
  setwd("plots")
  htmlwidgets::saveWidget(pca.3d,paste0(title,"_PCA_3D.html"))
  setwd("..")
  
  pca.line = plot_ly(data = pca, x= ~PC1, y= ~PC2, z= ~PC3, 
                   color = ~line, colors = "Dark2",  
                   symbol = ~treatment,
                   symbols = c("circle","x"),
                   hoverinfo="text", 
                   text = ~paste0("<br>Sample:",Row.names,
                                  "<br>trait:",trait,
                                  "<br>lifespan:",lifespan,
                                  "<br>treatment:",treatment,
                                  "<br>Run.order:",runOrder,
                                  "<br>line.weight:",line.weight,
                                  "<br>batch:",block),
                   marker = list(size=6))%>%
    add_markers() %>%
    layout(title = "3D PCA - line type coded by color, treatment coded by symbol",
           legend = list(size =10),
           scene = list(xaxis = list(title = "PC1"),
                        yaxis = list(title = "PC2"),
                        zaxis = list(title = "PC3")))
  
  setwd("plots")
  htmlwidgets::saveWidget(pca.line,paste0(title,"_PCA_line.html"))
  setwd("..")
  
  pca.batch = plot_ly(data = pca, x= ~PC1, y= ~PC2, z= ~PC3, 
                      color = ~block,
                      hoverinfo="text", 
                      text = ~paste0("<br>Sample:",Row.names,
                                     "<br>trait:",trait,
                                     "<br>lifespan:",lifespan,
                                     "<br>treatment:",treatment,
                                     "<br>Run.order:",runOrder,
                                     "<br>line.weight:",line.weight,
                                     "<br>batch:",block),
                      marker = list(size=6))%>%
    add_markers() %>%
    layout(title = "3D PCA - batch coded by color ",
           legend = list(size =10),
           scene = list(xaxis = list(title = "PC1"),
                        yaxis = list(title = "PC2"),
                        zaxis = list(title = "PC3")))
  
  setwd("plots")
  htmlwidgets::saveWidget(pca.batch,paste0(title,"_PCA_batch.html"))
  setwd("..")
  
  
  pca.order = plot_ly(data = pca, x= ~PC1, y= ~PC2, z= ~PC3, 
                      color = ~runOrder,
                      symbol = ~comb,
                      symbols = c("circle","x","diamond","square"),
                      hoverinfo="text", 
                      text = ~paste0("<br>Sample:",Row.names,
                                     "<br>trait:",trait,
                                     "<br>lifespan:",lifespan,
                                     "<br>treatment:",treatment,
                                     "<br>Run.order:",runOrder,
                                     "<br>line.weight:",line.weight,
                                     "<br>batch:",block),
                      marker = list(size=6))%>%
    add_markers() %>%
    layout(title = "3D PCA - run order coded by color",
           legend = list(size =10),
           scene = list(xaxis = list(title = "PC1"),
                        yaxis = list(title = "PC2"),
                        zaxis = list(title = "PC3")))
  
  setwd("plots")
  htmlwidgets::saveWidget(pca.order, paste0(title,"_PCA_Runorder.html"))
  setwd("..")
  
  pca.sampleWeight = plot_ly(data = pca, x= ~PC1, y= ~PC2, z= ~PC3, 
                      color = ~sampleWeight,
                      symbol = ~comb,
                      symbols = c("circle","x","diamond","square"),
                      hoverinfo="text", 
                      text = ~paste0("<br>Sample:",Row.names,
                                     "<br>trait:",trait,
                                     "<br>lifespan:",lifespan,
                                     "<br>treatment:",treatment,
                                     "<br>Run.order:",runOrder,
                                     "<br>line.weight:",line.weight,
                                     "<br>batch:",block),
                      marker = list(size=6))%>%
    add_markers() %>%
    layout(title = "3D PCA - sample weight coded by color",
           legend = list(size =10),
           scene = list(xaxis = list(title = "PC1"),
                        yaxis = list(title = "PC2"),
                        zaxis = list(title = "PC3")))
  
  setwd("plots")
  htmlwidgets::saveWidget(pca.sampleWeight, paste0(title,"_PCA_sampleweight.html"))
  setwd("..")
  
  pca.linewg = plot_ly(data = pca, x= ~PC1, y= ~PC2, z= ~PC3, 
                       color = ~line.weight,
                       symbol = ~trait,
                       symbols = c("circle","x","diamond","square"),
                       hoverinfo="text", 
                       text = ~paste0("<br>Sample:",Row.names,
                                      "<br>trait:",trait,
                                      "<br>lifespan:",lifespan,
                                      "<br>treatment:",treatment,
                                      "<br>Run.order:",runOrder,
                                      "<br>line.weight:",line.weight,
                                      "<br>batch:",block),
                       marker = list(size=6))%>%
    add_markers() %>%
    layout(title = "3D PCA - line weight coded by color, treatment & genotype coded by symbol",
           legend = list(size =10),
           scene = list(xaxis = list(title = "PC1"),
                        yaxis = list(title = "PC2"),
                        zaxis = list(title = "PC3")))
  
  setwd("plots")
  htmlwidgets::saveWidget(pca.linewg, paste0(title,"_PCA_linewg.html"))
  setwd("..")
  
  pca.sampmiss = plot_ly(data = pca, x= ~PC1, y= ~PC2, z= ~PC3, 
                       color = ~percent_missings,
                       #symbol = ~comb,
                       #symbols = c("circle","x","diamond","square"),
                       hoverinfo="text", 
                       text = ~paste0("<br>Sample:",Row.names,
                                      "<br>trait:",trait,
                                      "<br>lifespan:",lifespan,
                                      "<br>treatment:",treatment,
                                      "<br>Run.order:",runOrder,
                                      "<br>line.weight:",line.weight,
                                      "<br>batch:",block),
                       marker = list(size=6))%>%
    add_markers() %>%
    layout(title = "3D PCA - percent missings coded by color, treatment & genotype coded by symbol",
           legend = list(size =10),
           scene = list(xaxis = list(title = "PC1"),
                        yaxis = list(title = "PC2"),
                        zaxis = list(title = "PC3")))
  
  setwd("plots")
  htmlwidgets::saveWidget(pca.sampmiss, paste0(title,"_PCA_sampmiss.html"))
  setwd("..")
  

  outlist = list("pca_df"=pca,
                 "pve"=pve[1:3],
                 "pca.3d"=pca.3d,
                 "pca.line"=pca.line,
                 "pca.batch"=pca.batch,
                 "pca.order"=pca.order,
                 "pca.sampleWeight"=pca.sampleWeight,
                 "pca.linewg"=pca.linewg,
                 "pca.sampmiss"=pca.sampmiss)
  return(outlist)
}

