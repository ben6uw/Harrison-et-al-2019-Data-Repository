pca_set2 <- function(input, sample_df, title){
  pca <- prcomp(t(na.omit(input)),scale. = TRUE)
  var_exp = pca$sdev^2
  pve= round((var_exp/sum(var_exp))*100,1)
  print(pve[1:3])
  
  pca <- as.data.frame(pca$x[,1:5])
  pca <- merge(pca, sample_df, by.x ="row.names", by.y = "Sample.ID")
  pca = pca[match(sample_df$Sample.ID, pca$Row.names),]
  
  pca.3d = plot_ly(data = pca, x= ~PC1, y= ~PC2, z= ~PC3, 
                   color = ~line, 
                   symbol = ~treatment,
                   symbols = c("circle","circle-open"),
                   size = ~runOrder,
                   sizes= c(10,70),
                   hoverinfo="text", 
                   text = ~paste0("<br>Sample:",Row.names,
                                  "<br>trait:",trait,
                                  "<br>treatment:",treatment,
                                  "<br>Run.order:",runOrder,
                                  "<br>line.weight:",line.weight,
                                  "<br>batch:",block),
                   marker = list(sizemode = "diameter", opacity=0.6))%>%
    add_markers() %>%
    layout(title = "3D PCA - line coded by color, size by runOrder",
           legend = list(size =10),
           scene = list(xaxis = list(title = "PC1"),
                        yaxis = list(title = "PC2"),
                        zaxis = list(title = "PC3")))
  
  setwd("plots")
  htmlwidgets::saveWidget(pca.3d,paste0(title,"_PCA_3D.html"))
  setwd("..")
  
  return(list("pca_df"=pca))
}
