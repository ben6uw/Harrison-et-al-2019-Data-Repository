
  
  
```{r setup, include=FALSE}
  library(knitr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(BiocStyle)
  library(limma)
  library(edgeR)
  library(threejs)
  library(htmlwidgets)
  library(xtable)
  library(ReportingTools)
  library(plotly)
  library(openxlsx)
  library(DT)
  library(affycoretools)
  library(Glimma)
  library(ComplexHeatmap)

```


```{r parameter, include=FALSE}
opts_chunk$set(echo=FALSE, fig_path="plots/tmp")
# set parameters
fdrcutoff=0.1
missingCutoff =5
datamode = "positive"
removecol="negLabel"
datapath="../../../Metabolomics2016/Newdatafiles/pid1652_Positive_sample/XCMS-Posdata-93samples-Promislow.xlsx"
```

# load unique functions:
source("./pca_set.R",local = T) 
source("./pca_set2.R",local = T)
source("./pca_set3.R",local = T)
source("./heatmaps_set.R",local = T)

```{r main, include=FALSE}
setwd(paste0("/data5/analyses/Daniel_Promislow/Ben_Harrison_pid2005/",datamode,"_mode/bioconductor"))
# read in sample information.
samps= read.csv("../../../Metabolomics2016/OlddataFiles/sample.info.noqc.csv",stringsAsFactors = F)
samps= select(samps, -removecol)
samps$treatment <- as.factor(ifelse(samps$treatment =="B-", "Ctr", "H2O2"))
samps$trait <- factor(ifelse(samps$lifespan >= 80, "resistant", "sensitive"),
                      levels = c("sensitive","resistant"))

colnames(samps)[1] <- "Sample.ID"
# re-order by the run order
samps = samps[order(samps$runOrder),]
# Some summary tables
table(samps$treatment)
table(samps$trait)
table(samps$treatment,samps$trait)
# 16 lines
samps$line=gsub("-","_", samps$line)
length(unique(samps$line))
samps$line = as.factor(samps$line)
# 3 BATCH!!!
table(samps$line,samps$block)
samps$block = as.factor(samps$block)
samps$comb = factor(paste(samps$trait,samps$treatment, sep="_"),
                       levels=c("sensitive_Ctr", "sensitive_H2O2", "resistant_Ctr", "resistant_H2O2"))
#-----------------------------------------------------
raw.data = read.xlsx(datapath)
colnames(raw.data) = gsub("name","metabolites", colnames(raw.data))
any(duplicated(raw.data$metabolites)) #FALSE
sum(! is.na(raw.data$CV))
row.names(raw.data) = raw.data$metabolites
colnames(raw.data) <- gsub("\\.mzdata", "", colnames(raw.data))
colnames(raw.data) <- gsub("_", "", colnames(raw.data))
# one sample named "neg46.2"; not sure what dose'.2' means;
(colnames(raw.data) <- gsub("\\.2","", colnames(raw.data)))
stopifnot(all(samps$Sample.ID %in% colnames(raw.data)))
#features info
feature_info = select(raw.data, featureidx:postHoc)
# do log2 transformation, replace zero with NA.
cleanup <- function(input, selected_column){
  dat= dplyr::select(input, starts_with(selected_column))

  print(colnames(dat))
  dat[dat==0] <-NA
  dat = as.matrix(log2(dat))

  print(paste0("missing%=",round(mean(is.na(dat))*100)))
  print(apply(na.omit(dat),2, range))
  
  outlist = list("dat" = dat)
  return(outlist)
}
#----------------------------------
# QC sample
sample.qc =  cleanup(raw.data, "QC")$dat
#QCname different from negative mode, only QC1 not QCpos1 
o=order(as.integer(gsub("QC","",colnames(sample.qc)))) 
# re-order as the runorder
sample.qc = sample.qc[,o]
colnames(sample.qc)
# create a CV% for QC sample.
CV= data.frame(metabolites = row.names(sample.qc),
               CV=(apply(sample.qc,1,function(x) 
                 sd(x, na.rm = T))/rowMeans(sample.qc,na.rm = T))*100,
                                stringsAsFactors = F)
#----------------------------------
# column names in raw.data is UWMRC
substr(datamode,1,3)
raw.matrix =  cleanup(raw.data,substr(datamode,1,3))$dat
raw.matrix = raw.matrix[, match(samps$Sample.ID,colnames(raw.matrix))]
stopifnot(all.equal(colnames(raw.matrix), samps$Sample.ID)) #TRUE
stopifnot(dim(raw.matrix)[2]== nrow(samps))
samps$percent_missings = round(colMeans(is.na(raw.matrix))*100,2)
stopifnot(all.equal(row.names(sample.qc), row.names(raw.matrix)))
#----------------------------------
# EDA
if(!dir.exists("plots")) dir.create("plots")
mycol <- c("#00B2EE", "#FFC125", "#A8A8A8")

pdf(paste0("./plots/Density_plot_",datamode,".pdf"))
plot(density(raw.matrix[,1], na.rm=TRUE), col=mycol[1], lty=2,   main = "Desnsity plot") 
apply(raw.matrix[,2:ncol(raw.matrix)], 2, function(x) lines(density(x, na.rm=TRUE), col=mycol[1], lty=2))
apply(sample.qc, 2, function(x) lines(density(x, na.rm=TRUE), col=mycol[2]))
legend("topright", col=mycol, legend = c("Samples", "QC"), lty =1)
dev.off()
#-----------------------------------
# BOX plot
pdf(paste0("plots/boxplots_",datamode,".pdf"), width = 12)
boxplot(raw.matrix,las=2,col=mycol[1],ylab="log2 abundance", main="Sample")
boxplot(sample.qc,las=2,col=mycol[2],ylab="log2 abundance", main="QC")
dev.off()
#---------------------------------------------
pdf(paste0("./plots/sample_hist_",datamode,".pdf"))
par(mfrow = c(2,2))
hist(colSums(raw.matrix,na.rm = T),breaks = 40, col="light blue",  main = "Sum of log2(abundance) by sample")
hist(colMeans(raw.matrix,na.rm = T),breaks = 40, col="yellow",main = "Mean of log2(abundance) by sample ")
hist(rowSums(raw.matrix,na.rm = T),breaks = 40, col="light blue", main = "Sum of log2(abundance) by feature ")
hist(rowMeans(raw.matrix,na.rm = T),breaks = 40,col="yellow", main = "Mean of log2(abundance) by feature ")
hist(colMeans(is.na(raw.matrix)),breaks = 40, 
     main = "Percentage of missings by sample")

hist(rowMeans(is.na(raw.matrix)),breaks = 40,
     main = "Percentage of missings by feature")
dev.off()
#-----------------------------------------
# Missing data
# overall percent of missing data 
mean(is.na(raw.matrix))*100 #2%
# percent of features have < 5% missing data
mean(rowMeans(is.na(raw.matrix))*100 < 5)*100
# percent missing per sample
colMeans(is.na(raw.matrix))*100
#-----------------------------------
# missing by feature, across all samples  
missing_by_feature = data.frame(metabolites = row.names(raw.matrix),
                                mean_log2abundance_noNA = rowMeans(raw.matrix,na.rm = T),
                                percent_missing = rowMeans(is.na(raw.matrix))*100,
                                stringsAsFactors = F)

missing_by_feature = dplyr::left_join(missing_by_feature,CV,by="metabolites")
# plot the percent_missing over mean_mz
pdf(paste0("./plots/percent_missing_meanabundance_",datamode,".pdf"))
plot(missing_by_feature$mean_log2abundance_noNA, missing_by_feature$percent_missing, 
     pch=19, col="light blue", xlab="mean log2(abundance)", ylab="Percent of missingness", 
     main = datamode)
points(missing_by_feature$mean_log2abundance_noNA[missing_by_feature$percent_missing>=missingCutoff],
       missing_by_feature$percent_missing[missing_by_feature$percent_missing>=missingCutoff], pch=19, col="orange")
legend("topright",
       legend=c(paste0("<",missingCutoff,"% missingness"), 
                paste0(">=",missingCutoff,"% missingness")),pch=19, cex=0.8,
       col=c("light blue", "orange"))
dev.off()
#----------------------------------------------
table(missing_by_feature$percent_missing)
samps %>% group_by(comb) %>% 
  summarise(mean(percent_missings))

ms.plot = ggplot(samps, aes(x= reorder(Sample.ID, percent_missings), y= percent_missings,       
                 fill= comb)) +
                 geom_bar(stat="identity") + 
                 labs(title = paste0("Sample missing data "), 
                 x = "",  y= "Percent of missingness", fill="") +
                 theme(legend.position = "top",
                        axis.text.x=element_text(angle=90, hjust=1)) 

pdf("./plots/Percent_missing_sample.pdf")
ms.plot
dev.off()
#-----------------------------------------------
# 3D-pca
source("pca_set.R",local = T)
pca_set(input = raw.matrix,sample_df = samps,title = paste0("Raw_data_",datamode))
#----------------------------------------------
# median normalization 
apply(raw.matrix,2, function(x) median(x,na.rm = T))
(median.global <- median(raw.matrix, na.rm = TRUE))
dat.matrix <- apply(raw.matrix, 2, function(x) x - median(x, na.rm = TRUE) + median.global)

pdf(paste0("./plots/normalization_boxplots_",datamode,".pdf"))
boxplot(raw.matrix,las=2,col=mycol[1],ylab="log2(abundance)", main="Before normalization", cex.axis =0.6)
boxplot(dat.matrix,las=2,col=mycol[2],ylab="log2(abundance)", main="After normalization", cex.axis = 0.6)
dev.off()

# select features with missingness <5
select_feature = missing_by_feature[missing_by_feature$percent_missing <missingCutoff & missing_by_feature$CV<30,]$metabolites

dat.matrix.filter = dat.matrix[row.names(dat.matrix) %in% select_feature,]
dim(dat.matrix.filter) 
mean(is.na(dat.matrix.filter))*100
range(apply(dat.matrix.filter, 2, function(x) mean(is.na(x))*100))
limma::plotDensities(dat.matrix.filter,legend = F)
# use KNN imputation
impute.matrix = impute::impute.knn(dat.matrix.filter)$data
pca_set(input = impute.matrix,sample_df = samps,title = paste0("Imputed_",datamode))
#limma::plotDensities(impute.matrix,legend = F)
#---------------------------------------------------------
stopifnot(all.equal(colnames(impute.matrix),samps$Sample.ID))
#----------------------------------------------
# check drifting
zz.long = data.frame(raw.matrix) %>%
  mutate(metabolites = row.names(.))%>%
  gather(., Sample.ID, log2abundance,-metabolites) %>%
  left_join(., samps,by="Sample.ID")

pdf(paste0("plots/",datamode,"_metabolite_dotplot_before_normalization.pdf"),height = 180,width = 18)
ggplot(zz.long[zz.long$metabolites %in% select_feature[1:200],],aes(x=runOrder, y =log2abundance,color=comb))+
  geom_point()+
  facet_wrap(~metabolites, ncol=2)+
  theme(legend.position = "top", legend.direction = "horizontal")
dev.off()
#---------------------------------------------
# check what other covariates need to be adjusted besides batch.
noBatch = removeBatchEffect(impute.matrix,design = model.matrix(~comb,samps),
                        batch = samps$block)
pca_set(noBatch, samps, title= paste0("noBatch_",datamode))

#---------------------------------------------------------
# use limma pipeline
#----------------------------------------
levels(samps$comb)
design <- model.matrix(~0 + comb + block + sampleWeight, data = samps)
colnames(design) = gsub("comb","", colnames(design))
colnames(design)
contrast = makeContrasts(resistant_H2O2 - resistant_Ctr,
                         sensitive_H2O2 - sensitive_Ctr,
                         resistant_Ctr - sensitive_Ctr,
                         resistant_H2O2 - sensitive_H2O2,
                         (resistant_H2O2 - resistant_Ctr) - (sensitive_H2O2 - sensitive_Ctr),
                         levels = design)
colnames(contrast) = gsub(" - ","_vs_", colnames(contrast))

fit <- lmFit(impute.matrix, design = design)
fit <- contrasts.fit(fit, contrast)
fit <- eBayes(fit,robust = T,trend = T)
plotSA(fit)
out <- lapply(1:ncol(contrast), function(x) topTable(fit, coef = x, n= Inf, adjust = "BH",p.value=fdrcutoff, sort.by = "p")[,-6])
sapply(out,nrow) 
#-------------------
out2 <- lapply(1:ncol(contrast), function(x) topTable(fit, coef = x, n= Inf, adjust = "BH",p.value=1, sort.by = "none")[,-6])
pdf("plots/hist_pvalue.pdf")
par(mfrow = c(2,2))
for (i in 1:length(out2)){
   hist(out2[[i]]$P.Value, main=paste0("out2_",colnames(contrast)[i]))
}
dev.off()
#---------------------------------------------
sapply(out2,nrow)
names(out2) = colnames(contrast)
outputname = colnames(contrast)
outputname[ncol(contrast)] <- "Interaction"
#---------------------------------------------
# remove batch and sample weight for plotting purpose only.
noCV = removeBatchEffect(impute.matrix,design = model.matrix(~comb,samps),
                        batch = samps$block, covariates = samps$sampleWeight)
pca_set(noCV, samps, title= paste0("noCV_",datamode))
for (i in 1:length(unique(samps$treatment))){
  subset = samps[samps$treatment ==levels(samps$treatment)[i],]
  pca_set2(noCV[, match(subset$Sample.ID,colnames(noCV))], subset, 
           title= paste0("noCV_2",datamode,levels(samps$treatment)[i]))
}

noCV_pca= pca_set(noCV, samps, title=paste0("noCV_",datamode))$pca_df
pve=pca_set(noCV, samps, title=paste0("noCV_",datamode))$pve
#---------------------------------------------
if(!file.exists("glimma_plots")) dir.create("glimma_plots")
for (i in 1:length(out2)){
  out2[[i]]$metabolites= row.names(out2[[i]])
  out2[[i]]= select(out2[[i]], metabolites, everything())
 
  status = ifelse((out2[[i]]$adj.P.Val<fdrcutoff & out2[[i]]$logFC>0),1,0)
  status = ifelse((out2[[i]]$adj.P.Val<fdrcutoff & out2[[i]]$logFC<0),-1,status)
  nameregrx = gsub("_vs_", "|",colnames(contrast)[i])
  stopifnot(all.equal(names(fit$coefficients[,i]), out2[[i]]$metabolites))
  glMDPlot(fit[,i], 
           counts =noCV[,grepl(nameregrx, samps$comb)], 
           status = status,
           anno = as.data.frame(out2[[i]]$metabolites), 
           samples = samps$Sample.ID[grepl(nameregrx, samps$comb)],
           main = outputname[i], 
           groups = samps$comb[grepl(nameregrx, samps$comb)],
           transform=FALSE,
           side.main="metabolites",
           display.columns= "metabolites",
           xlab="Average log2 abundance",
           side.ylab = "log2 abundance",
           html = paste0(outputname[i],"_MDplot"),
           launch = FALSE, folder = "glimma_plots")
}

MDplots = hwriter::hwrite(colnames(contrast),
               link = paste0("./glimma_plots/",outputname,"_MDplot.html"), table=FALSE)
#----------------------------------
if(!file.exists("spreadsheets")) dir.create("spreadsheets")
for (i in 1:ncol(contrast)){
     if (nrow(out[[i]])>0){
  colnames(out[[i]]) <- gsub("adj.P.Val", "FDR",colnames(out[[i]]))
  out[[i]]$metabolites = row.names(out[[i]])
  out[[i]]= dplyr::left_join(out[[i]],    
                            feature_info[,c("metabolites","mzmed","rtmed")],by="metabolites")
  out[[i]]= select(out[[i]], metabolites, everything())
  write.xlsx(out[[i]], paste0("spreadsheets/",datamode,"_",outputname[i],".xlsx"), row.names=F)
     }
}
#------------------------------------------
# heatmap for interaction
source("heatmaps_set.R",local = T)
if (nrow(out[[ncol(contrast)]])>0){
    heatmaps_set(input = noCV,selected_metabolite = out[[ncol(contrast)]]$metabolites,
             selected_sample = samps,title = paste0(datamode,"_treatment_trait_interactions"))
}
#---------------------------------------------
# ggplot2's stat_ellipse() is using car::dataEllipse(), which draws data ellipse not confidence ellipse.
#github.com/tidyverse/ggplot2/issues/2776
mycol2=c("#87CEFA", "#104E8B","#EEC900", "#CD6600")
tiff(paste0("plots/",datamode,"_noCV_2DPCA_Confidence95_ellipse.tiff"), units="in",width = 8, height=8,res =600, compression ="lzw")
ggpubr::ggscatter(data=noCV_pca,x="PC1",y ="PC2", 
                  xlab = paste0("PC1 (", pve[1],"%)"),
                  ylab = paste0("PC2 (", pve[2],"%)"), color="comb",
                  palette=mycol2,ellipse=TRUE, mean.point=FALSE,star.plot = TRUE,
                  ellipse.type = "confidence", ellipse.level = 0.95,show.legend.text=F,
                  ellipse.alpha=0.1)+
  labs(title = paste0("PCA plot on removed Batch & Sampleweight",datamode,"data - 95% confidence Ellipse"), colour = NULL,fill=NULL)+
  theme(axis.text = element_text(size = 12, face = "bold"))
dev.off()
#---------------------------------------------
tiff(paste0("plots/",datamode,"_noCV_2DPCA_Data_ellipse.tiff"), units="in",     
     width = 8, height=8,res =600, compression ="lzw")
ggplot(data= noCV_pca, aes(x=PC1, y=PC2, color = comb))+
  geom_point()+
  scale_colour_manual(values=mycol2) +  
  stat_ellipse(type = "t", level = 0.5)+
  labs(title = paste0("PCA plot on removed Batch & Sampleweight",datamode,"data - 0.5data Ellipse"),
       x = paste0("PC1 (", pve[1],"%)"),
       y = paste0("PC2 (", pve[2],"%)"),
       colour= NULL)+
  ggpubr::theme_pubr()+
  ggpubr::labs_pubr()
dev.off()
#---------------------------------------------
###Venn
if(!file.exists("venns")) dir.create("venns")
collist = list(c(1,2),c(3,4))

venns <- makeVenn(fit, contrast= contrast[1:4,],
                  design= design[,1:4], collist = collist,
                  p.value=fdrcutoff, method="both", affy=FALSE,
                  probecol="metabolites")

vennlink = vennPage(venns, pagename = "venn", shift.title = T,
                    pagetitle = paste0("Venn Diagram ", datamode))
#------------------------------------------
names(out2) = outputname
all = do.call(cbind, out2)
colnames(all) = gsub("adj.P.Val|adj.P.value", "FDR", colnames(all))
colnames(all)[1] = sapply(strsplit(colnames(all)[1],"\\."),"[",2)
all = select(all, -ends_with(".AveExpr"))
all = select(all, -ends_with(".metabolites"))
all=dplyr::left_join(all,select(feature_info,  metabolites,mzmed, rtmed),by= "metabolites")
all.equal(all$metabolites,row.names(noCV))
all = cbind(all,noCV)
write.xlsx(all, paste0("Ben_Harrison_",datamode,"_analysis_all.xlsx"), row.names=F)
#---------------------------------------------
if(!file.exists("mummichog_input")) dir.create("mummichog_input")
#output for mummichog
out3=out2
for (i in 1:length(out3)){
  out3[[i]]= dplyr::left_join(out3[[i]], feature_info,by="metabolites") %>%
            select(., mzmed, rtmed, adj.P.Val,t)
  print(colnames(out3[[i]]) )
    if (nrow(out[[i]])> 0) {
     stopifnot(sum(out3[[i]]$adj.P.Val<fdrcutoff)==nrow(out[[i]]))
     write.table(out3[[i]],paste0("mummichog_input/",datamode,"_",outputname[i],".txt"),
                 sep="\t",quote = F,row.names = F)
    }  
}
#------------------------------------------
save.image(paste0("Ben_Harrison_pid2005_",datamode,".Rdata"))
#------------------------------------------
if(!file.exists("mummichog_output")) dir.create("mummichog_output")
setwd("mummichog_output")
(mummichog_inputfiles = list.files("../mummichog_input",full.names = T))
(mummichog_OUTfiles = list.files("../mummichog_input",full.names = F))
mummichog_OUTfiles = gsub(".txt","",mummichog_OUTfiles)

#Single file format (all features included)
for (i in seq_along(mummichog_inputfiles)){
cat(paste0("~/miniconda2/envs/mummichog/bin/mummichog1 -f ",mummichog_inputfiles[i],
           " -n fly",
           " -o ",mummichog_OUTfiles[i],
           " -m ",datamode, # mode
           " -c ",fdrcutoff, #FDR cutoff,
           " -p 100", "\n"#num of permutation
           ))
}
setwd("..")

```


