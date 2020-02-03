---
  title: "Promislow_pid1652_negative_metabolite_analysis"
output:
  BiocStyle::html_document:
  toc: true
css: style.css
bibliography: Statistical_analysis.bib
---
  
  ```{r setup, echo=FALSE, results="hide", warning=FALSE}
suppressPackageStartupMessages({
  library(knitr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(limma)
  library(threejs)
  library(htmlwidgets)
  library(xtable)
  library(ReportingTools)
  library(plotly)
  library(openxlsx)
  library(impute)
  library(missForest)
  library(DT)
  library(VIM)
  library(caret)
  library(VennDiagram)
  library(gridExtra)
  library(pROC)
  library(grid)
  library(ComplexHeatmap)
})

options(digits = 3)
```


```{r, include=FALSE}
# set column class
setwd("/Users/.. ") # Add file path to workign directory

dir()
data.neg.orig <- read.csv("../XCMS-Negdata-93samples-Promislow.csv")
# no duplicated names
length(unique(data.neg.orig$name))
# value 0 is missing values
data.neg.orig[data.neg.orig==0]=NA
names(data.neg.orig) <- gsub(".mzdata", "", names(data.neg.orig))
names(data.neg.orig) <- gsub("_", "", names(data.neg.orig))
# one sample named "neg46.2"; not sure what dose'.2' means;
names(data.neg.orig) <- gsub("neg46.2","neg46", names(data.neg.orig))
#------------------------------------

# read in sample information.
cc <- c("character","character","integer","factor","facotr","fa")
sample.info = read.csv('../sample.info.csv',header=TRUE)
sample.info =sample.info[,-2]
sample.info[,1] <- as.character(sample.info[,1])
sample.info$treatment <- as.factor(ifelse(sample.info$treatment =="B-", "Ctr", "H2O2"))
sample.info$trait <- as.factor(ifelse(sample.info$lifespan >= 80, "long_live", "short_live"))
names(sample.info)[1] <- gsub("Label", ".sample", names(sample.info)[1])

# Some summary tables
table(sample.info$treatment)
table(sample.info$trait)
table(sample.info$treatment,sample.info$trait)

# 16 lines
length(unique(sample.info$line))
table(sample.info$line,sample.info$trait) 
#----------------------------------------
# EDA
# Use density plot to compare overall distribution across sample.
mycol <- c("#8B8378", "#EE9A00") 
if(!dir.exists("plots")) dir.create("plots")

pdf("./plots/lifespan.pdf")
hist(sample.info$lifespan,breaks=40, col="grey")
dev.off()

pdf("./plots/density_plot_neg.pdf")
plot(density(log2(data.neg.orig[,30]), na.rm=TRUE), col=mycol[1], lty=2,   main = "Negative sample", ylim=c(0,0.25)) 
apply(data.neg.orig[,31:122], 2, function(x) lines(density(log2(x), na.rm=TRUE), col=mycol[1], lty=2))
apply(data.neg.orig[,21:29], 2, function(x) lines(density(log2(x), na.rm=TRUE), col=mycol[2]))
legend("topright", col=mycol, legend = c("Samples", "QC"), lty =2)
dev.off()
#-----------------------------------
# BOX plot
boxplot(log2(data.neg.orig[,21:29]),las=2,col=mycol[1],ylab="log2_mz", main="Negative sample QC")
boxplot(log2(data.neg.orig[,30:122]),las=2,col=mycol[2],ylab="log2_mz", main="Negative sample")

data.neg.orig.long <- tidyr::gather(data.neg.orig[,-c(1,3:20)],neg.sample,mz,QCneg6:neg93)

data.neg.orig.long <- left_join(data.neg.orig.long,select(sample.info, c(neg.sample,treatment ,trait)), by= "neg.sample")
data.neg.orig.long$treatment = as.character(data.neg.orig.long$treatment)
data.neg.orig.long$treatment[is.na(data.neg.orig.long$treatment)] <- "QC"

mycol2 <- c(mycol,"#483D8B")
neg.box <- ggplot(data.neg.orig.long, aes(x=neg.sample, y=log2(mz), fill=factor(treatment)))+
  geom_boxplot()+ 
  scale_fill_manual(values=mycol2)+ 
  theme(legend.position = "top") +
  labs(title = "Negative samples", y = 
         "log2_mz", fill="Treatment") + 
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(angle = 90))

pdf("plots/negsam.boxplots.pdf", width = 18,height = 15)
neg.box
dev.off()

rm(data.neg.orig.long)
#---------------------------------------------
pdf("./plots/Neg.sample.hist.pdf")
par(mfrow = c(2,2))
hist(colSums(log2(data.neg.orig[,30:122]),na.rm = T),breaks = 40, col="light blue",  main = "Sum of log2(mz) by sample")
hist(colMeans(log2(data.neg.orig[,30:122]),na.rm = T),breaks = 40, col="yellow",main = "Mean of log2(mz) by sample")
hist(rowSums(log2(data.neg.orig[,30:122]),na.rm = T),breaks = 40, col="light blue", main = "Sum of log2(mz) by feature")
hist(rowMeans(log2(data.neg.orig[,30:122]),na.rm = T),breaks = 40,col="yellow", main = "Mean of log2(mz) by feature")
hist(colMeans(is.na(data.neg.orig[,30:122])),breaks = 40, 
     main = "Percentage of missings by sample")
hist(rowMeans(is.na(data.neg.orig[,30:122])),breaks = 40,
     main = "Percentage of missings by feature")
dev.off()
#--------------------------------------

# subset only neg data
data.neg <- dplyr::select(data.neg.orig, c(name, neg5:neg93))
# RE-ORDER the column, to make it match with sample.info
re.order <- colnames(data.neg[-1])[match(sample.info$neg.sample,colnames(data.neg[-1]) )]
re.order[!is.na(re.order)] -> re.order
data.neg <- data.neg[ ,c("name", re.order)]

#-----------------------------------
# long format
data.neg.long <- data.neg %>% 
  tidyr::gather(.,neg.sample,mz,neg1:neg93)

data.neg.long <- left_join(select(sample.info, c(neg.sample,treatment,line.weight, trait)), data.neg.long, by= "neg.sample")
#----------------------------------
box_bytreatment <- ggplot(data = data.neg.long, 
                          aes(x=trait, y= log2(mz),fill=factor(treatment))) +
  geom_boxplot(width = 0.5) + 
  scale_fill_manual(values=mycol)+ 
  theme(legend.position = "top") +
  labs(title = "Negative samples", 
       y = "log2_mz", fill="Treatment") 

pdf("./plots/box_bytreatment.pdf")
box_bytreatment
dev.off()


#--------------------------------------------
# PCA
neg.matrix <- as.matrix(data.neg[,-1])
row.names(neg.matrix) <- data.neg$name

#2417 features have no missing data
nrow(na.omit(neg.matrix)) 
# use complete data for PCA
pca.neg <- prcomp(t(na.omit(log2(neg.matrix))),scale. = TRUE)
pca.neg <- as.data.frame(pca.neg$x[,1:3])

pca.neg <- merge(pca.neg, sample.info[,c("neg.sample","treatment","trait","sampleInfo","line")], by.x="row.names", by.y="neg.sample")

pca.neg.plot1 <- ggplot(pca.neg, aes(x=PC1,y=PC2))+
  geom_point(aes(color=treatment,
                 text = sampleInfo))+ 
  labs(title ="Neg PCA plot by treatment",
       x="PC 1", y="PC 2", colour = "")+
  scale_colour_manual(values=mycol)

ggplotly(pca.neg.plot1)

setwd("./plots")
saveWidget(scatterplot3js(pca.neg[,2:4], 
                          labels = paste(pca.neg$sampleInfo,pca.neg$trait), 
                          color=mycol[as.numeric(pca.neg$treatment)],
                          grid = FALSE, renderer ="canvas"),
           "pca.neg.by.treatment.html")
setwd("..")
#--------------------
pca.neg.plot2 <- ggplot(pca.neg, aes(x=PC1,y=PC2))+
  geom_point(aes(color=trait,
                 text = sampleInfo))+ 
  labs(title ="Neg PCA plot by trait",
       x="PC 1", y="PC 2", colour = "")

ggplotly(pca.neg.plot2)

pdf("./plots/PCA_bytraint_nomissing.pdf")
pca.neg.plot2
dev.off()


# 3D PCA
mycol3 = c("#FF4040", "#00CDCD")
setwd("./plots")
saveWidget(scatterplot3js(pca.neg[,2:4], 
                          labels = paste(pca.neg$sampleInfo,pca.neg$trait),
                          color=mycol3[as.numeric(pca.neg$trait)],
                          grid = FALSE, renderer ="canvas"),
           "pca.neg.by.traits.html")

setwd("..")

#----------------------------
# add PCA plot color by trait and treatment.
pca.neg$combine <- paste(pca.neg$trait, pca.neg$treatment,sep="_")
pca.neg$combine <- as.factor(pca.neg$combine)
mycol4 = c("#EEC900", "#CD6600", "#87CEFA", "#104E8B")

pca.neg.plot3 <- ggplot(pca.neg, aes(x=PC1,y=PC2))+
  geom_point(aes(color=combine,
                 text = sampleInfo))+ 
  labs(title ="Neg PCA plot by trait and treatment",
       x="PC 1", y="PC 2", colour = "")+
  scale_colour_manual(values=mycol4)

setwd("./plots")
saveWidget(ggplotly(pca.neg.plot3),"2D.pca.neg.by.traits.treatment.html")

saveWidget(scatterplot3js(pca.neg[,2:4], 
                          labels = paste(pca.neg$combine),
                          color=mycol4[as.numeric(pca.neg$combine)],
                          grid = FALSE, renderer ="canvas"),
           "pca.neg.by.traits.treatment.html")

setwd("..")

#------------------------------
# Request from Ben Harrison:
# Would you mind modifying the PC plots that you made from the analysis of the negative mode data so that the points are colored by line and that the treatment/control status is depicted as closed/open symbols?
#----------------------------
write.csv(pca.neg,"pca.neg.csv", quote = F,row.names = F)
pca.neg.plot.5 <- ggplot(pca.neg, aes(x=PC1,y=PC2,color=line,shape=treatment))+
  geom_point(aes(text =sampleInfo))+ 
  labs(title ="Neg PCA plot by line and treatment",
       x="PC 1", y="PC 2", colour = "")+
  scale_shape_manual(values =c(1,19))


setwd("./plots")
saveWidget(ggplotly(pca.neg.plot.5),"2D.pca.neg.by.line.treatment.html")
setwd("..")
#----------------------------------
# Missing data
# overall percent of missing data (total=6028 missing data)
mean(is.na(data.neg[,-1]))
# percent of features have < 5% missing data
mean(rowMeans(is.na(data.neg[,-1]))*100  < 5)
# percent of features have < 10% missing data
mean(rowMeans(is.na(data.neg[,-1]))*100  < 10)
#from 92.1% to 95.5%, not gain much from 5 to 10% cutoff

# percent of missingness by feature
table(round(rowMeans(is.na(data.neg[,-1]))*100,digits = 0))
# percent of missingness by sample
table(round(colMeans(is.na(data.neg[,-1]))*100,digits = 0))

#-------------------------------
# contingency table for missing data:
tb1 = data.neg.long %>% 
  group_by(trait,treatment) %>%
  summarise(N_missing = sum(is.na(mz)),
            percent_missing = mean(is.na(mz))*100)%>% ungroup()
tb1 = tb1[c(2,1,4,3),]
tb1 = tb1%>% group_by(trait)%>%
  mutate(label_y = cumsum(N_missing))%>% ungroup()


tb1.plot <- ggplot(tb1, aes(x=trait, y= N_missing,         
                            fill= treatment))+
  geom_bar(stat="identity") + 
  scale_fill_manual(values= mycol)+
  geom_text(aes(y=label_y*0.8, label=N_missing), 
            colour="white")+
  labs(title = "Negative sample missing data count", 
       x = "", 
       y= "Number of missingness",
       fill="Treatment") 

pdf("./plots/Freq_missing_contingency.pdf")
tb1.plot
dev.off()

#-------------------------------------
# missing by feature, across all samples  
missing_by_feature = data.neg.long %>% 
  group_by(name) %>%
  summarise(mean_log2mz = mean(log2(mz), na.rm=TRUE),
            percent_missing = mean(is.na(mz))*100)%>% ungroup()
#--------------------------------------------
# plot the percent_missing over mean_mz
pdf("./plots/Neg_percent_missing_meanMZ.pdf")
plot(missing_by_feature$mean_log2mz, missing_by_feature$percent_missing, 
     pch=19, col="light blue", xlab="mean log2(mz)", ylab="Percent of missingness", main = "Negative samples")
points(missing_by_feature$mean_log2mz[missing_by_feature$percent_missing>5],
       missing_by_feature$percent_missing[missing_by_feature$percent_missing>5], pch=19, col="orange")
legend("topright", legend=c("<=5% missingness", "> 5% missingness"),
       pch=19, col=c("light blue", "orange"))
dev.off()

#--------------------------------------------
# normalize by sample: subtract the median value of each sample (centering)
# then add the global median.
neg.matrix <- as.matrix(data.neg[,-1])
row.names(neg.matrix) <- data.neg$name
neg.matrix = log2(neg.matrix)
(median.global <- median(neg.matrix, na.rm = TRUE))
neg.matrix <- apply(neg.matrix, 2, function(x) x - median(x, na.rm = TRUE) + median.global)

# long format for after normalization
data.neg.norm <- as.data.frame(neg.matrix)
data.neg.norm$name <- row.names(data.neg.norm) 
data.neg.norm <- data.neg.norm[,c(94,1:93)]
data.neg.norm.long <- data.neg.norm %>% 
  tidyr::gather(.,neg.sample,log2_mz,neg1:neg93)

#---------------------------------
neg.box.bef <- ggplot(data.neg.long, aes(x=neg.sample, y=log2(mz), fill=mycol[1]))+
  geom_boxplot(outlier.size=0.3, outlier.shape = 21)+ 
  scale_fill_manual(values=mycol)+ 
  theme(legend.position = "top") +
  labs(title = "Before normalization", y = 
         expression(paste(log[2],"(mz)")), 
       xlab="Negative sample") + 
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(angle = 90),
        legend.position = "none")

neg.box.aft <- ggplot(data.neg.norm.long, aes(x=neg.sample, y=log2_mz, fill=mycol[2]))+
  geom_boxplot(outlier.size=0.3, outlier.shape = 21)+ 
  scale_fill_manual(values=mycol[2])+ 
  theme(legend.position = "top") +
  labs(title = "After normalization", y = 
         expression(paste(log[2],"(mz)")), 
       xlab="Negative sample") + 
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(angle = 90),
        legend.position = "none")

pdf("./plots/normalization_boxplots.pdf")
neg.box.bef
neg.box.aft
dev.off()

#plot_ly(data.neg.long, y=~log2(mz),
#        =~neg.sample, color ="black",
#        alpha=0.1,type ="box")%>%
#    layout(title ="Before normalization")

#-----------------------------
# if only keep feature has < 5% overall missingness
select_feature = missing_by_feature[missing_by_feature$percent_missing<5,]$name
neg.matrix.filter = neg.matrix[row.names(neg.matrix) %in% select_feature,]
dim(neg.matrix.filter)

# double check if missingness < 5%
range(apply(neg.matrix.filter, 1, function(x) mean(is.na(x))*100))

#---------------------------------------------------
# Missingness patterns after filetering:
pdf("./plots/missingness.pattern.lessthan5.pdf")
VIM::aggr(neg.matrix.filter,sortVars = TRUE,only.miss = TRUE)
dev.off()

# Use Imputation methods:
# 1. small value (sv)replacement:
# for every missing values replace with half of minimum value of entire dataset
sv.imp <- neg.matrix.filter
sv.imp[is.na(sv.imp)] <- min(sv.imp,na.rm=TRUE)/2

#---------------------------------------------------
# 2. Mean replacement: for every missing values replace with mean value of that feature across all samples (exclude missing in calculation)

impute.mean <- function(x) {
  z <- apply(x,1, function(x) mean(x, na.rm=TRUE))
  for (i in seq(nrow(x)))
    x[i,][is.na(x[i,])] <- z[i]
  return(x)
}

mean.imp <- neg.matrix.filter
mean.imp = impute.mean(mean.imp)

#----------------------------------------------------------------------
# 3. knn (use 'impute' package)
# use expression matrix as input data. defaulthe rowmax is 50%, > 50% missingness, will use overall mean per sample for imputed value.
# An expression matrix with metabolites in the rows, samples in the columns
knn.imp <- impute::impute.knn(neg.matrix.filter,k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)$data

# Question: should i impute based on entire data then filter or
# filter out >5% then impute.
#---------------------------------
# 4. RF (use 'missForest' package)
# variables in the columns and observations in the rows
if(!file.exists("rf.imp.Rdata")) {
  set.seed(32173)
  Sys.setenv(R_THREADS = 20)
  rf.imp <- missForest(t(neg.matrix.filter),verbose = TRUE)
  rf.OOBerror <- rf.imp$OOBerror
  rf.imp <- t(rf.imp$ximp)
  save(list = "rf.imp", file = "rf.imp.Rdata")
}else{
  load("rf.imp.Rdata")
}

#---------------------------------
# compare the PCA on non-imputed data vs imputed data,
# see any improvement on classify trait.
pca.imp <- function(x, method) {
  imp <- prcomp(t(x),scale. = TRUE)
  pve <- (imp$sdev^2)/sum(imp$sdev^2)
  imp.df <- as.data.frame(imp$x[,1:3])
  imp.df <- merge(imp.df, sample.info[,c("neg.sample","treatment","trait","sampleInfo")], by.x="row.names", by.y="neg.sample")
  
  
  imp.pcaplot <- ggplot(imp.df, aes(x=PC1,y=PC2))+
    geom_point(aes(color=trait,
                   text = paste0(imp.df$trait,imp$sampleInfo)))+ 
    labs(title =paste0("Neg PCA plot by trait- ",
                       method," imputation"),
         x="PC 1", y="PC 2", colour = "")
  
  
  imp.pcaplot.3D <- scatterplot3js(imp.df[,2:4], 
                                   labels = paste(imp.df$sampleInfo,imp.df$trait),
                                   color=mycol3[as.numeric(imp.df$trait)],
                                   grid = FALSE, renderer ="canvas")
  
  out <- list(imp,imp.df, imp.pcaplot,imp.pcaplot.3D, pve)
  names(out) <- c("pca.out","pca.df", "pca.plot","pcaplot.3D", "PVE")
  return(out)
}

noimp.pca = pca.imp(na.omit(neg.matrix.filter), method = "completed")
# 1. use sv.imp data for PCA
sv.pca = pca.imp(sv.imp, method = "SV")
# 2. mean.imp
mean.pca = pca.imp(mean.imp, method = "Mean")
# 3. KNN
knn.pca = pca.imp(knn.imp, method = "KNN")
# 4. RF
rf.pca = pca.imp(rf.imp, method ="Random Forest")

pdf("./plots/PCA_imputations.pdf")
sv.pca$pca.plot 
mean.pca$pca.plot
knn.pca$pca.plot
rf.pca$pca.plot
dev.off()


setwd("plots")
saveWidget(sv.pca$pcaplot.3D, "pca.neg.sv_imp.html")
saveWidget(mean.pca$pcaplot.3D, "pca.neg.mean_imp.html")
saveWidget(knn.pca$pcaplot.3D, "pca.neg.knn_imp.html")
saveWidget(rf.pca$pcaplot.3D, "pca.neg.rf_imp.html")
setwd("..")
#---------------------------
# compare the number of the PCs needed to reach 80% variance in complete data vs imputed data.
pve = data.frame(noimp.pca$PVE, sv.pca$PVE, mean.pca$PVE, knn.pca$PVE, rf.pca$PVE)
pve_20PCs = data.frame(PVE = colSums(pve[1:20,])*100)

#Plot the cumulative variance
col5 = c("#000080","#00BFFF", "#EE7600", "#76EEC6", "#EE2C2C")
pdf("./plots/cumsum_var_pca.pdf")
plot(cumsum(noimp.pca$PVE), pch =19, col=col5[1], main= "Cumulative variance from PCA")
points(cumsum(sv.pca$PVE), pch =18, col=col5[2])
points(cumsum(mean.pca$PVE), pch =17, col=col5[3])
points(cumsum(knn.pca$PVE), pch =16, col=col5[4])
points(cumsum(rf.pca$PVE), pch =15, col=col5[5])
legend("bottomright", legend = c("completed","SV", "mean","KNN", "RF"),
       col=col5[1:5], pch = 19:15)
dev.off()
#------------------------------
# transfer to long format for visulization
reformat <- function (x) {
  df <- as.data.frame(x)
  df$name <- row.names(df)
  df <- select(df,name,neg1:neg93)
  df <- df %>% tidyr::gather(.,neg.sample,log2_mz,neg1:neg93)
  df <- left_join(select(sample.info, c(neg.sample,treatment,line.weight, trait)), df, by= "neg.sample")
}

# non-imputed data
data.neg.noimp <- reformat(neg.matrix.filter)
impute_ind <- factor(ifelse(is.na(data.neg.noimp$log2_mz), "Yes", "No"))

all.equal(data.neg.noimp[,1:5], reformat(sv.imp)[,1:5])
all.equal(data.neg.noimp[,1:5], reformat(mean.imp)[,1:5])
all.equal(data.neg.noimp[,1:5], reformat(knn.imp)[,1:5])
all.equal(data.neg.noimp[,1:5], reformat(rf.imp)[,1:5])

# check the imputed data distribution
impute_plot <- function(imp,method) {
  tmp <- reformat(imp)
  tmp$impute_ind <- impute_ind 
  
  p <- ggplot(data=tmp, aes(x=factor(impute_ind),y=log2_mz,fill=trait))+
    geom_boxplot()+
    labs(title = paste0(method," Imputed"), x="Imputed or not")+
    theme(legend.position = "top")
  print(p)
  
}

if(!file.exists("./plots/imputedata_distri.pdf")) {
  pdf("./plots/imputedata_distri.pdf")
  impute_plot(sv.imp, method= "SV")
  impute_plot(mean.imp,method= "Mean" )
  impute_plot(knn.imp, method="KNN")
  impute_plot(rf.imp, method = "RF")
  dev.off()
}

#--------------------------
# knn imputed value change over selection of k(1 to 30)
tmp.knn = list()
for (d in 1:30){
  tmp.knn[[d]] <-impute::impute.knn(neg.matrix.filter,k = d, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)$data
}


knn_nrmse <- function(x){
  df <- reformat(x)
  df$impute_ind <- impute_ind 
  imp <- filter(df, impute_ind=="Yes")
  noimp <- filter(df, impute_ind=="No")
  
  NRMSE = sqrt((mean(noimp$log2_mz)-mean(imp$log2_mz))^2/var(noimp$log2_mz))
  return(NRMSE)
}

NRMSE = sapply(tmp.knn, knn_nrmse)
knn.df <- data.frame(k=1:30, NRMSE = NRMSE)
plot(knn.df$k, knn.df$NRMSE, type="b", xlab="K", ylab="NRMSE",
     main ="NRMSE between nonimputed vs knn imputed")

#----------------------------------------------
# transverse to each row is a sample, and each column is a predictor.
row.names(sample.info) <- sample.info$neg.sample

# transverse. add sample info.
reformat.model <- function (x) {
  df <- as.data.frame(t(x))
  df <- merge(select(sample.info, c(treatment,line.weight, trait)), df, by="row.names")
  return(df)
}

#-----------------------------
# with line.weight as a covariate (full model)
univariate_logistic <- function (x, group){
  # filter based on treatment.
  tmp = reformat.model(x)
  tmp2 = filter(tmp,treatment==paste0(group))
  
  model = list()
  Z = coef = p_values = rep(NA,(ncol(tmp2)-4))
  for (i in 1:(ncol(tmp2)-4)) {
    model[[i]] = summary(glm(relevel(trait, ref = "short_live") ~ tmp2[,(4+i)] + line.weight, data = tmp2, family=binomial))
    Z[i] = model[[i]]$coef[2,3]
    coef[i] = model[[i]]$coef[2,1]
    p_values[i] = model[[i]]$coef[2,4]
  }
  
  model_summary = data.frame(metabolite = row.names(x), Z, coef, p_values)
  model_summary$FDR <- p.adjust(model_summary$p_values,method = "BH")
  
  out <- list(model_summary, model)
  names(out) <- c("model_summary","model")
  return(out)
}

# no covariate (null model)
univariate_logistic.nocov <- function (x, group){
  # filter based on treatment.
  tmp = reformat.model(x)
  tmp2 = filter(tmp,treatment==paste0(group))
  
  model = list()
  Z = coef = p_values= rep(NA,(ncol(tmp2)-4))
  for (i in 1:(ncol(tmp2)-4)) {
    model[[i]] = summary(glm(relevel(trait, ref = "short_live") ~ tmp2[,(4+i)], data = tmp2, family=binomial))
    Z[i] = model[[i]]$coef[2,3]
    coef[i] = model[[i]]$coef[2,1]
    p_values[i] = model[[i]]$coef[2,4]
  }
  
  model_summary = data.frame(metabolite = row.names(x),Z, coef, p_values)
  model_summary$FDR <- p.adjust(model_summary$p_values,method = "BH")
  
  out <- list(model_summary, model)
  names(out) <- c("model_summary","model")
  return(out)
}

#------------------------------
# 1. Univariate logistic model for each data generated by different the imputation methods with line.weight as covariate.
# use 10% FDR  as cutoff
fdr = 0.1
#--------------
# no imputation
noimp_uni_Ctr <- univariate_logistic(neg.matrix.filter, group ="Ctr")$model_summary
noimp_uni_H2O2 <- univariate_logistic(neg.matrix.filter, group ="H2O2")$model_summary

sig_Ctr_noimp = noimp_uni_Ctr[noimp_uni_Ctr$FDR < fdr,]
sig_H2O2_noimp = noimp_uni_H2O2[noimp_uni_H2O2$FDR < fdr,]

# small value imputed
sv_uni_Ctr <- univariate_logistic(sv.imp, group ="Ctr")$model_summary
sv_uni_H2O2 <- univariate_logistic(sv.imp, group ="H2O2")$model_summary
#Warning message:
#glm.fit: fitted probabilities numerically 0 or 1 occurred
sig_Ctr_sv_imp = sv_uni_Ctr[sv_uni_Ctr$FDR < fdr,]
sig_H2O2_sv_imp = sv_uni_H2O2[sv_uni_H2O2$FDR < fdr,]

# mean value imputed
mean_uni_Ctr <- univariate_logistic(mean.imp, group ="Ctr")$model_summary
mean_uni_H2O2 <- univariate_logistic(mean.imp, group ="H2O2")$model_summary
sig_Ctr_mean_imp = mean_uni_Ctr[mean_uni_Ctr$FDR < fdr,]
sig_H2O2_mean_imp = mean_uni_H2O2[mean_uni_H2O2$FDR < fdr,]

# 10KNN imputed
knn_uni_Ctr <- univariate_logistic(knn.imp, group ="Ctr")$model_summary
knn_uni_H2O2 <- univariate_logistic(knn.imp, group ="H2O2")$model_summary
sig_Ctr_knn_imp = knn_uni_Ctr[knn_uni_Ctr$FDR < fdr,]
sig_H2O2_knn_imp = knn_uni_H2O2[knn_uni_H2O2$FDR < fdr,]

write.table(knn_uni_Ctr, "Mummichog_Results/Raw_files/knn_uni_Ctr.txt", quote = F, row.names = F)
write.table(knn_uni_H2O2, "Mummichog_Results/Raw_files/knn_uni_H2O2.txt", quote = F, row.names = F)

# RF Imputed
rf_uni_Ctr <- univariate_logistic(rf.imp, group ="Ctr")$model_summary
rf_uni_H2O2 <- univariate_logistic(rf.imp, group ="H2O2")$model_summary
sig_Ctr_rf_imp = rf_uni_Ctr[rf_uni_Ctr$FDR < fdr,]
sig_H2O2_rf_imp = rf_uni_H2O2[rf_uni_H2O2$FDR < fdr,]

# KNN and RF significant feature are the same in control group.
all.equal(sig_Ctr_knn_imp$metabolite, sig_Ctr_rf_imp$metabolite)

# 7 features in rf but not in knn; 0 rev.
setdiff(sig_H2O2_rf_imp$metabolite,sig_H2O2_knn_imp$metabolite)
setdiff(sig_H2O2_knn_imp$metabolite, sig_H2O2_rf_imp$metabolite)
#-------------------------
# output significant gene lists.
ls.ctr <- list(sig_Ctr_noimp,sig_Ctr_sv_imp,sig_Ctr_mean_imp,sig_Ctr_knn_imp,sig_Ctr_rf_imp)
names(ls.ctr) <- c("control samples, no imputation",
                   "control samples, SV imputation",
                   "control samples, Mean imputation",
                   "control samples, KNN imputation",
                   "control samples, RF imputation")


ls.H2O2 <- list(sig_H2O2_noimp,sig_H2O2_sv_imp,sig_H2O2_mean_imp,sig_H2O2_knn_imp,sig_H2O2_rf_imp)
names(ls.H2O2) <- c("H2O2 samples, no imputation",
                    "H2O2 samples, SV imputation",
                    "H2O2 samples, Mean imputation",
                    "H2O2 samples, KNN imputation",
                    "H2O2 samples, RF imputation")

# Generate spreadsheets for each contrast.
if(!file.exists("spreadsheets")) dir.create("spreadsheets")
if(!file.exists("reports")) dir.create("reports")

ctab <- lapply(1:length(ls.ctr), function(x) if(nrow(ls.ctr[[x]]) > 0){
  tmp <- ReportingTools::HTMLReport(names(ls.ctr)[x], names(ls.ctr)[x],"./reports",".")
  publish(ls.ctr[[x]], tmp)
  write.xlsx(ls.ctr[[x]], paste0("spreadsheets/",names(ls.ctr)[x], ".xlsx"), quote= FALSE,row.names= FALSE)
  finish(tmp)
  return(tmp)
})


htab <- lapply(1:length(ls.H2O2), function(x) if(nrow(ls.H2O2[[x]]) > 0){
  tmp <- ReportingTools::HTMLReport(names(ls.H2O2)[x], names(ls.H2O2)[x], "./reports",".")
  publish(ls.H2O2[[x]], tmp)
  write.xlsx(ls.H2O2[[x]], paste0("spreadsheets/",names(ls.H2O2)[x], ".xlsx"), quote= FALSE,row.names= FALSE)
  finish(tmp)
  return(tmp)
})


#---------------------------------
# Compare the results of using 4 imputation,
# none of imputation data the H2O2 group had sig Features.
# use Venn diagrams to compare the results in ctr group.
venncol= c("#8DEEEE", "#FFB90F","#483D8B", "#838B8B")

# Control
venn.plot.c <- venn.diagram(list(sig_Ctr_noimp$metabolite,
                                 sig_Ctr_sv_imp$metabolite,
                                 sig_Ctr_mean_imp$metabolite,
                                 sig_Ctr_knn_imp$metabolite), NULL,
                            width = 3200, fill=venncol, 
                            alpha=c(0.5,0.5,0.5,0.5), 
                            cex = 1.5, cat.fontface=8,  
                            category.names=c("No imputation", "SV imputed", "Mean imputed", "KNN imputed"),main="Control group")

pdf("./plots/Venn_Ctr_sig_features.pdf")
grid.draw(venn.plot.c)
dev.off()

# difference between imputed vs no imputed results.
# 6 feature in noimputed but not in sv imputed
setdiff(sig_Ctr_noimp$metabolite,sig_Ctr_sv_imp$metabolite)
# 2 feature in mean, knn or rf but not in no-imputed.
setdiff(sig_Ctr_mean_imp$metabolite,sig_Ctr_noimp$metabolite)
setdiff(sig_Ctr_knn_imp$metabolite,sig_Ctr_noimp$metabolite)
setdiff(sig_Ctr_rf_imp$metabolite,sig_Ctr_noimp$metabolite)

#-----------------------------
# H2O2
venn.plot.h1 <- venn.diagram(list(sig_H2O2_noimp$metabolite,sig_H2O2_sv_imp$metabolite,sig_H2O2_knn_imp$metabolite,sig_H2O2_rf_imp$metabolite), NULL,
                             width=3200, fill=venncol,main=expression('H'[2]*'O'[2]*' group'), alpha=c(0.5,0.5,0.5,0.5), cex = 1.5, cat.fontface=8,  category.names=c("NO imputation", "SV imputed", "KNN imputed","RF imputed"))
pdf("./plots/Venn_H2O2_sig_features_1.pdf")
grid.draw(venn.plot.h1)
dev.off()

venn.plot.h2 <- venn.diagram(list(sig_H2O2_noimp$metabolite,sig_H2O2_mean_imp$metabolite,sig_H2O2_knn_imp$metabolite,sig_H2O2_rf_imp$metabolite), NULL, width=3200, fill=venncol, alpha=c(0.5,0.5,0.5,0.5), cex = 1.5, cat.fontface=8,  main=expression('H'[2]*'O'[2]*' group'), category.names=c("NO imputation", "Mean imputed", "KNN imputed","RF imputed"))

pdf("./plots/Venn_H2O2_sig_features_2.pdf")
grid.draw(venn.plot.h2)
dev.off()
#-------------------------------
# add Venn diagrams comparing control and H2O2 significant features.

venn.plot.join1 <- venn.diagram(list(sig_Ctr_noimp$metabolite,
                                     sig_H2O2_noimp$metabolite), NULL, width=3200, fill=venncol[3:4], alpha=c(0.5,0.5), cex = 1.5, cat.fontface=8,  main="No imputation", category.names=c("Control", "H2O2"))

venn.plot.join2 <- venn.diagram(list(sig_Ctr_sv_imp$metabolite,
                                     sig_H2O2_sv_imp$metabolite), NULL, width=3200, fill=venncol[3:4], alpha=c(0.5,0.5), cex = 1.5, cat.fontface=8,  main="SV imputation", category.names=c("Control", "H2O2"))

venn.plot.join3 <- venn.diagram(list(sig_Ctr_mean_imp$metabolite,
                                     sig_H2O2_mean_imp$metabolite), NULL, width=3200, fill=venncol[3:4], alpha=c(0.5,0.5), cex = 1.5, cat.fontface=8,  main="Mean imputation", category.names=c("Control", "H2O2"))


venn.plot.join4 <- venn.diagram(list(sig_Ctr_knn_imp$metabolite,
                                     sig_H2O2_knn_imp$metabolite), NULL, width=3200, fill=venncol[3:4], alpha=c(0.5,0.5), cex = 1.5, cat.fontface=8,  main="KNN imputation", category.names=c("Control", "H2O2"))


venn.plot.join5 <- venn.diagram(list(sig_Ctr_rf_imp$metabolite,
                                     sig_H2O2_rf_imp$metabolite), NULL, width=3200, fill=venncol[3:4], alpha=c(0.5,0.5), cex = 1.5, cat.fontface=8,  main="RF imputation", category.names=c("Control", "H2O2"))

plot.new()
pdf("plots/Venn_uni_control_h2o2.pdf", width = 14, height = 12)
# setup layout
gl <- grid.layout(nrow=2, ncol=2)
# setup viewports
vp.1 <- viewport(layout.pos.col=1, layout.pos.row=1) 
vp.2 <- viewport(layout.pos.col=2, layout.pos.row=1) 
vp.3 <- viewport(layout.pos.col=1, layout.pos.row=2) 
vp.4 <- viewport(layout.pos.col=2, layout.pos.row=2) 

# init layout
pushViewport(viewport(layout=gl))
# access the first position
pushViewport(vp.1)
# start new base graphics in first viewport
par(new=TRUE, fig=gridBase::gridFIG())
grid.draw(venn.plot.join2)
# done with the first viewport
popViewport()
# move to the next viewport
pushViewport(vp.2)
grid.draw(venn.plot.join3)
# done with this viewport
popViewport(1)

pushViewport(vp.3)
grid.draw(venn.plot.join4)
# done with the first viewport
popViewport(1)

# move to the next viewport
pushViewport(vp.4)
grid.draw(venn.plot.join5)
# done with this viewport
popViewport(1)
dev.off()

#----------------
# difference between imputed vs no imputed results.
# 8 feature in noimputed but not in sv imputed; 3 rev
setdiff(sig_H2O2_noimp$metabolite,sig_H2O2_sv_imp$metabolite)
setdiff(sig_H2O2_sv_imp$metabolite,sig_H2O2_noimp$metabolite)
# 1 feature in mean but not in no-imputed;1 rev
setdiff(sig_H2O2_mean_imp$metabolite,sig_H2O2_noimp$metabolite)
setdiff(sig_H2O2_noimp$metabolite,sig_H2O2_mean_imp$metabolite)

# 5 feature in KNN but not in no-imputed;2 rev
setdiff(sig_H2O2_knn_imp$metabolite,sig_H2O2_noimp$metabolite)
setdiff(sig_H2O2_noimp$metabolite,sig_H2O2_knn_imp$metabolite)
# 11 feature in RF but not in no-imputed;1 rev
setdiff(sig_H2O2_rf_imp$metabolite,sig_H2O2_noimp$metabolite)
setdiff(sig_H2O2_noimp$metabolite,sig_H2O2_rf_imp$metabolite)

#----------------------------------
# no covariate in model.
# no imputation
noimp_uni_Ctr.nocov <- univariate_logistic.nocov(neg.matrix.filter, group ="Ctr")$model_summary
noimp_uni_H2O2.nocov <- univariate_logistic.nocov(neg.matrix.filter, group ="H2O2")$model_summary
sum(noimp_uni_Ctr.nocov$FDR < fdr);sum(noimp_uni_H2O2.nocov$FDR < fdr)

# small value imputed
sv_uni_Ctr.nocov <- univariate_logistic.nocov(sv.imp, group ="Ctr")$model_summary
sv_uni_H2O2.nocov <- univariate_logistic.nocov(sv.imp, group ="H2O2")$model_summary
#Warning message:
#glm.fit: fitted probabilities numerically 0 or 1 occurred
sum(sv_uni_Ctr.nocov$FDR < fdr);sum(sv_uni_H2O2.nocov$FDR < fdr)


# mean value imputed
mean_uni_Ctr.nocov <- univariate_logistic.nocov(mean.imp, group ="Ctr")$model_summary
mean_uni_H2O2.nocov <- univariate_logistic.nocov(mean.imp, group ="H2O2")$model_summary
sum(mean_uni_Ctr.nocov$FDR < fdr);sum(mean_uni_H2O2.nocov$FDR < fdr)

# 10KNN imputed
knn_uni_Ctr.nocov <- univariate_logistic.nocov(knn.imp, group ="Ctr")$model_summary
knn_uni_H2O2.nocov <- univariate_logistic.nocov(knn.imp, group ="H2O2")$model_summary
sum(knn_uni_Ctr.nocov$FDR < fdr);sum(knn_uni_H2O2.nocov$FDR < fdr)

# RF Imputed
rf_uni_Ctr.nocov <- univariate_logistic.nocov(rf.imp, group ="Ctr")$model_summary
rf_uni_H2O2.nocov <- univariate_logistic.nocov(rf.imp, group ="H2O2")$model_summary
sum(rf_uni_Ctr.nocov$FDR < fdr);sum(rf_uni_H2O2.nocov$FDR < fdr)

#----------------------
# check if line.weight is a confounder
# compare the coef of features in full vs null models see if change >10%.
confounder.df <- function(full_model, null_model) {
  full_model <- select(full_model, coef)
  null_model <- select(null_model, coef)
  combined <- merge(full_model,null_model, by ="row.names")
  colnames(combined)[2:3] <- c("coef_full", "coef_null")
  combined$perc_diff <- abs((combined$coef_full - combined$coef_null))/combined$coef_null
  # IF coeff change > 10% called confounder
  combined$confounder <- ifelse(combined$perc_diff > 0.1 ,"YES","NO") 
  print(table(combined$confounder))
  return(combined)
}
#---------------------------------------
confounder.noimp_ctr <- confounder.df(noimp_uni_Ctr,noimp_uni_Ctr.nocov)
confounder.noimp_H2O2 <- confounder.df(noimp_uni_H2O2,noimp_uni_H2O2.nocov)

confounder.sv_ctr <- confounder.df(sv_uni_Ctr,sv_uni_Ctr.nocov)
confounder.sv_H2O2 <- confounder.df(sv_uni_H2O2,sv_uni_H2O2.nocov)

confounder.mean_ctr <- confounder.df(mean_uni_Ctr,mean_uni_Ctr.nocov)
confounder.mean_H2O2 <- confounder.df(mean_uni_H2O2,mean_uni_H2O2.nocov)

confounder.knn_ctr <- confounder.df(knn_uni_Ctr,knn_uni_Ctr.nocov)
confounder.knn_H2O2 <- confounder.df(knn_uni_H2O2,knn_uni_H2O2.nocov)

confounder.rf_ctr <- confounder.df(rf_uni_Ctr,rf_uni_Ctr.nocov)
confounder.rf_H2O2 <- confounder.df(rf_uni_H2O2,rf_uni_H2O2.nocov)
#---------------------
ls_c <- list(confounder.noimp_ctr,confounder.sv_ctr,confounder.mean_ctr,
             confounder.knn_ctr,confounder.rf_ctr)
# around 42% are confounders in control
num_confounder_Ctr <- sapply(ls_c, function(x) nrow(filter(x,confounder=="YES")))*100/sapply(ls_c, function(x) nrow(filter(x,confounder=="NO")))

ls_h <- list(confounder.noimp_H2O2,confounder.sv_H2O2,confounder.mean_H2O2,
             confounder.knn_H2O2,confounder.rf_H2O2)
# around 32% are confounders in control
num_confounder_H2O2 <-sapply(ls_h, function(x) nrow(filter(x,confounder=="YES")))*100/sapply(ls_h, function(x) nrow(filter(x,confounder=="NO")))
#------------------------------------------

# conclusion: line.weight is considered as confounder in ~25-30% univariate models.
# so keep it there.

#--------------------------------------------
# Test the interaction between treatment and trait using linear regression for each feature, adjusting for line.weight.
univariate_linear.inter <- function (x){
  
  tmp = reformat.model(x)
  tmp$trait = relevel(tmp$trait, ref="short_live")
  model = model2 = residual = list()
  trait.coef = trait.t.stat = trait.p.values = treat.coef = treat.t.stat = treat.p.values = int.coef = int.t.stat = int.p.values = rep(NA,(ncol(tmp)-4))
  
  
  for (i in 1:(ncol(tmp)-4)) {
    # use interaction term between mz and treatment
    model[[i]] = summary(lm(tmp[,(4+i)] ~ trait*treatment + line.weight, data = tmp))
    model2[[i]] = summary(lm(tmp[,(4+i)] ~ line.weight, data = tmp))
    residual[[i]] = residuals(model2[[i]])
    
    trait.coef[i] = model[[i]]$coef[2,1]
    trait.t.stat[i] = model[[i]]$coef[2,3]
    trait.p.values[i] = model[[i]]$coef[2,4]
    
    treat.coef[i] =model[[i]]$coef[3,1]
    treat.t.stat[i] =model[[i]]$coef[3,3]
    treat.p.values[i] = model[[i]]$coef[3,4]
    
    int.coef[i] =model[[i]]$coef[5,1]
    int.t.stat[i] = model[[i]]$coef[5,3]
    int.p.values[i] = model[[i]]$coef[5,4]
  }
  names(residual) <- rownames(x)
  residual = do.call(rbind, residual)
  colnames(residual) <- colnames(x)
  
  model_summary = data.frame(metabolite = row.names(x), trait.coef,trait.t.stat,trait.p.values,treat.coef,treat.t.stat,treat.p.values,int.coef,int.t.stat,int.p.values)
  row.names(model_summary) = row.names(x)
  model_summary$trait.FDR <- p.adjust(model_summary$trait.p.values,method = "BH")
  model_summary$treat.FDR <- p.adjust(model_summary$treat.p.values,method = "BH")
  model_summary$int.FDR <- p.adjust(model_summary$int.p.values,method = "BH")
  
  out <- list(model_summary, model,residual)
  names(out) <- c("model_summary","model","residual")
  return(out)
}


#-------------------------
# no imputation
noimp_uni_int <- univariate_linear.inter(neg.matrix.filter)$model_summary
sum(noimp_uni_int$int.FDR < fdr)
sig_noimp_uni_int <- noimp_uni_int[noimp_uni_int$int.FDR < fdr,]

# sv imputation
sv_imp_uni_int <- univariate_linear.inter(sv.imp)$model_summary
sum(sv_imp_uni_int$int.FDR < fdr)
sig_sv_imp_uni_int <- sv_imp_uni_int[sv_imp_uni_int$int.FDR < fdr,]

# mean imputation
mean_imp_uni_int <- univariate_linear.inter(mean.imp)$model_summary
sum(mean_imp_uni_int$int.FDR < fdr)
sig_mean_imp_uni_int <- mean_imp_uni_int[mean_imp_uni_int$int.FDR < fdr,]

# 10KNN imputed
knn_imp_uni_int <- univariate_linear.inter(knn.imp)$model_summary
sum(knn_imp_uni_int$int.FDR < fdr)
sig_knn_imp_uni_int <- knn_imp_uni_int[knn_imp_uni_int$int.FDR < fdr,]

write.table(knn_imp_uni_int, "Mummichog_Results/Raw_files/knn_imp_uni_int.txt", sep="\t", row.names = FALSE,quote = FALSE)


# RF Imputed
rf_imp_uni_int <- univariate_linear.inter(rf.imp)$model_summary
sum(rf_imp_uni_int$int.FDR < fdr)
sig_rf_imp_uni_int <- rf_imp_uni_int[rf_imp_uni_int$int.FDR < fdr,]

#-------------------------
# Venn diagrams
venn.plot.int <- venn.diagram(list(sig_noimp_uni_int$metabolite,sig_mean_imp_uni_int$metabolite,sig_knn_imp_uni_int$metabolite,sig_rf_imp_uni_int$metabolite), NULL, width=3200, fill=venncol, alpha=c(0.5,0.5,0.5,0.5), cex = 1.5, cat.fontface=8,  main= "Linear Regression with Trait x Treatment interaction", category.names=c("NO imputation", "Mean imputed", "KNN imputed","RF imputed"))

pdf("./plots/Venn_lm.interaction_sig_features.pdf")
grid.draw(venn.plot.int)
dev.off()
#------------------------------
# Heatmap of the significant features (use noimp_uni_int)
# long format for after normalization
sig.uni.int.long <- data.frame(knn.imp, name = row.names(knn.imp)) %>% 
  tidyr::gather(.,neg.sample,log2_mz,neg1:neg93)%>%
  filter(name %in% sig_knn_imp_uni_int$metabolite)
sig.uni.int.long <- left_join(sig.uni.int.long, sample.info[,c("neg.sample","treatment","trait","lifespan")], by ="neg.sample")
# Scatter plot of metabolite abundance vs life span, color by treatment for each feature.
splot <- ggplot(sig.uni.int.long, aes(x= trait, y =log2_mz, col = treatment))+  geom_point(size=0.8, alpha=0.7,position=position_dodge(width=0.3)) + facet_wrap(~ name)+ labs(title ="Scatter plot by treatment",
                                                                                                                                                                              x="Trait", y="log2 (mz)", colour = "")+
  scale_colour_manual(values=mycol)+theme(legend.position = "top", legend.direction = "horizontal",axis.text.x = element_text(angle = 90))

pdf("plots/scatter.plot.inter.pdf", width =15, height =18)
splot
dev.off()

#---------------------
# heatmap of significant features
# Plot with residuals removing effect of line.weight.
knn.residual = univariate_linear.inter(knn.imp)$residual
sig.uni.int.long2 <- data.frame(knn.residual, name = row.names(knn.residual)) %>% 
  tidyr::gather(.,neg.sample, residual,neg1:neg93)%>%
  filter(name %in% sig_knn_imp_uni_int$metabolite)
sig.uni.int.long2 <- left_join(sig.uni.int.long2, sample.info[,c("neg.sample","treatment","trait","lifespan")], by ="neg.sample")
sig.uni.int.long2$comb <- paste(sig.uni.int.long2$trait, sig.uni.int.long2$treatment,sep="_")
sig.uni.int.long2 = select(sig.uni.int.long2, name, residual, comb) %>%
  group_by(name, comb) %>% 
  mutate(mean_residual = mean(residual, na.rm=T)) %>% ungroup()

sig.uni.int.long2 = select(sig.uni.int.long2, name, comb, mean_residual)
sig.uni.int.long2 = sig.uni.int.long2[!duplicated(sig.uni.int.long2),]

sig.uni.int = tidyr::spread(sig.uni.int.long2, key= comb, value= mean_residual)
sig.uni.int = as.data.frame(sig.uni.int)
rownames(sig.uni.int) = sig.uni.int$name
sig.uni.int = as.matrix(sig.uni.int[,-1])

#-----------------------
# center and scale sig.uni.int row-wise
sig.uni.int2 = t(apply(sig.uni.int,1, function(x) (x-mean(x))/sd(x)))

setwd("plots")
saveWidget(d3heatmap::d3heatmap(sig.uni.int2,dendrogram="row",cexCol =0.3,
                                cexRow = 0.7, color=rev(RColorBrewer::brewer.pal(5,"RdBu"))),"heatmap_sig_knn_inter_scaled.html")
setwd("..")

heatmap_sig_int = Heatmap(sig.uni.int2, col=rev(RColorBrewer::brewer.pal(5,"RdBu")),
                          cluster_rows = TRUE, cluster_columns = FALSE,row_names_gp = gpar(cex=0.3),
                          column_names_gp = gpar(cex=0.8),
                          column_title = "Metabolites with significant interaction between trait and treatment",name ="z-score")

pdf("plots/heatmap_sig_int_scaled.pdf")
heatmap_sig_int
dev.off()


#------------------------
# output significant interaction features.
ls.int <- list(sig_noimp_uni_int,sig_sv_imp_uni_int,sig_mean_imp_uni_int,sig_knn_imp_uni_int,sig_rf_imp_uni_int)

#eset = list(neg.matrix.filter, sv.imp, mean.imp, knn.imp, rf.imp)
#rmbatch <- lapply(eset, function(x) removeBatchEffect(x,batch = sample.info$line.weight))

names(ls.int) <- c("Trait & Treatment interaction, No imputation",
                   "Trait & Treatment interaction, SV imputation",
                   "Trait & Treatment interaction, Mean imputation",
                   "Trait & Treatment interaction, KNN imputation",
                   "Trait & Treatment interaction, RF imputation")

inter.tab <- lapply(1:length(ls.int), function(x) if(nrow(ls.int[[x]]) > 0){
  tmp <- ReportingTools::HTMLReport(names(ls.int)[x], names(ls.int)[x],"./reports",".")
  publish(ls.int[[x]][,c(1,8:10,13)], tmp)
  write.xlsx(ls.int[[x]][,c(1,8:10,13)], paste0("spreadsheets/",names(ls.int)[x], ".xlsx"), quote= FALSE,row.names= FALSE)
  finish(tmp)
  return(tmp)
})

#-----------------------------------------
# 2. Multivariate Predictive model.
# since sample size is very small, use entire dataset as training set.
# use the knn imputed results for model build

#---------------------------------------
# Resampleing method use Repeated 10-Fold CV 10 times.
# 10-fold is the default setting.
# Use ROC as the metric
# AUC is computed using ModelMetircs::auc, and the 
# ROC from result is the average of AUC value across 10 CV.
#(e.g. 10-Fold x 10 repeated =100 sets held-out sample)
# within each set of held-out sample, AUC is computed using ModelMetircs::auc, and the use modelfit$reample will show the 100 set of results for optimal parameter chosen. SD is just the sd(fit$resample$ROC)
# auc(ob_class, pred) #pred is the probability being in the long-live

ctrl <- trainControl(method = "repeatedcv",
                     repeats = 10, 
                     savePredictions = TRUE,
                     returnResamp="final",
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)
#-----------------------------------------
# Build the same models for control group
# A) Fit SVM, PLSDA, KNN, GLMNET, NSC in Ctr group using knn imputed data.
training.c = filter(reformat.model(knn.imp),treatment=="Ctr")[,-c(1:3)]
nzv.c <- nearZeroVar(training.c[,-1], saveMetrics= TRUE)
table(nzv.c$zeroVar) #none

#-----------------------------------
if(!file.exists("Ctr.fit.knn.Rdata")) {
  
  # SVM with Radial Basis Function Kernel
  set.seed(1237)
  svmFit.c <- train(trait ~ .,
                    data = training.c,
                    method = "svmRadial",
                    tuneLength = 10,
                    trControl = ctrl,
                    metric ="ROC",
                    preProcess = c("center","scale"))
  
  #svmFit.c
  #plot(svmFit.c, scales = list(x = list(log = 2)))
  #svmFit.c$finalModel
  
  ## norm = "overall", which is equivalento to "average" but in percentages.
  #confusionMatrix.train(svmFit.c)
  #dotPlot(varImp(svmFit.c, scale = FALSE))
  
  #-----------------------------
  # SVM with Linear Kernel
  set.seed(1237)
  svmFit.c2 <- train(trait ~ .,
                     data = training.c,
                     method = "svmLinear",
                     tuneLength = 10,
                     trControl = ctrl,
                     metric ="ROC",
                     preProcess = c("center","scale"))
  
  #-------------------------------------
  #PLSDA
  set.seed(1237)
  plsFit.c <- train(trait ~ .,
                    data = training.c,
                    method = "pls",
                    tuneLength = 15,
                    trControl = ctrl,
                    metric="ROC",
                    preProcess = c("center","scale"))
  
  #----------------------------
  # KNN
  set.seed(1237)
  knnFit.c <- train(trait ~ .,
                    data = training.c,
                    method = "knn",
                    tuneLength = 15,
                    trControl = ctrl,
                    metric="ROC",
                    preProcess = c("center","scale"))
  
  #--------------------------------
  # penalized logistic model
  glmGrid <- expand.grid(.alpha= c(0,0.1,0.2,0.4,0.6,0.8,1),
                         .lambda = seq(0.01, 0.2, length= 40))
  
  set.seed(1237)
  glmFit.c <- train(trait ~ .,
                    data = training.c,
                    method = "glmnet",
                    tuneGrid = glmGrid,
                    trControl = ctrl,
                    metric="ROC",
                    preProcess = c("center","scale"))
  
  #-----------------------
  ## Nearest Shrunken Centroids
  ## chose the specific range of tuning parameters
  nscGrid <- data.frame(.threshold = 0:25)
  set.seed(1237)
  nscFit.c <- train(trait ~ .,
                    data = training.c,
                    method = "pam",
                    tuneGrid = nscGrid,
                    trControl = ctrl,
                    metric="ROC",
                    preProcess = c("center","scale"))
  
  ## varImp will return the variable importance based on the
  ## distance between the class centroid and the overall centroid
  #dotPlot(varImp(nscFit.c))
  save(list = c("svmFit.c", "svmFit.c2","plsFit.c","knnFit.c",
                "glmFit.c", "nscFit.c"), file = "Ctr.fit.knn.Rdata")
}else{
  load("Ctr.fit.knn.Rdata")
}
#-------------------------
df.ctr.fit <- data.frame(Models=c("SVM (Radial Basis Function kernel)",
                                  "SVM (Linear kernel)",
                                  "PLSDA", "KNN", "PLM", "NSC"),
                         Tuning_paramters = c("cost, sigma","cost=1",
                                              "num comp","k",
                                              "alpha, lambda",
                                              "threshold"),
                         Optimal_values = c(paste0("cost=",svmFit.c$bestTune$C,
                                                   ", sigma=",round(svmFit.c$bestTune$sigma,4)),
                                            "cost=1", plsFit.c$bestTune$ncomp,
                                            knnFit.c$bestTune$k,
                                            paste0("alpha=",glmFit.c$bestTune$alpha,
                                                   ", lambda=",round(glmFit.c$bestTune$lambda,4)),
                                            nscFit.c$bestTune$threshold),
                         num_predictors = c(length(predictors(svmFit.c)),
                                            length(predictors(svmFit.c2)),
                                            length(predictors(plsFit.c)),
                                            length(predictors(knnFit.c)),
                                            length(predictors(glmFit.c)),
                                            length(predictors(nscFit.c))))

# comparing the models interms of their resampling results.
resamps.c <- resamples(list(SVM.r =svmFit.c, SVM.l =svmFit.c2,
                            PLSDA = plsFit.c,KNN = knnFit.c, 
                            PLM = glmFit.c,NSC = nscFit.c))

summary(resamps.c)
# Visualizing the resamples
splom(resamps.c, metric ="ROC", pscales = 0)
dotplot(resamps.c, metric = "ROC", main="Control group")

# comparing models
rocDiffs.c <- diff(resamps.c, metric = "ROC")
summary(rocDiffs.c)
dotplot(rocDiffs.c, metric= "ROC", main="Control group")

# clustering the models
plot(caret:::cluster.resamples(resamps.c), main = "Control group")

#-------------------------------------
# ROC curves
svm.c.roc = svmFit.c$pred[svmFit.c$pred$C ==svmFit.c$bestTune$C,]
svm.c2.roc = svmFit.c2$pred[svmFit.c2$pred$C ==svmFit.c2$bestTune$C,]
pls.c.roc = plsFit.c$pred[plsFit.c$pred$ncomp== plsFit.c$bestTune$ncomp,]
knn.c.roc = knnFit.c$pred[knnFit.c$pred$k== knnFit.c$bestTune$k,]
glm.c.roc = glmFit.c$pred[(glmFit.c$pred$alpha == glmFit.c$bestTune$alpha & 
                             glmFit.c$pred$lambda==glmFit.c$bestTune$lambda),]
nsc.c.roc = nscFit.c$pred[nscFit.c$pred$threshold==nscFit.c$bestTune$threshold,]
#Note that pROC::roc default using second level as interest (case);need to reverse the level.
#caret function treat the first level as 'positive' level.
ROC.C = list(svm.c.roc,svm.c2.roc,pls.c.roc,knn.c.roc,glm.c.roc,nsc.c.roc)
ROCCURVE.C = lapply(ROC.C, function(x) pROC::roc(x$obs, x$long_live,levels = rev(levels(x$obs))))

roccol= c("#00688B", "#EE3B3B", "#CDAD00", "#5F9EA0", "#00008B", "#EE7600")

pdf("plots/ROC_curves_control.knn.pdf")
plot(ROCCURVE.C[[1]], legacy.axes = TRUE, main="ROC curves- Control group",col=roccol[1])
for (i in 2:length(ROCCURVE.C)){
  plot(ROCCURVE.C[[i]], legacy.axes = TRUE, add=TRUE,col=roccol[i])
}
legend("bottomright", legend= paste0(df.ctr.fit$Models,": auc=",
                                     round(sapply(ROCCURVE.C, pROC::auc),4)),
       col=roccol,lty=1, cex=0.6)
dev.off()
#---------------------------------------
# B) Fit SVM, PLSDA, KNN, GLMNET, NSC in H2O2 group using knn imputed data.
training.h = filter(reformat.model(knn.imp),treatment=="H2O2")[,-c(1:3)]
nzv.h <- nearZeroVar(training.h[,-1], saveMetrics= TRUE)
table(nzv.h$zeroVar) #none
# Variable importance code can extract using:
#plsFit.h$modelInfo$varImp

#-----------------------------------
## SVM with Radial Basis Function Kernel
if(!file.exists("H2O2.fit.knn.Rdata")) {
  set.seed(1237)
  svmFit.h <- train(trait ~ .,
                    data = training.h,
                    method = "svmRadial",
                    tuneLength = 10,
                    trControl = ctrl,
                    metric ="ROC",
                    preProcess = c("center","scale"))
  
  #------------------------------
  ## SVM with Linear Kernel (caret fored C =1 ONLY)
  set.seed(1237)
  svmFit.h2 <- train(trait ~ .,
                     data = training.h,
                     method = "svmLinear",
                     tuneLength = 10,
                     trControl = ctrl,
                     metric ="ROC",
                     preProcess = c("center","scale"))
  
  #-----------------------------
  ##PLSDA
  set.seed(1237)
  plsFit.h <- train(trait ~ .,
                    data = training.h,
                    method = "pls",
                    tuneLength = 15,
                    trControl = ctrl,
                    metric="ROC",
                    preProcess = c("center","scale"))
  #----------------------------
  ## KNN
  set.seed(1237)
  knnFit.h <- train(trait ~ .,
                    data = training.h,
                    method = "knn",
                    tuneLength = 15,
                    trControl = ctrl,
                    metric="ROC",
                    preProcess = c("center","scale"))
  
  #--------------------------------
  ## penalized logistic model
  set.seed(1237)
  glmFit.h <- train(trait ~ .,
                    data = training.h,
                    method = "glmnet",
                    tuneGrid = glmGrid,
                    trControl = ctrl,
                    metric="ROC",
                    preProcess = c("center","scale"))
  
  #-----------------------
  ## Nearest Shrunken Centroids
  ## chose the specific range of tuning parameters
  set.seed(1237)
  nscFit.h <- train(trait ~ .,
                    data = training.h,
                    method = "pam",
                    tuneGrid = nscGrid,
                    trControl = ctrl,
                    metric="ROC",
                    preProcess = c("center","scale"))
  
  #----------------------------------
  save(list = c("svmFit.h", "svmFit.h2","plsFit.h","knnFit.h",
                "glmFit.h", "nscFit.h"), file = "H2O2.fit.knn.Rdata")
}else{
  load("H2O2.fit.knn.Rdata")
}

#-------------------------

df.h2o2.fit <- data.frame(Models=c("SVM (Radial Basis Function kernel)",
                                   "SVM (Linear kernel)",
                                   "PLSDA", "KNN", "PLM", "NSC"),
                          Tuning_paramters = c("cost, sigma","cost=1",
                                               "num comp","k",
                                               "alpha, lambda",
                                               "threshold"),
                          Optimal_values = c(paste0("cost=",svmFit.h$bestTune$C,
                                                    ", sigma=",round(svmFit.h$bestTune$sigma,4)),
                                             "cost=1", plsFit.h$bestTune$ncomp,
                                             knnFit.h$bestTune$k,
                                             paste0("alpha=",glmFit.h$bestTune$alpha,
                                                    ", lambda=",round(glmFit.h$bestTune$lambda,4)),
                                             nscFit.h$bestTune$threshold),
                          num_predictors = c(length(predictors(svmFit.h)),
                                             length(predictors(svmFit.h2)),
                                             length(predictors(plsFit.h)),
                                             length(predictors(knnFit.h)),
                                             length(predictors(glmFit.h)),
                                             length(predictors(nscFit.h))))

# comparing the models interms of their resampling results.
resamps.h <- resamples(list(SVM.r =svmFit.h, SVM.l =svmFit.h2,
                            PLSDA = plsFit.h,KNN = knnFit.h, 
                            PLM = glmFit.h, NSC = nscFit.h))

summary(resamps.h)
# Visualizing the resamples
splom(resamps.h, metric ="ROC", pscales = 0,pch=19)
dotplot(resamps.h, metric = "ROC", main="H2O2 group")

# comparing models
rocDiffs.h <- diff(resamps.h, metric = "ROC")
summary(rocDiffs.h)
dotplot(rocDiffs.h, metric= "ROC", main="H2O2 group")

# clustering the models
plot(caret:::cluster.resamples(resamps.h))
#-----------------------------------
# ROC curves
svm.h.roc = svmFit.h$pred[svmFit.h$pred$C ==svmFit.h$bestTune$C,]
svm.h2.roc = svmFit.h2$pred[svmFit.h2$pred$C ==svmFit.h2$bestTune$C,]
pls.h.roc = plsFit.h$pred[plsFit.h$pred$ncomp== plsFit.h$bestTune$ncomp,]
knn.h.roc = knnFit.h$pred[knnFit.h$pred$k== knnFit.h$bestTune$k,]
glm.h.roc = glmFit.h$pred[(glmFit.h$pred$alpha == glmFit.h$bestTune$alpha & 
                             glmFit.h$pred$lambda==glmFit.h$bestTune$lambda),]
nsc.h.roc = nscFit.h$pred[nscFit.h$pred$threshold==nscFit.h$bestTune$threshold,]

ROC.h = list(svm.h.roc, svm.h2.roc, pls.h.roc,
             knn.h.roc, glm.h.roc, nsc.h.roc)
ROCCURVE.h = lapply(ROC.h, function(x) pROC::roc(x$obs, x$long_live,levels = rev(levels(x$obs))))

pdf("plots/ROC_curves_H2O2.knn.pdf")
plot(ROCCURVE.h[[1]], legacy.axes = TRUE, 
     main="ROC curves- H2O2 group",col=roccol[1])
for (i in 2:length(ROCCURVE.h)){
  plot(ROCCURVE.h[[i]], legacy.axes = TRUE, add=TRUE,col=roccol[i])
}
legend("bottomright", legend= paste0(df.h2o2.fit$Models,": auc=",
                                     round(sapply(ROCCURVE.h, pROC::auc),4)),
       col=roccol,lty=1, cex=0.6)
dev.off()

#---------------------------------------------------

# use RF imputed data to build models.
#---------------------------------------------

#---------------------------------
# Build the same models for control group
# C) Fit SVM, PLSDA, KNN, GLMNET, NSC in Ctr group using RF imputed data.
training.c.rf = filter(reformat.model(rf.imp),treatment=="Ctr")[,-c(1:3)]
nzv.c.rf <- nearZeroVar(training.c.rf[,-1], saveMetrics= TRUE)
#-----------------------------------
if(!file.exists("Ctr.fit.rf.Rdata")) {
  
  # SVM with Radial Basis Function Kernel
  set.seed(1237)
  svmFit.c.rf <- train(trait ~ .,
                       data = training.c.rf,
                       method = "svmRadial",
                       tuneLength = 10,
                       trControl = ctrl,
                       metric ="ROC",
                       preProcess = c("center","scale"))
  
  #-----------------------------
  # SVM with Linear Kernel
  set.seed(1237)
  svmFit.c.rf2 <- train(trait ~ .,
                        data = training.c.rf,
                        method = "svmLinear",
                        tuneLength = 10,
                        trControl = ctrl,
                        metric ="ROC",
                        preProcess = c("center","scale"))
  
  #-------------------------------------
  #PLSDA
  set.seed(1237)
  plsFit.c.rf <- train(trait ~ .,
                       data = training.c.rf,
                       method = "pls",
                       tuneLength = 15,
                       trControl = ctrl,
                       metric="ROC",
                       preProcess = c("center","scale"))
  
  #----------------------------
  # KNN
  set.seed(1237)
  knnFit.c.rf <- train(trait ~ .,
                       data = training.c.rf,
                       method = "knn",
                       tuneLength = 15,
                       trControl = ctrl,
                       metric="ROC",
                       preProcess = c("center","scale"))
  
  #--------------------------------
  # penalized logistic model
  set.seed(1237)
  glmFit.c.rf <- train(trait ~ .,
                       data = training.c.rf,
                       method = "glmnet",
                       tuneGrid = glmGrid,
                       trControl = ctrl,
                       metric="ROC",
                       preProcess = c("center","scale"))
  
  #-----------------------
  ## Nearest Shrunken Centroids
  ## chose the specific range of tuning parameters
  
  set.seed(1237)
  nscFit.c.rf <- train(trait ~ .,
                       data = training.c.rf,
                       method = "pam",
                       tuneGrid = nscGrid,
                       trControl = ctrl,
                       metric="ROC",
                       preProcess = c("center","scale"))
  
  save(list = c("svmFit.c.rf", "svmFit.c.rf2","plsFit.c.rf","knnFit.c.rf",
                "glmFit.c.rf", "nscFit.c.rf"), file = "Ctr.fit.rf.Rdata")
}else{
  load("Ctr.fit.rf.Rdata")
}
#-------------------------
df.ctr.fit.rf <- data.frame(Models=c("SVM (Radial Basis Function kernel)",
                                     "SVM (Linear kernel)",
                                     "PLSDA", "KNN", "PLM", "NSC"),
                            Tuning_paramters = c("cost, sigma","cost=1",
                                                 "num comp","k",
                                                 "alpha,lambda",
                                                 "shrinkage threshold"),
                            Optimal_values = c(paste0("cost=",svmFit.c.rf$bestTune$C,
                                                      ", sigma=",round(svmFit.c.rf$bestTune$sigma,4)),
                                               "cost=1", plsFit.c.rf$bestTune$ncomp,
                                               knnFit.c.rf$bestTune$k,
                                               paste0("alpha=",glmFit.c.rf$bestTune$alpha,
                                                      ", lambda=",round(glmFit.c.rf$bestTune$lambda,4)),
                                               nscFit.c.rf$bestTune$threshold),
                            num_predictors = c(length(predictors(svmFit.c.rf)),
                                               length(predictors(svmFit.c.rf2)),
                                               length(predictors(plsFit.c.rf)),
                                               length(predictors(knnFit.c.rf)),
                                               length(predictors(glmFit.c.rf)),
                                               length(predictors(nscFit.c.rf))))

# comparing the models interms of their resampling results.
resamps.c.rf <- resamples(list(SVM.r =svmFit.c.rf, SVM.l =svmFit.c.rf2,
                               PLSDA = plsFit.c.rf,KNN = knnFit.c.rf, 
                               PLM = glmFit.c.rf,NSC = nscFit.c.rf))

# comparing models
rocDiffs.c.rf <- diff(resamps.c.rf, metric = "ROC")
summary(rocDiffs.c.rf)
dotplot(rocDiffs.c.rf, metric= "ROC", main="Control group")

#-------------------------------------
# ROC curves
svm.c.rf.roc = svmFit.c.rf$pred[svmFit.c.rf$pred$C ==svmFit.c.rf$bestTune$C,]
svm.c.rf2.roc = svmFit.c.rf2$pred[svmFit.c.rf2$pred$C ==svmFit.c.rf2$bestTune$C,]
pls.c.rf.roc = plsFit.c.rf$pred[plsFit.c.rf$pred$ncomp== plsFit.c.rf$bestTune$ncomp,]
knn.c.rf.roc = knnFit.c.rf$pred[knnFit.c.rf$pred$k== knnFit.c.rf$bestTune$k,]
glm.c.rf.roc = glmFit.c.rf$pred[(glmFit.c.rf$pred$alpha == glmFit.c.rf$bestTune$alpha & 
                                   glmFit.c.rf$pred$lambda==glmFit.c.rf$bestTune$lambda),]
nsc.c.rf.roc = nscFit.c.rf$pred[nscFit.c.rf$pred$threshold==nscFit.c.rf$bestTune$threshold,]

ROC.c.rf = list(svm.c.rf.roc,svm.c.rf2.roc,pls.c.rf.roc,
                knn.c.rf.roc,glm.c.rf.roc,nsc.c.rf.roc)
ROCCURVE.c.rf = lapply(ROC.c.rf, function(x) pROC::roc(x$obs, x$long_live,levels = rev(levels(x$obs))))

pdf("plots/ROC_curves_control.RF.pdf")
plot(ROCCURVE.c.rf[[1]], legacy.axes = TRUE, main="ROC curves- Control group (RF)",col=roccol[1])
for (i in 2:length(ROCCURVE.c.rf)){
  plot(ROCCURVE.c.rf[[i]], legacy.axes = TRUE, add=TRUE,col=roccol[i])
}
legend("bottomright", legend= paste0(df.ctr.fit.rf$Models,": auc=",
                                     round(sapply(ROCCURVE.c.rf, pROC::auc),4)),
       col=roccol,lty=1, cex=0.6)
dev.off()


#-------------------------------------
# D) Fit SVM, PLSDA, KNN, GLMNET, NSC in H2O2 group using RF imputed data.
training.h.rf = filter(reformat.model(rf.imp),treatment=="H2O2")[,-c(1:3)]
#-----------------------------------
## SVM with Radial Basis Function Kernel
if(!file.exists("H2O2.fit.rf.Rdata")) {
  
  set.seed(1237)
  svmFit.h.rf <- train(trait ~ .,
                       data = training.h.rf,
                       method = "svmRadial",
                       tuneLength = 10,
                       trControl = ctrl,
                       metric ="ROC",
                       preProcess = c("center","scale"))
  
  #------------------------------
  ## SVM with Linear Kernel (caret fored C =1 ONLY)
  set.seed(1237)
  svmFit.h.rf2 <- train(trait ~ .,
                        data = training.h.rf,
                        method = "svmLinear",
                        tuneLength = 10,
                        trControl = ctrl,
                        metric ="ROC",
                        preProcess = c("center","scale"))
  
  #-----------------------------
  ##PLSDA
  set.seed(1237)
  plsFit.h.rf <- train(trait ~ .,
                       data = training.h.rf,
                       method = "pls",
                       tuneLength = 15,
                       trControl = ctrl,
                       metric="ROC",
                       preProcess = c("center","scale"))
  
  #----------------------------
  ## KNN
  set.seed(1237)
  knnFit.h.rf <- train(trait ~ .,
                       data = training.h.rf,
                       method = "knn",
                       tuneLength = 15,
                       trControl = ctrl,
                       metric="ROC",
                       preProcess = c("center","scale"))
  
  #--------------------------------
  ## penalized logistic model
  set.seed(1237)
  glmFit.h.rf <- train(trait ~ .,
                       data = training.h.rf,
                       method = "glmnet",
                       tuneGrid = glmGrid,
                       trControl = ctrl,
                       metric="ROC",
                       preProcess = c("center","scale"))
  
  #-----------------------
  ## Nearest Shrunken Centroids
  ## chose the specific range of tuning parameters
  
  set.seed(1237)
  nscFit.h.rf <- train(trait ~ .,
                       data = training.h.rf,
                       method = "pam",
                       tuneGrid = nscGrid,
                       trControl = ctrl,
                       metric="ROC",
                       preProcess = c("center","scale"))
  
  
  save(list = c("svmFit.h.rf", "svmFit.h.rf2","plsFit.h.rf","knnFit.h.rf",
                "glmFit.h.rf", "nscFit.h.rf"), file = "H2O2.fit.rf.Rdata")
}else{
  load("H2O2.fit.rf.Rdata")
}

#-------------------------

df.h2o2.fit.rf <- data.frame(Models=c("SVM (Radial Basis Function kernel)",
                                      "SVM (Linear kernel)",
                                      "PLSDA", "KNN", "PLM", "NSC"),
                             Tuning_paramters = c("cost, sigma","cost=1",
                                                  "num comp","k",
                                                  "alpha, lambda",
                                                  "shrinkage threshold"),
                             Optimal_values = c(paste0("cost=",svmFit.h.rf$bestTune$C,
                                                       ", sigma=",round(svmFit.h.rf$bestTune$sigma,4)),
                                                "cost=1", plsFit.h.rf$bestTune$ncomp,
                                                knnFit.h.rf$bestTune$k,
                                                paste0("alpha=",glmFit.h.rf$bestTune$alpha,
                                                       " , lambda=",round(glmFit.h.rf$bestTune$lambda,4)),
                                                nscFit.h.rf$bestTune$threshold),
                             num_predictors = c(length(predictors(svmFit.h.rf)),
                                                length(predictors(svmFit.h.rf2)),
                                                length(predictors(plsFit.h.rf)),
                                                length(predictors(knnFit.h.rf)),
                                                length(predictors(glmFit.h.rf)),
                                                length(predictors(nscFit.h.rf))))

# comparing the models interms of their resampling results.
resamps.h.rf <- resamples(list(SVM.r =svmFit.h.rf, SVM.l =svmFit.h.rf2,
                               PLSDA = plsFit.h.rf,KNN = knnFit.h.rf, 
                               PLM = glmFit.h.rf, NSC = nscFit.h.rf))


# comparing models
rocDiffs.h.rf <- diff(resamps.h.rf, metric = "ROC")
summary(rocDiffs.h.rf)
dotplot(rocDiffs.h.rf, metric= "ROC" , main ="H2O2 group")

# clustering the models
plot(caret:::cluster.resamples(resamps.h.rf))

#-------------------------------------
# ROC curves
svm.h.rf.roc = svmFit.h.rf$pred[svmFit.h.rf$pred$C ==svmFit.h.rf$bestTune$C,]
svm.h.rf2.roc = svmFit.h.rf2$pred[svmFit.h.rf2$pred$C ==svmFit.h.rf2$bestTune$C,]
pls.h.rf.roc = plsFit.h.rf$pred[plsFit.h.rf$pred$ncomp== plsFit.h.rf$bestTune$ncomp,]
knn.h.rf.roc = knnFit.h.rf$pred[knnFit.h.rf$pred$k== knnFit.h.rf$bestTune$k,]
glm.h.rf.roc = glmFit.h.rf$pred[(glmFit.h.rf$pred$alpha == glmFit.h.rf$bestTune$alpha & 
                                   glmFit.h.rf$pred$lambda==glmFit.h.rf$bestTune$lambda),]
nsc.h.rf.roc = nscFit.h.rf$pred[nscFit.h.rf$pred$threshold==nscFit.h.rf$bestTune$threshold,]

ROC.h.rf = list(svm.h.rf.roc,svm.h.rf2.roc,pls.h.rf.roc,
                knn.h.rf.roc,glm.h.rf.roc,nsc.h.rf.roc)
ROCCURVE.h.rf = lapply(ROC.h.rf, function(x) pROC::roc(x$obs, x$long_live,levels = rev(levels(x$obs))))

pdf("plots/ROC_curves_H2O2.RF.pdf")
plot(ROCCURVE.h.rf[[1]], legacy.axes = TRUE, main="ROC curves- H2O2 group (RF)",col=roccol[1])
for (i in 2:length(ROCCURVE.h.rf)){
  plot(ROCCURVE.h.rf[[i]], legacy.axes = TRUE, add=TRUE,col=roccol[i])
}
legend("bottomright", legend= paste0(df.h2o2.fit.rf$Models,": auc=",
                                     round(sapply(ROCCURVE.h.rf, pROC::auc),4)),
       col=roccol,lty=1, cex=0.6)
dev.off()
#-------------------------------
# add new request: PCA plot by treatment and trait
pca.neg2 <- prcomp(t(na.omit(neg.matrix.filter)),scale. = TRUE)
pca.neg2 <- as.data.frame(pca.neg2$x[,1:3])

pca.neg2 <- merge(pca.neg2, sample.info[,c("neg.sample","treatment","trait","sampleInfo")], by.x="row.names", by.y="neg.sample")
pca.neg2$combine <- paste(pca.neg2$trait, pca.neg2$treatment,sep="_")
pca.neg2$combine <- as.factor(pca.neg2$combine)
mycol4 = c("#EEC900", "#CD6600", "#87CEFA", "#104E8B")

pca.neg.plot4 <- ggplot(pca.neg2, aes(x=PC1,y=PC2))+
  geom_point(aes(color=combine,
                 text = sampleInfo))+ 
  labs(title ="Neg PCA plot by trait&treatment filtered features post-normalization",
       x="PC 1", y="PC 2", colour = "")+
  scale_colour_manual(values=mycol4)


setwd("./plots")
saveWidget(ggplotly(pca.neg.plot4),"2D.pca.neg.by.traits.treatment.postfilter.html")

saveWidget(scatterplot3js(pca.neg2[,2:4], 
                          labels = paste(pca.neg2$combine),
                          color=mycol4[as.numeric(pca.neg2$combine)],
                          grid = FALSE, renderer ="canvas"),
           "pca.neg.by.traits.treatment.postfilter.html")

setwd("..")
#-------------------------------------------
# Use knn imputed data, build PLSDA models to predict  treatment group
# models build in long vs short group separately.
training.long = filter(reformat.model(knn.imp),trait=="long_live")[,-c(1,3,4)]

set.seed(1237)
plsFit.long <- train(relevel(treatment,ref="H2O2") ~ .,
                     data = training.long,
                     method = "pls",
                     tuneLength = 15,
                     trControl = ctrl,
                     metric="ROC",
                     preProcess = c("center","scale"))

plsFit.long.roc = plsFit.long$pred[plsFit.long$pred$ncomp== plsFit.long$bestTune$ncomp,]
confusionMatrix.train(plsFit.long)
#-----------------------
training.short = filter(reformat.model(knn.imp),trait=="short_live")[,-c(1,3,4)]
set.seed(1237)
plsFit.short <- train(relevel(treatment,ref="H2O2") ~ .,
                      data = training.short,
                      method = "pls",
                      tuneLength = 15,
                      trControl = ctrl,
                      metric="ROC",
                      preProcess = c("center","scale"))

confusionMatrix.train(plsFit.short)
#---------------------------
# Plot Roc curve
plsFit.short.roc = plsFit.short$pred[plsFit.short$pred$ncomp== plsFit.short$bestTune$ncomp,]


ROC.treatment = list(plsFit.long.roc, plsFit.short.roc)
ROCCURVE.treatment = lapply(ROC.treatment, function(x) pROC::roc(x$obs, x$Ctr,levels =rev(levels(plsFit.long.roc$obs))))

roccol2= c("#00688B", "#EE3B3B")

pdf("plots/ROC_curves_treatment.knn.pdf")
plot(ROCCURVE.treatment[[1]], legacy.axes = TRUE,
     main="ROC curves- PLSDA treatment group",col=roccol2[1])
plot(ROCCURVE.treatment[[2]], legacy.axes = TRUE, add=TRUE,col=roccol2[2])
legend("bottomright", legend= paste0(c("long_lived", "short_lived"),": Average auc=",
                                     c("0.962", "0.998")),
       col=roccol2,lty=1, cex=0.6)
dev.off()


#-------------------------------
save.image("Promislow_pid1652_negative_metabolite_analysis.RData")
#------------------------------------
```


```{r, echo = FALSE, results = "asis"}

BiocStyle::markdown()

```

**Author**: Lu Wang <br/>
  **Compiled**: `r date()`<br/>
  
  
  # Introduction
  
  The objective of this project is to compare missing value imputation strategies, and build univariate and multivariate predictive model for metabolomic data. In this dataset, there are `r ncol(neg.matrix)` samples, and `r nrow(neg.matrix)` metabolites. `r length(unique(sample.info$line))` genetic fly lines are treated with either control or $H_2O_2$. The lifespan (trait) which measured in the $H_2O_2$ group is dichotomized into two classes using 80 hr as the cutoff (46 samples in long live and 47 samples in short live group). Table below shows the sample information and stratification by treatment and trait.

```{r, echo = FALSE, results = "asis"}
DT::datatable(select(sample.info,neg.sample,line,treatment, lifespan,line.weight,sampleInfo,trait), 
              rownames= FALSE,
              caption = "Negative samples used in this study.")%>% formatRound(c(4,5),digit=2)

kable(table(sample.info$treatment,sample.info$trait),
      caption = "2 x 2 table of number of samples")

```



# Objectives and Results

## Objective 1. Exploratory data analysis (EDA).


Based on the PCA plots, we can see that samples are well separated by trait (long live vs. short live), but not by the treatment group (control vs. $H_2O_2$). 

<figure>
  ```{r, echo = FALSE,warning=FALSE,  results = "asis"}
ggplotly(pca.neg.plot2) 
ggplotly(pca.neg.plot1)
ggplotly(pca.neg.plot3)
ggplotly(pca.neg.plot.5)

plot_ly(data = pca.neg, x= ~PC1, y= ~PC2, z= ~PC3, 
        color = ~line, colors = rainbow(16), symbol =~treatment,
        symbols= 1:2, 
        marker = list(size=6))%>%
  add_markers() %>%
  layout(title = "3D PCA by line and treatment",
         scene = list(xaxis = list(title = "PC1"),
                      yaxis = list(title = "PC2"),
                      zaxis = list(title = "PC3")))

```
<figcaption>
  Fig1. 2D PCA plot using `r nrow(na.omit(neg.matrix))` features with no missing values. If you hover the cursor over the symbol, it will show the sample name. You can also click [3D PCA plot by trait](plots/pca.neg.by.traits.html) or [3D PCA plot by treatment](plots/pca.neg.by.treatment.html) or [3D PCA plot by both trait and treatment](plots/pca.neg.by.traits.treatment.html), these links will bring up a 3D plot that you can click on with your mouse, and rotate to inspect from different angles. 

The symbols are colored depending on trait (or treatment), and if you hover your mouse over a given symbol, the sample name will appear at the top left of the page.

NOTE: The four-way colored PCA plot showing above is generated with pre-filtering and sample normalized data. You can click [here](plots/2D.pca.neg.by.traits.treatment.postfilter.html) to view the one generated post-filtering and sample normalization.
</figcaption></figure>
  
  
  Missing values are known to be a problematic issue for metabolomic data. Currently, three types of missing value mechanisms have been identified [@Rubin:1976]. And a previous study has summarized them in terms of mass-spectrometry (MS)-based analysis [@Lazar:2016]: 
  **Missing Completely At Random (MCAR)**, which corresponds to the combination and propagation of multiple minor errors or stochastic fluctuations. Each missing value cannot be directly explained by the nature of the feature or by its measured intensity. MCAR is often combined with a more general class **Missing At Random (MAR)** when selecting imputation methods for MS analysis. **Missing Not At Random (MNAR)** which corresponds to chemical species whose abundances are close enough to the limit of detection of the instrument record a higher rate of missing values. 

Overall, there is `r round(mean(is.na(data.neg[,-1]))*100,2)`% of missingness found in the data, and `r round(mean(rowMeans(is.na(data.neg[,-1]))*100< 5)*100,2)`% of the metabolites have less than 5% missingness. Figure below shows the number of missing values stratified by treatments and traits.

<figure>
  ```{r, echo = FALSE, warning=FALSE, results = "asis"}
tb1.plot
```
<figcaption>
  Fig2. The number of missing values in data stratified by trait and treatment.
</figcaption></figure>
  
  
  In order to identify the missing value mechanisms in the data, we plot the percent of missingness versus the average  $log_2(mz)$ values for each metabolite across all the samples. As shown in the figure below, the metabolites with low average intensities are intended to have higher missingness, which suggested missing not at random (MNAR) mechanism maybe involved. However, for the metabolites with high intensities, MCAR/MAR maybe involved.


<figure>
  ```{r, echo = FALSE, results = "asis",fig.height=9}
plot(missing_by_feature$mean_log2mz, missing_by_feature$percent_missing, 
     pch=19, col="light blue", xlab=expression(paste("mean ", log[2],"(mz)")), ylab="Percent of missingness", main = "Negative samples")
points(missing_by_feature$mean_log2mz[missing_by_feature$percent_missing>5],
       missing_by_feature$percent_missing[missing_by_feature$percent_missing>5], pch=19, col="orange")
legend("topright", legend=c("<=5% missingness", "> 5% missingness"),
       pch=19, col=c("light blue", "orange"))

```
<figcaption>
  Fig3. Scatter plot of the percent of missingness versus the average $log_2(mz)$ values for each metabolite across all the samples. 
</figcaption></figure>
  
  
  
  ## Objective 2. Data preprocessing.
  
  ### Data normalization (by sample).
  To remove systematic variation between experimental conditions unrelated to the biological differences (i.e. dilutions, mass) we normalized the sample using median normalization method such that all the samples ended up with the same median value. Briefly, all the samples subtract their own median value (post $log_2$ transformation) then add the global median across all the data.  

NOTE: The original proposed centering and scaling by sample will force all the sample to follow the same distribution, which is based on a strong assumption. Further exploration maybe needed to test the impact of different normalization methods.

<figure>
  ```{r, echo = FALSE,warning=FALSE, fig.width=12, fig.height=8}
neg.box.bef
neg.box.aft

```
<figcaption> 
  Fig4. Boxplots of the $log_2$ transformed data before and after sample normalization. After normalization, the boxes line up on the horizontal line (median).
</figcaption></figure>
  
  
  
  ### Missing data imputation.
  
  We selected `r nrow(neg.matrix.filter)` features with less than 5% missingness for imputation. And percentage missingness after filtering is `r round(mean(is.na(neg.matrix.filter))*100,2)`%. Imputation methods are selected based on a previous study [@UHPLCMS:2016].

Four different imputation methods of missing values are applied: 
  
  a) Small value (SV) replacement: for every metabolite the missing values are replaced by the half of the minimum value of the entire dataset. 

b) Mean replacement: for every metabolite the missing values are replaced by the average value of the corresponding metabolite across all the samples.

c) K-nearest neighbor (KNN) imputation: for each metabolite with missing values, we find the k (here used 10) nearest neighbors using a Euclidean distance, confined to the columns (samples) for which that metabolite is not missing. For every metabolite the missing values are imputed by taking the average of those observed values of its neighbors. This is implemented in the Bioconductor impute package. The author of this package suggested KNN impute is robust towards the exact parameter used (number of k), and k=10 had shown good performance in the data introduced with 1-20% missingness [@impute:2001].


d) Random forest (RF) imputation: This is implemented in the R missForest package [@missForest:2012]. For each variable missForest fits a
random forest on the observed part and then predicts the missing part. The algorithm running iteratively, continuously updating the imputed matrix variable-wise, and is assessing its
performance between iterations. This assessment is done by considering the difference(s) between the previous imputation result and the new imputation result. As soon as this difference increases the algorithm stops.


In order to evaluate the difference performance of imputation methods, we perform PCA and compare the proportion of variance explained (PVE) by top 20 principal components (PCs).
```{r, echo = FALSE, results = "asis"}
kable(pve_20PCs, caption ="PVE by Top 20 PCs")
```

<figure>
  ```{r, echo = FALSE, warning= FALSE, fig.width=18,fig.height=15}
grid.arrange(sv.pca$pca.plot, 
             knn.pca$pca.plot, 
             mean.pca$pca.plot,
             rf.pca$pca.plot, ncol=2)

```
<figcaption> Fig5. 2D PCA plots for data after imputation. You can click [here](plots/PCA_imputations.pdf) to download the pdf.
</figcaption></figure>
  You can click here to view the 3D-PCA plots for [sv.imputation.PCA](plots/pca.neg.sv_imp.html), [mean.imputation.PCA](plots/pca.neg.mean_imp.html), [knn.imputation.PCA](plots/pca.neg.knn_imp.html), [rf.imputation.PCA](plots/pca.neg.rf_imp.html).

Based on the PCA results, the differences among different imputation methods are very small. 
A recent study [@Lazar:2016] suggested that choosing a method adapted to the nature of the missing value is more important than choosing a method itself. For example, KNN is designed for MCAR, while SV (or similar) method is the most naive method to deal with MNAR. The author also advised that in the absence of knowledge about the nature of the missing values, we should rely on a MCAR/MAR imputation method. Since we excluded the metabolites with > 5% missingness, the remaining ones are more likely to be related to the MCAR/MAR mechanism. In the multivariate section, we only select KNN and RF for model building.


## Objective 3. Univariate logistic models.


We fit logistic models for each metabolite with trait (long live vs short live) as the outcome variable, and the metabolite abundance ($log_2(mz)$) as the independent variable. This is done for control and $H_2O_2$ treated groups separately. We compare the coefficient change of the metabolite abundance with or without adjusting for the line weight, and  around `r round(num_confounder_Ctr[1],1)`% (control group) and `r round(num_confounder_H2O2[1],1)`% ($H_2O_2$ group) of metabolites have coefficient change > 10% after adjusting for the line weight, so we decided to keep it as a potential confounder for all the metabolites. 


Rather than modeling this trait directly, a logistic regression model allows us to establish a relationship between a binary outcome variable (e.g. trait) and a group of predictor variables (e.g. a metabolite's abundance, and line weight). It describes the relationship between a metabolite's abundance and the log odds of being in the long live group.
                                                                                                                                                                                               $$p(X) = Pr(trait=long live|X) = \frac{e^{\beta_0 + \beta_iX}}{1+e^{\beta_0 + \beta_iX}}$$ this can transform to:
                                                                                                                                                                                                 $$\frac{p(X)}{1-p(X)} = {e^{\beta_0 + \beta_iX}}$$
                                                                                                                                                                                                 $\frac{p(X)}{1-p(X)}$ called the odds, which is the probability of an event happening over an event not happening. If we take log of both sides, will become:
                                                                                                                                                                                                 $$log(\frac{p(X)}{1-p(X)}) = {\beta_0 + \beta_iX}$$
                                                                                                                                                                                                 The left-hand side is called the *log-odds* or *logit*, which is linear in *X*. This transformation is an attempt to get around the restricted range problem.  It maps probability ranging between 0 and 1 to log odds ranging from negative infinity to positive infinity.
                                                                                                                                                                                               
                                                                                                                                                                                               The interpretation of the logistic regression coefficients is somewhat tricky. 
                                                                                                                                                                                               For example, in control group fly, coefficient of M161T11 is 1.34 (no imputation data), let's fix the abundance of M161T11 at 14,  then the log-odds of being in long live trait when M161T11 =14 is:
                                                                                                                                                                                               $$log(\frac{p}{1-p}) |(M161T11= 14) =  \beta_0 + 1.34*14 + \beta_2*line.weight$$
                                                                                                                                                                                               If we increase one-unit of M161T11 to 15, then the log-odds of being in long live trait is:
                                                                                                                                                                                               $$log(\frac{p}{1-p}) |(M161T11= 15) =  \beta_0 + 1.34*15 + \beta_2*line.weight$$
                                                                                                                                                                                               Taking the difference of the two equations, we have:
                                                                                                                                                                                               $$log(\frac{p}{1-p}) |(M161T11= 15) - log(\frac{p}{1-p}) |(M161T11= 14) = 1.34$$
                                                                                                                                                                                               
                                                                                                                                                                                               which means for a one-unit increase in the abundance of metabolite M161T11, the expected change in log odds of being in long live fly controlling for the line weight is 1.34, and if we take exponentiation, we will have $\frac{odds(M161T11=15)}{odds(M161T11=14)}$ = exp(1.34) = `r round(exp(1.34),2)`. That is, for a one-unit increase in M161T11 abundance, we expect to see about 282% increase in the odds of being in a long live trait given the same value of line weight.
                                                                                                                                                                                               
                                                                                                                                                                                               Tables below list the number of metabolites significantly associated with log odds of being in the long live trait relative to the short live trait adjusting for line weight at a FDR cutoff of 0.1, which basically means that at most 10% of the metabolites we have selected could be false positives. We test with the completed data (no imputation) and 4 different imputation methods.
                                                                                                                                                                                               
                                                                                                                                                                                               ```{r, echo = FALSE, results = "asis"}
                                                                                                                                                                                               
                                                                                                                                                                                               sc <- function(x) x
                                                                                                                                                                                               
                                                                                                                                                                                               ctr.uni.df <- data.frame(Comparison =  sapply(ctab, function(x) XML::saveXML(Link(x))),
                                                                                                                                                                                               "Number significant" = sapply(ls.ctr, nrow), check.names = FALSE)
                                                                                                                                                                                               print(xtable(ctr.uni.df, caption = "Number of significant metabolites using univariate logistic models at 10% FDR in Control group"),
                                                                                                                                                                                               type = "html", include.rownames = FALSE, sanitize.text.function = sc, sanitize.colnames.function = sc)
                                                                                                                                                                                               
                                                                                                                                                                                               ```
                                                                                                                                                                                               
                                                                                                                                                                                               
                                                                                                                                                                                               ```{r, echo = FALSE, results = "asis"}
                                                                                                                                                                                               
                                                                                                                                                                                               sc2 <- function(x) x
                                                                                                                                                                                               
                                                                                                                                                                                               h2o2.uni.df <- data.frame(Comparison =  sapply(htab, function(x) XML::saveXML(Link(x))),
                                                                                                                                                                                               "Number significant" = sapply(ls.H2O2, nrow), check.names = FALSE)
                                                                                                                                                                                               print(xtable(h2o2.uni.df, caption = "Number of significant metabolites using univariate logistic models at 10% FDR in H2O2 group"),
                                                                                                                                                                                               type = "html", include.rownames = FALSE, sanitize.text.function = sc2, sanitize.colnames.function = sc2)
                                                                                                                                                                                               
                                                                                                                                                                                               ```
                                                                                                                                                                                               NOTE: The link in the above table will load an HTML table containing the
                                                                                                                                                                                               metabolites that are significantly different in this analysis. By default,
                                                                                                                                                                                               the table shows only 10 rows, and is sorted on the first column. You
                                                                                                                                                                                               can change the number of rows using the drop-down at the top left of
                                                                                                                                                                                               the page, and can sort on any column by simply clicking on the header
                                                                                                                                                                                               of the column you want to sort by. 
                                                                                                                                                                                               We also output the same table in Excel format - you can find it in the 
                                                                                                                                                                                               spreadsheets sub-directory that accompanies this document.
                                                                                                                                                                                               
                                                                                                                                                                                               
                                                                                                                                                                                               In addition, we generated four-way Venn diagrams to compare the results between completed data (no imputation) and different imputation methods. 
                                                                                                                                                                                               
                                                                                                                                                                                               <figure>
                                                                                                                                                                                               ```{r, echo = FALSE, results = "asis"}
                                                                                                                                                                                               grid.draw(venn.plot.c)
                                                                                                                                                                                               ```
                                                                                                                                                                                               
                                                                                                                                                                                               ```{r, echo = FALSE, results = "asis"}
                                                                                                                                                                                               grid.draw(venn.plot.h2)
                                                                                                                                                                                               ```
                                                                                                                                                                                               <figcaption>
                                                                                                                                                                                               Fig6. Results are presented in four-way Venn diagrams to compare different imputation methods and completed data. The results of KNN and RF in control group are identical, therefore, we only include KNN in the Venn diagram in control group. For $H_2O_2$ group, figure above excludes SV replacement result, you can click [here](plots/Venn_H2O2_sig_features_1.pdf) to view Venn diagram with SV but not mean replacement result. 
                                                                                                                                                                                               </figcaption></figure>
                                                                                                                                                                                               
                                                                                                                                                                                               In order to compare the significant features detected from two treatment groups, we generated 2-way Venn diagrams between two sets of results.
                                                                                                                                                                                               <figure>
                                                                                                                                                                                               ```{r, echo = FALSE, results = "asis"}
                                                                                                                                                                                               grid.draw(venn.plot.join4)
                                                                                                                                                                                               ```
                                                                                                                                                                                               <figcaption>
                                                                                                                                                                                               Fig7. The comparison of results of control group and $H_2O_2$ using KNN imputed data. 
                                                                                                                                                                                               You can also click [here](plots/Venn_uni_control_h2o2.pdf) to view Venn diagrams with other imputation methods, results are similar.
                                                                                                                                                                                               
                                                                                                                                                                                               
                                                                                                                                                                                               ## Objective 4. Multivariate predictive models.
                                                                                                                                                                                               
                                                                                                                                                                                               The goal is to use supervised learning methods to classify whether the fly belongs to long live or short live trait given their metabolite abundances (predictors). We selected six commonly used classification models: support vector machines (SVMs) with radial basis function (RBF) kernel, SVMs with linear kernel, partial least squares discriminant analysis (PLSDA), KNN, penalized logistic model (PLM), and nearest shrunken centroid model (NSC, a.k.a. PAM, for predictive analysis for microarrays). These predictors were modeled using either a normal distribution (e.g. PLSDA, NSC, PLM) or a non-parametric density (e.g. KNN, SVMs).
                                                                                                                                                                                               
                                                                                                                                                                                               The rational of our selection is based on a general strategy: Starting with models that are the least interpretable and most flexible, which have a high likelihood of producing the empirically optimum results (e.g. SVMs). Then investigate simpler models that are not complete black boxes, such as PLSDA [@Kuhn:2002]. Also, since this is a high-dimensional problem (p>>n), we considered models with build-in feature selection methods, such as PLM and NSC. 
                                                                                                                                                                                               
                                                                                                                                                                                               
                                                                                                                                                                                               **SVMs**: In classification, SVMs separate the different classes of data by a hyper-plane. In terms of classification performance, hyper-plane is the one with the maximal margin of separation between the classes, where "margin" is defined as the minimum distance from sample points to the hyperplane.  A subset of sample point(s) that lay on the margin is called support vectors. The cost parameter C of the SVM controls the penalty paid by the SVM for misclassifying a training point.  A high cost value C 
                                                                                                                                                                                               will force the SVM to create a prediction function that is complex enough to misclassify as few training points as possible. SVMs mapping the input data into a high-dimensional feature space by a kernel function (a function that quantifies the similarity of two observations). For example, the linear kernel is the inner products between the images of two data points in the feature space. In addition to linear kernel, we also tested a very popular RBF kernel, which can accommodate a non-linear boundary between class. There are two tuning parameters in the RBF kernel, cost and sigma (which is part of The RBF kernel function). Tuning parameter 'C' is held constant at a value of 1 in the linear kernel ([kernlab](https://cran.r-project.org/web/packages/kernlab/vignettes/kernlab.pdf) package). The predictor data should be centered and scaled prior to fitting so that attributes whose values are large in magnitude
                                                                                                                                                                                               do not dominate the calculations.
                                                                                                                                                                                               
                                                                                                                                                                                               
                                                                                                                                                                                               **PLSDA**: PLSDA which is an extension of the traditional PLS regression model, developed for classification setting ([pls](https://cran.r-project.org/web/packages/pls/vignettes/pls-manual.pdf) package). PLS finds latent variables that can maximize correlation with a continuous response value. Applying PLS in the classification setting with (a two-group problem), we would expect that the latent variables would be minimizing misclassification error. There is one tuning parameter: the number of latent variables to be retained.
                                                                                                                                                                                               
                                                                                                                                                                                               
                                                                                                                                                                                               **KNN**: KNN model for classification predicts a new sample using the K-closest samples from the training data, and the 'closeness' is determined via the distance metric (e.g. Euclidean distance) ([knn](https://cran.r-project.org/web/packages/kknn/kknn.pdf) package).
                                                                                                                                                                                               To allow each predictor to contribute equally to the distance calculation, we center and scale all the predictors prior to model fit.
                                                                                                                                                                                               Class probability estimates for the new sample are calculated as the proportion
                                                                                                                                                                                               of training set neighbors in each class. The new sample’s predicted class
                                                                                                                                                                                               is the class with the highest probability estimate. 
                                                                                                                                                                                               There is one tuning parameter: the number of neighbors (K) to be retained.
                                                                                                                                                                                               
                                                                                                                                                                                               
                                                                                                                                                                                               **Penalized logistic model**: PLM uses ridge and lasso penalties simultaneously on the binomial likelihood function ([glmnet](https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html) package). There are two tuning parameters: alpha and lambda.
                                                                                                                                                                                               Here, the alpha value is the “mixing proportion” that toggles between the pure lasso penalty (when alpha = 1) and a pure ridge-regression-like penalty (alpha = 0). The other tuning parameter lambda controls the total amount of penalization. In the situations where there are a large number of predictors and a small training set sample, the penalty term can stabilize the logistic regression
                                                                                                                                                                                               model coefficients. Like the lasso, this results in regression coefficients with values of absolute 0,
                                                                                                                                                                                               thus simultaneously accomplishing regularization and feature selection at the
                                                                                                                                                                                               same time. As with ridge regression, adding a penalty can
                                                                                                                                                                                               also provide a countermeasure against highly correlated predictors.
                                                                                                                                                                                               
                                                                                                                                                                                               
                                                                                                                                                                                               **Nearest Shrunken Centroids**: NSC was originally developed for DNA microarray data ([pamr](https://cran.r-project.org/web/packages/pamr/pamr.pdf) package), where the number of predictors is large and the number of samples is small. It is a linear classification model that is well suited for high-dimensional problems. For
                                                                                                                                                                                               each class, the centroid of the data is found by taking the average value of
                                                                                                                                                                                               each predictor (per class) in the training set. The overall centroid is computed
                                                                                                                                                                                               using the data from all of the classes.
                                                                                                                                                                                               If a predictor does not contain much information for a particular class, its
                                                                                                                                                                                               centroid for that class is likely to be close to the overall centroid.
                                                                                                                                                                                               
                                                                                                                                                                                               
                                                                                                                                                                                               
                                                                                                                                                                                               In general, when building a predictive model, we will split the data into a "training" set and a "test" set. The "training" data set is the general term for the samples used to create the model, while the "test" set is used to qualify performance. Ideally, the model should be evaluated on samples that are not used to build or fine-tune the model, so that they provide an unbiased sense of model effectiveness. However, our sample size is very small (46 or 47 samples per treatment group), we used resampling methods (a 10-fold cross-validation and repeated 10 times) to produce appropriate estimates of model performance from all the samples as the "training" set. For example, samples are randomly partitioned into 10 sets of roughly equal size, and the 9-sets (90%) of data were used to fit the model, and 1-set (10%) held-out data were used to validate the model performance [@Kuhn:2002].
                                                                                                                                                                                               
                                                                                                                                                                                               
                                                                                                                                                                                               Table below shows the confusion matrix example for generic classes “event” and “nonevent.” The top row of the table corresponds to samples predicted to be events. Some are predicted correctly (the true positives, or TP) while others are inaccurately classified (false positives or FP). Similarly, the second row contains the predicted negatives with
                                                                                                                                                                                               true negatives (TN) and false negatives (FN). Accuracy is defined as the fraction of correct classifications out of the total number of samples. 
                                                                                                                                                                                               Accuracy is calculated as: **(TP + TN)/Total**.
                                                                                                                                                                                               
                                                                                                                                                                                               ```{r, echo = FALSE, results = "asis"}
                                                                                                                                                                                               eg <- data.frame(" " = c(" ", "Predicted"," "),
                                                                                                                                                                                               " " = c(" ","event","no event"),
                                                                                                                                                                                               Reference = c("event","TP","FN"),
                                                                                                                                                                                               " " = c("non_event","FP","TN"))
                                                                                                                                                                                               names(eg) <- c("","","Reference", "")
                                                                                                                                                                                               
                                                                                                                                                                                               kable(eg)
                                                                                                                                                                                               ```
                                                                                                                                                                                               
                                                                                                                                                                                               
                                                                                                                                                                                               The sensitivity of the model is the rate that the event
                                                                                                                                                                                               of interest is predicted correctly for all samples having the event, or
                                                                                                                                                                                               **TP/(TP + FN)**. Conversely, the specificity is defined as
                                                                                                                                                                                               the rate that nonevent samples are predicted as nonevents, or **TN/(FP + TN)**.
                                                                                                                                                                                               
                                                                                                                                                                                               The Receiver operating characteristic (ROC) curve which combines these two parameters by plotting the sensitivity over the 1- specificity can also be used for a quantitative assessment of the model. A perfect model would have a ROC curve towards the left corner with the area under the ROC curve (AUC) of one. A completely ineffective model would result in an ROC curve that closely follows the 45 degree diagonal line and would have an AUC of approximately 0.5. 
                                                                                                                                                                                               
                                                                                                                                                                                               
                                                                                                                                                                                               Using the same cross-validation scheme, and the same resampling data sets by setting the same random number seed, we are able to compared methodologies based on resampling results. For example, a resampled estimate of the training set can also be obtained using confusion matrix and/or ROC metrics (e.g. sensitivity, specificity and AUC) are used as the measures of model performance, which are computed from the held-out samples and these values can be aggregated to diagnose issues with the model fit. And the optimal model parameter(s) is selected with the largest ROC (the average AUC value from resampling).
                                                                                                                                                                                               
                                                                                                                                                                                               
                                                                                                                                                                                               Data are centered and scaled by metabolite prior to model fit, such that all the metabolites have common standard deviation of one. All the steps are implemented in the R [caret](http://topepo.github.io/caret/index.html) package [@caret:2008].
                                                                                                                                                                                               
                                                                                                                                                                                               We list model parameters and results comparison among six models in control and $H_2O_2$ group separately. **Results below are generated using KNN imputed data.**
                                                                                                                                                                                               
                                                                                                                                                                                               
                                                                                                                                                                                               **Control group**
                                                                                                                                                                                               
                                                                                                                                                                                               Table below summarized the tuning parameter(s) and their optimal value(s) selected from the resampling in the control group. Since some models have built-in feature selection property, we also list number of predictors in the final model.
                                                                                                                                                                                               
                                                                                                                                                                                               
                                                                                                                                                                                               ```{r, echo = FALSE, results = "asis"}
                                                                                                                                                                                               knitr::kable(df.ctr.fit, caption = "Models used in control group data.")
                                                                                                                                                                                               ```
                                                                                                                                                                                               
                                                                                                                                                                                               
                                                                                                                                                                                               ```{r,echo=FALSE}
                                                                                                                                                                                               ctr.perf <- data.frame(Average.accuracy = c(sum(diag(confusionMatrix.train(svmFit.c)$table))/sum(confusionMatrix.train(svmFit.c)$table),                        sum(diag(confusionMatrix.train(svmFit.c2)$table))/sum(confusionMatrix.train(svmFit.c2)$table),                  sum(diag(confusionMatrix.train(plsFit.c)$table))/sum(confusionMatrix.train(plsFit.c)$table),    sum(diag(confusionMatrix.train(knnFit.c)$table))/sum(confusionMatrix.train(knnFit.c)$table),
                                                                                                                                                                                               sum(diag(confusionMatrix.train(glmFit.c)$table))/sum(confusionMatrix.train(glmFit.c)$table),                                            
                                                                                                                                                                                               sum(diag(confusionMatrix.train(nscFit.c)$table))/sum(confusionMatrix.train(nscFit.c)$table)),
                                                                                                                                                                                               Avarage_ROC =summary(resamps.c)$statistics[[1]][,4])
                                                                                                                                                                                               
                                                                                                                                                                                               ```
                                                                                                                                                                                               
                                                                                                                                                                                               
                                                                                                                                                                                               After comparing all six methods, we select PLSDA model which has the highest average ROC values to be the best predictive model. 
                                                                                                                                                                                               
                                                                                                                                                                                               <figure>
                                                                                                                                                                                               ```{r, echo = FALSE, results = "asis"}
                                                                                                                                                                                               plot(roc(pls.c.roc$obs,pls.c.roc$long_live,levels=rev(levels(pls.c.roc$obs))), legacy.axes = TRUE, 
                                                                                                                                                                                               main="PLSDA ROC curve- Control group",col=roccol[1])
                                                                                                                                                                                               text(x=0.2,y=0.3,col= roccol[6],
                                                                                                                                                                                               paste0("Average AUC=",
                                                                                                                                                                                               round(ctr.perf[row.names(ctr.perf)=="PLSDA",]$Avarage_ROC,3),
                                                                                                                                                                                               "\n","Average accuracy=",
                                                                                                                                                                                               round(ctr.perf[row.names(ctr.perf)=="PLSDA",]$Average.accuracy,3)))
                                                                                                                                                                                               
                                                                                                                                                                                               ```
                                                                                                                                                                                               <figcaption>Fig7. A plot of the ROC curve and AUC value of PLSDA model in control group.
                                                                                                                                                                                               </figcaption></figure>
                                                                                                                                                                                               
                                                                                                                                                                                               
                                                                                                                                                                                               
                                                                                                                                                                                               Last, we look at which feature(s) contribute most to the top components in the PLSDA model.
                                                                                                                                                                                               The variable importance measure here is based on weighted mean of the absolute regression coefficients. The weights are a function of the reduction of the sums of squares across the number of PLS components and are computed separately for each outcome. This is implemented with varImp function. Note that this is computed based on the 'training' data', results may change once we fit on 'test' data. You can click the bar to set up a filtering cutoff. The cutoff is arbitrary, just to give some ideas about which feature plays more important role. You can also click the buttons on the top-left corner to download or print the table.
                                                                                                                                                                                               
                                                                                                                                                                                               ```{r, echo = FALSE, results = "asis"}
                                                                                                                                                                                               varimp.Ctr = varImp(plsFit.c, scale =F)[[1]]
                                                                                                                                                                                               varimp.Ctr = varimp.Ctr[order(varimp.Ctr$Overall, decreasing = TRUE), , drop=FALSE]
                                                                                                                                                                                               DT::datatable((varimp.Ctr),filter = 'top',
                                                                                                                                                                                                             extensions = 'Buttons', 
                                                                                                                                                                                                             options = list(
                                                                                                                                                                                                               dom = 'Bfrtip',
                                                                                                                                                                                                               buttons = c('copy', 'csv', 'pdf', 'print'))) %>%formatRound(1,5)
                                                                                                                                                                                               
                                                                                                                                                                                               
                                                                                                                                                                                               ```
                                                                                                                                                                                               
                                                                                                                                                                                               <figure>
                                                                                                                                                                                                 ```{r, echo = FALSE, results = "asis"}
                                                                                                                                                                                               dotPlot(varImp(plsFit.c, scale = F))
                                                                                                                                                                                               
                                                                                                                                                                                               ```
                                                                                                                                                                                               <figcaption> Fig8. Needle plot of the importance values where the predictors are sorted from most important to least. Only top 20 features are show here.
                                                                                                                                                                                               </figcaption></figure>
                                                                                                                                                                                                 
                                                                                                                                                                                                 **$H_2O_2$ group**
                                                                                                                                                                                                 
                                                                                                                                                                                                 Table below summarized the tuning parameter(s) and their optimal value(s) selected from the resampling in $H_2O_2$ group.  
                                                                                                                                                                                               ```{r, echo = FALSE, results = "asis"}
                                                                                                                                                                                               knitr::kable(df.h2o2.fit, caption = "Models used in H2O2 group data.")
                                                                                                                                                                                               ```
                                                                                                                                                                                               
                                                                                                                                                                                               Table below shows the average percentage accuracy and average AUC value. These values are estimated from the **optimal** tuning parameter(s) determined from resampling.
                                                                                                                                                                                               ```{r,echo=FALSE}
                                                                                                                                                                                               h2o2.perf <- data.frame(Average.accuracy = c(sum(diag(confusionMatrix.train(svmFit.h)$table))/sum(confusionMatrix.train(svmFit.h)$table),                        sum(diag(confusionMatrix.train(svmFit.h2)$table))/sum(confusionMatrix.train(svmFit.h2)$table),                  sum(diag(confusionMatrix.train(plsFit.h)$table))/sum(confusionMatrix.train(plsFit.h)$table),    sum(diag(confusionMatrix.train(knnFit.h)$table))/sum(confusionMatrix.train(knnFit.h)$table),
                                                                                                                                                                                                                                            sum(diag(confusionMatrix.train(glmFit.h)$table))/sum(confusionMatrix.train(glmFit.h)$table),                                            
                                                                                                                                                                                                                                            sum(diag(confusionMatrix.train(nscFit.h)$table))/sum(confusionMatrix.train(nscFit.h)$table)),
                                                                                                                                                                                                                       Avarage_ROC =summary(resamps.h)$statistics[[1]][,4])
                                                                                                                                                                                               
                                                                                                                                                                                               ```
                                                                                                                                                                                               
                                                                                                                                                                                               
                                                                                                                                                                                               After comparing all six methods, we select PLSDA model which has the highest average ROC values to be the best predictive model. 
                                                                                                                                                                                               
                                                                                                                                                                                               <figure>
                                                                                                                                                                                                 ```{r, echo = FALSE, results = "asis"}
                                                                                                                                                                                               
                                                                                                                                                                                               plot(pROC::roc(pls.h.roc$obs,pls.h.roc$long_live,levels=rev(levels(pls.h.roc$obs))), legacy.axes = TRUE, 
                                                                                                                                                                                                    main="PLSDA ROC curve- H2O2 group",col=roccol[1])
                                                                                                                                                                                               text(x=0.2,y=0.3,col= roccol[6],
                                                                                                                                                                                                    paste0("Average AUC=",
                                                                                                                                                                                                           round(h2o2.perf[row.names(h2o2.perf)=="PLSDA",]$Avarage_ROC,3),
                                                                                                                                                                                                           "\n","Average accuracy=",
                                                                                                                                                                                                           round(h2o2.perf[row.names(h2o2.perf)=="PLSDA",]$Average.accuracy,3)))
                                                                                                                                                                                               
                                                                                                                                                                                               ```
                                                                                                                                                                                               <figcaption>Fig9. A plot of the ROC curve and AUC value of PLSDA model in $H_2O_2$ group.
                                                                                                                                                                                               </figcaption></figure>
                                                                                                                                                                                                 
                                                                                                                                                                                                 
                                                                                                                                                                                                 Last, we look at which feature(s) contribute most to the top components in the PLSDA model.
                                                                                                                                                                                               The variable importance measure here is based on weighted sums of the absolute regression coefficients.
                                                                                                                                                                                               
                                                                                                                                                                                               ```{r, echo = FALSE, results = "asis"}
                                                                                                                                                                                               varimp.H2O2 = varImp(plsFit.h, scale =F)[[1]]
                                                                                                                                                                                               DT::datatable((varimp.H2O2),filter = 'top',
                                                                                                                                                                                                             extensions = 'Buttons', 
                                                                                                                                                                                                             options = list(
                                                                                                                                                                                                               dom = 'Bfrtip',
                                                                                                                                                                                                               buttons = c('copy', 'csv', 'pdf', 'print'))) %>%formatRound(1,5)
                                                                                                                                                                                               
                                                                                                                                                                                               
                                                                                                                                                                                               ```
                                                                                                                                                                                               
                                                                                                                                                                                               <figure>
                                                                                                                                                                                                 ```{r, echo = FALSE, results = "asis"}
                                                                                                                                                                                               dotPlot(varImp(plsFit.h, scale = F))
                                                                                                                                                                                               
                                                                                                                                                                                               ```
                                                                                                                                                                                               <figcaption> Fig10. Needle plot of the importance values where the predictors are sorted from most important to least. Only top 20 features are show here.
                                                                                                                                                                                               </figcaption></figure>
                                                                                                                                                                                                 
                                                                                                                                                                                                 ### Summary
                                                                                                                                                                                                 
                                                                                                                                                                                                 We also test the same set of models using RF imputed data, though the optimal tuning parameter not exactly the same, the performance of the models are very similar, you can click [RF.results](RF.fit.neg.html) to view the detailed summary. Briefly, in control group, PLSDA, PLM and SVM (linear kernel) perform the best, and there is no significant difference in their ROC values from resampling results. Similar results are found in the $H_2O_2$ group. 
                                                                                                                                                                                               
                                                                                                                                                                                               Interestingly, SVM with linear kernel performs better than a more flexible RBF kernel. In this low n (sample), high p (feature) scenario,
                                                                                                                                                                                               the data probably cannot support a highly nonlinear model, and linear classification
                                                                                                                                                                                               boundaries probably are better choices.
                                                                                                                                                                                               
                                                                                                                                                                                               
                                                                                                                                                                                               
                                                                                                                                                                                               ### Limitations
                                                                                                                                                                                               
                                                                                                                                                                                               We are aware of the layers of uncertainties introduced from using different normalization methods, imputation methods, and different performance estimates (e.g. accuracy, kappa vs ROC). More important, in the current predictive modeling build, we used the entire data as our training data, we would need independent data sets to validate our results.
                                                                                                                                                                                               
                                                                                                                                                                                               ## Objective 5. Univariate linear regression models.
                                                                                                                                                                                               
                                                                                                                                                                                               
                                                                                                                                                                                               In order to determine if the differences between $H_2O_2$ treated and control group are different between the long lived and short lived flies, we built a linear regression model with both main effect and the interaction term between trait and treatment adjusting for the line weight. 
                                                                                                                                                                                               
                                                                                                                                                                                               
                                                                                                                                                                                               ```{r, echo = FALSE, results = "asis"}
                                                                                                                                                                                               
                                                                                                                                                                                               sc3 <- function(x) x
                                                                                                                                                                                               
                                                                                                                                                                                               uni.int <- data.frame(Comparison =  sapply(inter.tab, function(x) XML::saveXML(Link(x))),
                                                                                                                                                                                                                     "Number significant" = sapply(ls.int, nrow), check.names = FALSE)
                                                                                                                                                                                               print(xtable(uni.int, caption = "Number of metabolites with significant interactions between trait and treatment using univariate linear regression models at 10% FDR"),
                                                                                                                                                                                                     type = "html", include.rownames = FALSE, sanitize.text.function = sc3, sanitize.colnames.function = sc3)
                                                                                                                                                                                               
                                                                                                                                                                                               ```
                                                                                                                                                                                               NOTE: We also output the same table in Excel format - you can find it in the 
                                                                                                                                                                                               spreadsheets sub-directory that accompanies this document.
                                                                                                                                                                                               
                                                                                                                                                                                               We also plot the $log_2$ abundance of these metabolites stratified by trait and treatment group (knn imputed data), which will help you visualize the expression level of each feature. You can click [here](plots/scatter.plot.inter.pdf) to download the pdf version.
                                                                                                                                                                                               
                                                                                                                                                                                               In addition, we generated a four-way Venn diagram to compare the results between completed data (no imputation) and different imputation methods. 
                                                                                                                                                                                               
                                                                                                                                                                                               <figure>
                                                                                                                                                                                                 ```{r, echo = FALSE, results = "asis"}
                                                                                                                                                                                               grid.draw(venn.plot.int)
                                                                                                                                                                                               ```
                                                                                                                                                                                               <figcaption>
                                                                                                                                                                                                 Fig11. Results are presented in a four-way Venn diagram to compare different imputation methods and completed data. The results are very similar.
                                                                                                                                                                                               
                                                                                                                                                                                               In order to determine if the effect of $H_2O_2$ differs between long lived and short lived flies, we generated a heatmap of the `r nrow(sig_knn_imp_uni_int)` features that showed a significant interaction between trait and treatment adjusting for line weight (KNN imputed data). We use the residuals after removing the effect of line weight to compute a Z-score by centering and scaling across 4 groups for each metabolites. Features are clustered using hierarchical clustering method via the hclust function with the "complete agglomeration" method. Distance metrics for clustering are computed using Euclidean distance. You can also click [here](plots/heatmap_sig_knn_inter_scaled.html) to view an interactive version.
                                                                                                                                                                                               <figure>
                                                                                                                                                                                                 ```{r, echo = FALSE, results = "asis", fig.height=16,fig.width=8}
                                                                                                                                                                                               heatmap_sig_int
                                                                                                                                                                                               ```
                                                                                                                                                                                               </figure>
                                                                                                                                                                                                 
                                                                                                                                                                                                 
                                                                                                                                                                                                 
                                                                                                                                                                                                 ## Objective 6. Use PLSDA model to classify treatment groups.
                                                                                                                                                                                                 
                                                                                                                                                                                                 We also built two PLSDA models using metabolite abundance to classify treatment in the long lived and short lived groups, separately. Shown below are ROC curves from both PLSDA models.
                                                                                                                                                                                               
                                                                                                                                                                                               ```{r, echo = FALSE, results = "asis"}
                                                                                                                                                                                               plot(ROCCURVE.treatment[[1]], legacy.axes = TRUE,
                                                                                                                                                                                                    main="ROC curves- PLSDA treatment group",col=roccol2[1])
                                                                                                                                                                                               plot(ROCCURVE.treatment[[2]], legacy.axes = TRUE, add=TRUE,col=roccol2[2])
                                                                                                                                                                                               legend("bottomright", legend= paste0(c("long_lived", "short_lived"),": Average auc=",
                                                                                                                                                                                                                                    c("0.962", "0.998")),
                                                                                                                                                                                                      col=roccol2,lty=1, cex=0.6)
                                                                                                                                                                                               
                                                                                                                                                                                               ```
                                                                                                                                                                                               <figcaption>Fig12. A plot of the ROC curves and AUC value of PLSDA model in long lived and short lived group.
                                                                                                                                                                                               </figcaption></figure>
                                                                                                                                                                                                 
                                                                                                                                                                                                 
                                                                                                                                                                                                 ## Objective 7. Pathway analysis with Mummichog.
                                                                                                                                                                                                 
                                                                                                                                                                                                 We used  *Mummichog* (version 1.05) [@Mummichog:2013], to perform putative metabolite annotation and pathway analysis.
                                                                                                                                                                                               The metabolic model Drosophila melanogaster (BioCyc 16.5) was used to test the enrichment of pathways and networks. All the metabolites that fitted in the univariate logistic models, and univariate linear regression models were used as the input file (m/z values, retention time, FDR, and Z-statistics), FDR cutoff was set at 0.1 and 100 permutations were ran to estimate the null distribution. Mode (-m) select negative. For the univariate logistic model results, analysis were done for control and $H_2O_2$ group, separately. More detailed results are included in the Mummichog_Results subfolder.
                                                                                                                                                                                               
                                                                                                                                                                                               ## Objective 8. Visualizations in MetScape (CytoScape3 APP).
                                                                                                                                                                                               
                                                                                                                                                                                               Mummichog outputted .txt file under the 'sif/' directory are files intended for Cytoscape (cytoscape.org) visualization (File -> Import -> Network -> File). *.cys* files are saved under Cytoscape3 subfolders under each of the mummichog results folder.
                                                                                                                                                                                               
                                                                                                                                                                                               We also used *MetScape 3* to build pathway-based networks (compound, compound-gene) using Mummichog outputs (modules, ActivityNetworks). MetScape is a tool to visualize and interpret metabolomics and gene expression data in the context of human metabolic networks. It uses KEGG compound database and Edinburgh Human Metabolic Network as information sources. Network Type selected *Compound* or *Compund-Gene*, and organism selected *Human*.
                                                                                                                                                                                               
                                                                                                                                                                                               
                                                                                                                                                                                               I first convert the compound ID to compound names based ont the feature information generated from Mummichog (tsv/_tentative_featurematch_.xlsx). Please note that MetScape only recognize compound names or KEGG IDs, and it doesn't recognize the m/z values or the metabolite names from raw data file (e.g.  M184T5). And MetScape doesn't recognize all the input list of compounds, and I tried to convert those 'missing' compound names to KEGG IDs, and MetScape still failed to recognize those compounds. 
                                                                                                                                                                                               
                                                                                                                                                                                               You will find 6 files within each of the MetScape sub-folders, the list of compound names that used as input are included in *SOURCE.TARGET.xlsx. 
                                                                                                                                                                                               The unrecognized compounds are saved as *MetScape_missing.txt. Compound-gene relationship and gene discriptions are saved in _outputfile.xlsx.
                                                                                                                                                                                               Each of the secssion is saved as .cys file, and if you open it with CytoScape, you can view the generated *Compound* or *Compund-Gene* Networks, which are also saved as in pdf files.
                                                                                                                                                                                               
                                                                                                                                                                                               # R packages used
                                                                                                                                                                                               
                                                                                                                                                                                               
                                                                                                                                                                                               ```{r, echo = FALSE}
                                                                                                                                                                                               
                                                                                                                                                                                               sessionInfo()
                                                                                                                                                                                               
                                                                                                                                                                                               ```
                                                                                                                                                                                               
                                                                                                                                                                                               # References
                                                                                                                                                                                               