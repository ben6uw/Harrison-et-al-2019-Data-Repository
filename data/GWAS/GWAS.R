---
# Peroxide GWAS
---
library(tidyr)
library(splitstackshape) 
library(qqman) # GWAS analysis
library(devtools)
library(qvalue) # methods for genome-wide FDR
library(AssocTests) # calculates Tracy-Windom statistic for genotype PCs
library(corrplot)
library(AnnotationDbi) # pull annotations for genes
library(org.Dm.eg.db) # Drosophila-specific gene name conversion

mycol <- c("#8B8378", "#EE9A00") 
mycol2 <- c(mycol,"#483D8B")
mycol4 <- c("#EEC900", "#CD6600", "#87CEFA", "#104E8B")


## load PLink mapping data [see: ../GWAS_plink.log for details of PLINK call]
read.table("../linear-4-hiddencovars.maf0.05.assoc.linear", stringsAsFactors = F, header=T) -> plink
plink[plink$TEST == 'ADD', ] -> plink # removes results from the covariates (e.g. COV1, COV2...) for each SNP. What's left is ADD, which is the asociation test for the SNP + covariates on the phentoype in the additive model (in the DGRP this is 'two allele doses' vs. 'zero allele doses')  
plink[!is.na(plink$BP), ] -> plink 
plink[!is.na(plink$P), ] -> plink 
head(plink)


################################################################################################
## test the significance of each genotype PC using the Tracy-Windom statistic:
################################################################################################
read.table("../perox.eigenval", stringsAsFactors = F, header=F) -> eigenvalues # eigenvalues from PLINK'c --pca call to the bed file containing only lines measured in the peroxide trials at <30% missingess and >5% MAF
head(eigenvalues)
plot(eigenvalues$V1/sum(eigenvalues$V1), pch=19, las=1)

tw(eigenvalues$V1, length(eigenvalues$V1), criticalpoint = 0.9793) -> result # this 'critical point' corresponds to an alpha of 0.05.  I ran this liberal alpha to bias toward significance.
plot(result$statistic, pch=19, las=1, ylab='Tracy-Windom statistic', xlab='Genotype Principal Component')
abline(h=0, lty=2)
result$SigntEigenL # number of significant eigenvalues
################################################################################################


################################################################################################
## Q-Q plot
qqman::qq(plink$P, main=NULL, las=1, cex=0.8, cex.axis = 1)

quartz.save("../peroxide.QQ.plot.jpg", type = "jpg", device = dev.cur(), dpi = 300)
################################################################################################

p.adjust(plink$P) -> plink$FDR
table(plink$FDR)

significant.SNPS <- c('2L_530188_SNP', 
                      '2L_530218_SNP', 
                      '2L_531602_SNP',
                      '2L_536177_SNP',
                      '2L_536735_SNP',
                      '2L_547678_SNP',
                      '2L_6354477_SNP',
                      '2L_6869359_SNP',
                      '2L_6870147_SNP',
                      '2L_18439858_SNP',
                      '2R_12129030_SNP',
                      '3L_3143787_SNP',
                      '3L_10441133_SNP',
                      '3L_10441227_SNP') 

################################################################################################
# LD analysis among the SNPS:
################################################################################################

read.csv("../dgrp2.simplified.csv", stringsAsFactors = F, header=T) -> calls # reads a numeric df of the DGRP genotypes
row.names(calls) <- calls[ ,1]
calls[ ,-c(1)] -> calls

calls[rownames(calls) %in% significant.SNPS, ] -> genotypes

# get the identify of the DGRP lines measured in the study:
read.table("../peroxide.residuals.txt", header=T, stringsAsFactors = F) -> phenotype
cSplit(phenotype, 'FID', sep='_', drop=T) -> phenotype
paste0('Ral_', phenotype$FID_2) -> lines.measured
genotypes[ ,lines.measured] -> genotypes

# rm(calls)
# rm(phenotype)

t(genotypes) -> tgenotypes
cor(tgenotypes, use = 'complete.obs') -> c
abs(c) -> absc # get the absolute value of r 
corrplot.mixed(absc, upper='color', diag=NULL, lower.col=0, tl.pos = 'lt', tl.col=1)
corrplot.mixed(absc, lower='color', diag=NULL, upper.col=0, tl.pos = 'lt', tl.col=1)
corrplot(absc, method='color', diag=F, type='upper', tl.pos = 'lt', tl.col=1)
corrplot(absc, method='color', diag=T, tl.pos = 'lt', tl.col=1, tl.cex=1)
################################################################################################



################################################################################################
## Make a Manattan Plot:
################################################################################################
# reorder the chromosomes:
table(plink$CHR )   
plink$CHR + 1 -> plink$CHR
table(plink$CHR)
plink$CHR[plink$CHR == 6] <- 1 
plink$CHR[plink$CHR == 7] <- 6 
table(plink$CHR)
#

manhattan(plink, chr = "CHR", bp = "BP", p = "P", snp = as.character("SNP"), chrlabs = c('X', '2L', '2R','3L','3R','4'), genomewideline = F, highlight = NULL, logp = T, cex=0.5, cex.main = 1, cex.axis = 1,  cex.lab=1, ylim=c(2, 8), suggestiveline = F)
#
