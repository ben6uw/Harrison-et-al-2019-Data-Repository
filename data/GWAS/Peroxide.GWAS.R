---
title: "GWAS"
author: "Ben"
date: "1/18/2017" #with many revisions
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


## load PLink mapping data:
read.table("/Users/ben/Google_Drive/Applications/plink_mac/DGRP/four.pcs.assoc.linear", stringsAsFactors = F, header=T) -> plink
plink[plink$TEST == 'ADD', ] -> plink # removes results from the covariates (e.g. COV1, COV2...) for each SNP.  I can't remember what those are.  What's left is ADD, which is the asociation test for the SNP + covariates on the phentoype in the additive model (in the DGRP this is 'two allele doses' vs. 'zero allele doses')  
plink[!is.na(plink$BP), ] -> plink # removes 'unmapped?' SNPs
plink[!is.na(plink$P), ] -> plink 
head(plink)


################################################################################################
## test the significance of each genotype PC using the Tracy-Windom statistic:
################################################################################################
read.table("/Users/ben/Google_Drive/Applications/plink_mac/DGRP/perox.eigenval", stringsAsFactors = F, header=F) -> eigenvalues # eigenvalues from PLINK'c --pca call to the bed file containing only lines measured in the peroxide trials at <30% missingess and >5% MAF
head(eigenvalues)
plot(eigenvalues$V1/sum(eigenvalues$V1), pch=19, las=1)

tw(eigenvalues$V1, length(eigenvalues$V1), criticalpoint = 0.9793) -> result # this 'critical point' corresponds to an alpha of 0.05.  I ran this liberal alpha to bias toward significance.
plot(result$statistic, pch=19, las=1, ylab='Tracy-Windom statistic', xlab='Genotype Principal Component')
abline(h=0, lty=2)
result$SigntEigenL # number of significant eigenvalues
################################################################################################


################################################################################################
## Q-Q plot
################################################################################################
qqman::qq(plink$P,main=NULL,las=1)

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
                      '3L_10441227_SNP') # this list from the original GWAS


################################################################################################
# LD analysis among the SNPS:
################################################################################################

read.csv("/Users/ben/Google_Drive/Documents/generic.GWAS.analysis/dgrp2.simplified.csv", stringsAsFactors = F, header=T) -> calls # reads a numeric df of the DGRP genotypes
row.names(calls) <- calls[ ,1]
calls[ ,-c(1)] -> calls

head(calls)
calls[rownames(calls) %in% significant.SNPS, ] -> genotypes

significant.SNPS

# get the identify of the DGRP lines measured in the study:
read.table("/Users/ben/Google_Drive/Applications/plink_mac/DGRP/peroxide.residuals.txt", header=T, stringsAsFactors = F) -> phenotype
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
### first need to reorder the chromosomes:
table(plink$CHR )   
plink$CHR + 1 -> plink$CHR
table(plink$CHR)
plink$CHR[plink$CHR == 6] <- 1 
plink$CHR[plink$CHR == 7] <- 6 
table(plink$CHR)
################################################

par(mfrow=c(1,1))
quartz(height = 3, width = 5)
manhattan(plink, chr = "CHR", bp = "BP", p = "P", snp = as.character("SNP"), chrlabs = c('X', '2L', '2R','3L','3R','4'), genomewideline = F, highlight = NULL, logp = T, cex=0.5, cex.main = 1, cex.axis = 1,  cex.lab=1, ylim=c(2, 8), suggestiveline = F)

quartz.save("/Users/ben/Google_Drive/Documents/Peroxide/Results from online GWAS/peroxide.plink.manhattan.jpg", type = "jpg", device = dev.cur(), dpi = 400)

dev.off()
dev.off()
################################################################################################



################################################################################################
## read in and merge annotations if desired:
################################################################################################
read.table('/Users/ben/Google_Drive/Documents/generic.GWAS.analysis/DGRP.SNP_gene_annotations_all_genes_per_variant', header=T, stringsAsFactors = F, sep=',') -> annot
head(annot)

annot[annot$ID %in% plink$SNP, ] -> annot

annot[ ,c(1, 6:17)] -> x
head(x)
gather(x, key=gene.number, val=gene, -ID, factor_key = F) -> annot.long
head(annot.long)
annot.long[!is.na(annot.long$gene), ] -> annot.long

merge(annot.long, plink, by.x='ID', by.y='SNP') -> m
head(m)

aggregate(m, by=list(m$gene), function(x) min(x)) -> Pmins # formula version of aggregate - nice.
head(Pmins)

read.table("/Users/ben/Google_Drive/Documents/Peroxide/Peroxide Paper/Table1.txt", sep='\t', header=T, stringsAsFactors = F) -> table.1
head(table.1) 
merge(table.1, Pmins, by.x='FlyBase.ID', by.y='Group.1') -> temp
head(temp)
temp[ ,c(1:3, 6, 16, 4:5)] -> temp
names(temp) <- c("FlyBase.ID", "Gene", "Name", "Lead Variant", "Pmin", "Pgene", "Gene-Level FDR")
write.table(temp, "/Users/ben/Google_Drive/Documents/Peroxide/Peroxide Paper/Table_1.txt", sep='\t', row.names=F, quote=F)
################################################################################################

plink[plink$P < 10^-5, ] -> sigs
annot[annot$ID %in% sigs$SNP, ] -> x
x[c(1,4,6:17)] -> x
head(x)
gather(x, key=gene.number, val=gene, -ID, -GeneAnnotation_2, factor_key = F) -> p
head(p)
p[p$gene.number == 'Gene1' | !is.na(p$gene), ] -> pp
pp$GeneAnnotation_2[2]
cSplit(pp, 'GeneAnnotation_2', sep="(") -> ppp
head(ppp)
cSplit(ppp, 'GeneAnnotation_2_01', sep="[") -> p4
head(p4)
p4$GeneAnnotation_2_01_2 -> pp$location
str(pp)
pp[ ,c(1,4,5)] -> snp.annot.location

merge(sigs, pp, by.x = 'SNP', by.y= 'ID') -> snp.table
str(snp.table)
snp.table[ ,c(1, 9, 13, 14)]  -> snp.table
snp.table
rm(p)
rm(ppp)
rm(p4)

library(AnnotationDbi)
library(org.Dm.eg.db)
snp.table$symbol <- mapIds(org.Dm.eg.db, 
                           keys=snp.table$gene, 
                           column="SYMBOL", 
                           keytype="ENSEMBL",
                           multiVals="first")

snp.table$name =   mapIds(org.Dm.eg.db,
                          keys=snp.table$gene, 
                          column="GENENAME",
                          keytype="ENSEMBL",
                          multiVals="first")

snp.table
as.data.frame(snp.table) -> x
str(x)

write.csv(snp.table, file="/Users/ben/Google_Drive/Documents/Peroxide/Peroxide Paper/Table S1_from 11PCs.csv", quote=F, row.names=F) 

################################################################################################



################################################################################################
##### Plot Gene Region and Show Gene Diagrams ######
################################################################################################
head(plink)
plink[plink$CHR == 2 & plink$BP > 490000 & plink$BP < 560000, ] -> ush.range
head(ush.range)
ush.range$BP -> ush.coordinates
plot(-log(ush.range$P, 10) ~ BP, data=ush.range, pch=19, ylim=c(0,8))

## add gene diagrams
### you've got to download the exon coordinates from the UCSC genome browser tool (http://genome.ucsc.edu/cgi-bin/hgGateway?org=D.+melanogaster&db=0):  
### Brief Instructions:
# zoom to a region in the browser that includes all of the gene(s) you'd like diagrams for
# at the top:  Tools > Table Browser
# select 'position'
# under 'output format' select 'BED - browser extensible data'
# click 'get output'
# check 'include custom track header:'
# select 'exons plus'
# click 'get BED'
#### you'll need to reformat the resulting table.

read.table("/Users/ben/Google_Drive/Documents/Peroxide/Results from online GWAS/ush.region.exons.csv", header=T, sep=',', stringsAsFactors = F) -> exons
head(exons)
range(exons$stop)
table(exons$database.name, exons$orientation)

quartz(height = 3, width = 4)
plot(-log(ush.range$P, 10) ~ BP, data=ush.range, pch=19, ylim=c(-3,8), las=1, cex.axis=1, cex=0.5, xlab='Ch 2 - positon (nt)', ylab='P value (-log10)', cex.lab=1, cex.axis=0.5) 

## list of the genes in df.exons
unique(exons$database.name) -> genes.in.interval

### at this point, I resorted to pulling the translations of the ref.seq codes into gene names from the Flybase website
read.table("/Users/ben/Google_Drive/Documents/Peroxide/Results from online GWAS/ush.region.gene.codes.txt", header=T, stringsAsFactors = F) -> gene.codes
gene.codes

## add gene names to select mRNAs:
exons$gene.name <- NA
for(i in 1: length(gene.codes$ref.seq)) {
exons$gene.name[grepl(gene.codes$ref.seq[i], exons$database.name)] <- gene.codes$name[i]
} # done!

exons[!is.na(exons$gene.name), ]

plot(-log(ush.range$P, 10) ~ BP, data=ush.range, pch=19, ylim=c(-3,8), las=1, cex.axis=1, cex=0.5, xlab='Ch 2 - positon (nt)', ylab='P value (-log10)', cex.lab=1, cex.axis=0.5) 

## draw all genes:
for(g in 1:length(unique(gene.codes$name))) {
  gene.codes$name[g] -> gene 
  sample(mycol4, 1) -> color
exons[exons$gene.name == gene, ] -> target.exons
target.exons[!is.na(target.exons$start), ] -> target.exons
ifelse(target.exons$orientation[1] == '-', j <- 0.7, j <- 0)
target.exons$start -> starts
target.exons$stop -> stops
lines(c(starts[1], starts[nrow(target.exons)]), c(-1.2-j, -1.2-j), col=1, lwd=1)
for (i in 1:nrow(target.exons)) {
  rect(starts[i], -1.4-j, stops[i], -1-j , border = 1, col = color, lwd = 0)
} }
text(515000, -0.7, labels = 'ush', font=3, cex=0.5)
text(546000, -0.7, labels = 'AChRB3', font=3, cex=0.5)
text(544000, -2.5, labels = 'spp', font=3, cex=0.5)
text(541000, -2.7, labels = 'lwr', font=3, cex=0.5)
text(550000, -2.7, labels = 'Ets21C', font=3, cex=0.5)
text(558000, -2.5, labels = 'rampA', font=3, cex=0.5)

quartz.save("/Users/ben/Google_Drive/Documents/Peroxide/Results from online GWAS/ush.region.plot.jpg", type = "jpg", device = dev.cur(), dpi = 300)

##################################################################################################################################

