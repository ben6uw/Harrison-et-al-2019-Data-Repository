
read.table("/Users/ben/Google_Drive/Applications/GBLUP/individual.pheno.with.covars", header = T, stringsAsFactors = T) -> pheno # blocks with non-random selection of DGRP lines were excluded from hertiability analysis to conform with the asumptions of linear modeling

#####################################################################
# Broad Sense Heritabilty, estimated as the among line varaince over the total variance (using ANOVA)
#####################################################################
head(pheno)
lm(log(AgeH) ~ trial.wt, pheno)$residuals -> pheno$wt.residuals

length(unique(pheno$line)) - colSums(table(pheno$line, pheno$trial) == 0) # summary of lines measured per block

lines.per.block.threshold <- 13 # set minimum block size for H2 estimation
length(unique(pheno$line)) - colSums(table(pheno$line, pheno$trial) == 0) # lines per block
length(unique(pheno$line)) - colSums(table(pheno$line, pheno$trial) == 0) >= lines.per.block.threshold -> p
names(p[p]) -> blocks.for.heritability.estimation

H2 <- numeric()
for(i in 1:length(blocks.for.heritability.estimation)){
  pheno[pheno$trial == blocks.for.heritability.estimation[i], ] -> tmp
  anova(lm(wt.residuals ~ line, tmp))$"Sum Sq" -> afss
  (afss/sum(afss))[1] -> H2[i]}

H2 # block-wise estimates of broad sense heritability (H2)
mean(H2) # estimate of H2
sd(H2) /  sqrt(length(H2)) # standard error of H2  


#####################################################################
# estimate SNP-heritability using REML estimates of Vg and Ve based on random effects model of the genome relationship matrix
#####################################################################
library(NAM)


# set up an 'incidence matrix' of the fixed covariates:
fm <- log(AgeH) ~ trial + trial.wt + wolbachia + In.2L.t + In.2R.NS + In.3R.P + In.3R.Mo

X <- model.matrix(fm, data=pheno) # model.matrix takes a formula and makes an incidence matrix
dim(X)

# load relationship matrix 
# this one is based on LD-pruned variants, each at genotpye freq >70%, MAF >5%
read.table('/Users/ben/Google_Drive/Applications/plink_mac/DGRP/prune.rel', header=F, stringsAsFactors = F) -> rel
read.table('/Users/ben/Google_Drive/Applications/plink_mac/DGRP/prune.rel.id', header=F, stringsAsFactors = F) -> rel.id
colnames(rel) <- rel.id$V1
rownames(rel) <- rel.id$V1
as.matrix(rel) -> rel
rel[pheno$line, pheno$line] -> rel

# Use REML to estimate variance components:
start_time <- Sys.time()
reml(log(pheno$AgeH), K=rel, X=X) -> f_r_mod
end_time <- Sys.time()
end_time - start_time # displays the computation time (could take >20min)

f_r_mod$VC # variance components and heritability

#####################################################################

