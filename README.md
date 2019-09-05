# Harrison-et-al-2019-Data-Repository
Data for The metabolome as a link in the genotype-phenotype map for peroxide resistance in the fruit fly, Drosophila melanogaster (2019)

Authors: Benjamin R. Harrison, Lu Wang, Erika Gajda, Elise V. Hoffman, Brian Y. Chung, Scott D. Pletcher, Haiwei Gu, Dan Raftery and Daniel E.L. Promislow

requires R (https://www.r-project.org/)

List of folders and files in /data:
# DGRP_peroxide_trials
    Lifespan data for DGRP lines on peroxide food, with R code
# GAL4_geneswitch_UAS_RNAi
    Lifepan data for RNAi experiments on peroxide food, with R code
# GWAS
    R code for GWAS analysis, including purmutation approach to get P-gene

    requires: 
    PLINK (http://zzz.bwh.harvard.edu/plink/)
    BED files (http://dgrp2.gnets.ncsu.edu/data/website/dgrp2.bed) This is further modified with PLINK to make a BED file that only contains the lines in the study and the variants that pass the MAF and missingness thresholds. The plink call: plink -bfile (path to existing .bed file) --make-bed --maf 0.05 --geno 0.3 --keep-fam (list of FID to retain) --out 'name of output .bed file'
    VCF file for the DGRP Freeze 2.0 calls (http://dgrp2.gnets.ncsu.edu/data/website/dgrp2.vcf) These are used for LD analysis, you will benifit from converting the genotypes into a numeric matrix
    Variant annotation files (http://dgrp2.gnets.ncsu.edu/data/website/dgrp.fb557.annot.txt) for the DGRP 
# NPF.survival.H2O2
    Survival data for NPF mutant
# feeding_rate_experiments
    Feeding rate data for flies on peroxide food, with R code
# supplemental_carbohydrate_experiment
    Supplemental carbohydrate data, with R code
