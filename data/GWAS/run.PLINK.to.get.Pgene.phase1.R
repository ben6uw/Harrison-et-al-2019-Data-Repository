## R code to run PLINK and compile the results of GWAS permutations.   
## This is the first phase, to run 10,000 permutations to derive gene-level scores (Pgene) for genes at a low threshold (Pmin =< 0.01), the resulting Pgene values are used to cutoff a smaller list of genes to run 1,000,000 permutations of the same procedure on.  
## direct files to a working directory
## requires PLINK

read.table("/Users/ben/Google_Drive/Applications/plink_mac/DGRP/genes.to.test.01", stringsAsFactors = F) -> genes.to.test.01 # a list of genes that meet some criteria from a GWAS, in this case, these are genes containing at least one SNP that was significant at P >= 0.01 from a GWAS.
nperm=10000 
numeric() -> gene.score
numeric() -> perm.gene.score

for(i in 1:length(genes.to.test.01$V1)) {
  gene <- genes.to.test.01$V1[i]
  
  command = paste("plink --bfile dgrp2  --maf 0.05 --geno 0.3 --pheno line-matched.peroxide.resistance.GWAS.input.txt --linear --covar plink.eigenvec.txt --covar-number 1-4 --mperm ",
                  nperm,
                  "--mperm-save-all --set cutoff.sets --recode --gene",
                  gene,
                  sep = " ")
  system(command)
  
  read.table("/Users/ben/Google_Drive/Applications/plink_mac/DGRP/plink.mperm.dump.all")-> permutations 
  apply(permutations[1, ], 1, FUN=max, na.rm=T) -> Tmax
  
  permutations[2:nrow(permutations), 2:ncol(permutations)] -> perm.results
  if (is.data.frame(perm.results)) { apply(perm.results, 1, FUN=max, na.rm=T) -> Tmax.permutations }
  else { perm.results -> Tmax.permutations }
  length(Tmax.permutations[Tmax.permutations > Tmax])/length(Tmax.permutations) -> gene.score[i]
}

cbind(genes.to.test.01$V1, gene.score) -> P.gene.phase1.output
write.table( P.gene.phase1.output, file='/Users/ben/Google_Drive/Applications/plink_mac/DGRP/P.gene.phase1.output')

## apply cutoff based on these values and run 1,000,000 permutations on the smaller gene list.
# the code for running the 10^6 permutaitons runs on the terminal (mac) and is in the same GitHub folder as this R file.
