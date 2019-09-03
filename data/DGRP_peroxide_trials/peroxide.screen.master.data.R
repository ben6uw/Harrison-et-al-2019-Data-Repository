#Analyzing the master data frame made by William Gordon that contians sumary stats for the peroxide Trials: 1 through 20 (no #7 and #19) and 'Extremes1'.  This is a total of 19 trials:
library(dplyr)
library(plyr)
library(lme4)
install.packages('arm')
library(arm) # for display of modeling results
install.packages('lmerTest')
library(lmerTest) # this modifies the summery() function (adds a Pvalue, which may be incorrectly estimated).


read.csv("/Users/ben/Google_Drive/Documents/Peroxide/master.peroxide.data.csv", row.names=NULL) -> master
levels(master$treatment) <- c('control', 'control', 'H2O2') 
table(master$treatment)
table(master$trial)
master$trials <- factor(master$trial, levels = c('1','2','3','4','5','7','Extremes1','8','9','10','11','12','13','14','15','16','17','18','20'))
table(master$trial, master$trials)
master$trial <- master$trials

#I would like to test the effect of the switch from Debora Xi (DX) food to William Gordon (WG) Food.  I think that trials 13 through 20 were done with WG food, and all of the rest, including 'extremes' were done with DX food.

master$food <- as.factor(ifelse(master$trials == '1'| master$trials =='2'| master$trials =='3'| master$trials =='4'| master$trials =='5'| master$trials =='6'| master$trials =='Extremes'| master$trials =='8' | master$trials =='9'| master$trials =='10' | master$trials =='11' | master$trials =='12' , 'DX', 'WG'))
table(master$trials, master$food)
master[ ,-c(15)] -> master # gets rid of 'trials', leaving 'trial.
unique(master$food)
master[master$treatment == 'H2O2', ]-> treats
str(treats)

treats[treats$trial != 'Extremes1', ] -> treats # remove the 'extremes' trial prior to anova or mixed modeling of lifespan ~ _____, ______ etc.
treats[!is.na(treats$trial), ] -> treats 
table(treats$trial)

## use mixed model test to test the fixed effect of food: Lifespan ~ food + (1 | block) + (1 | genotype) + ϵ.
summary(lmer(log(meanLifespanHrs) ~ meanWeightThisTrial + food + (1|trial) + (1|line), treats) -> food.mod)
rand(food.mod)
plot(food.mod, pch=19, col=1, cex=0.5)

length(unique(treats$line))


# the df. master has data from all blocks except 'metabolomics' and 'trial 21'.  As far as I can tell, it also has the weight data for each line in EACH BLOCK!  I will do a lm for lifespan ~ weight + block:
par(mar=c(5, 4, 4, 2))
plot(meanLifespanHrs ~ meanWeightThisTrial, data=treats, pch=19, cex=0.5, las=1, cex.axis = 1.3, xlab='', ylab='')
summary(aov(meanLifespanHrs ~ meanWeightThisTrial + line, data=treats))

summary(lmer(meanLifespanHrs ~ meanWeightThisTrial + (1|line), treats) -> m1)
rand(m1)

# plot LS ~ weight for only those lines used for metabolomics:
head(treats)

sample.info <- read.csv("/Users/ben/Google_Drive/Documents/Peroxide/Lu Wangs Analysis/Bens.Analysis/sample.info.2.csv", header=TRUE, stringsAsFactors = F)
unique(sample.info$line) -> line.names
line.names[-1] -> line.names
gsub('ral-', 'Ral', line.names) -> metabolome.lines 

treats[treats$line %in% metabolome.lines, ] -> metabolome.lines
plot(meanLifespanHrs ~ meanWeightThisTrial, data=metabolome.lines, pch=19, cex=0.5, las=1, cex.axis = 1.3, xlab='', ylab='')

# to correct for weight effect, get residuals from lm(lifespan ~ weight):
treats[!is.na(treats$meanWeightThisTrial), ] -> treats
lm.wt <- lm(meanLifespanHrs ~ meanWeightThisTrial, data=treats)
summary(lm.wt)
treats$weight.residuals <- lm.wt$residuals
par(mfrow=c(2,2))
plot(treats$meanLifespanHrs ~ treats$meanWeightThisTrial)
abline(lm(treats$meanLifespanHrs ~ treats$meanWeightThisTrial))
plot(treats$weight.residuals ~ treats$meanWeightThisTrial)
abline(lm(treats$weight.residuals ~ treats$meanWeightThisTrial))

head(treats)

table(treats$trial)
#exclude trials with fewer than 7 lines from this analysis
treats[treats$trial != '13', ] -> treats
treats[treats$trial != '14', ] -> treats
treats[treats$trial != '16', ] -> treats
table(treats$trial)

par(mar=c(6,6,6,6))
par(mfrow=c(1,1))
boxplot(meanLifespanHrs ~ trials, data = treats, col = 4+c(1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2), main = 'Peroxide Lifespan by Block\ncolored by food recipe (as far as I can tell)\n trials are arranged from earliest to latest', xlab = 'trial', ylab = 'lifespan (h)', cex.lab=2, cex.main=2, outline=F)
stripchart(treats$meanLifespanHrs ~ treats$trials, vertical = TRUE, method = "jitter", add = TRUE, pch = 19, col = "black", cex=0.5)

########################################################################
## is there an effect of food on the weight.residuals?

boxplot(weight.residuals ~ food, data = treats, outline=F, ylim=c(-35, 55), las=1, ylab='weight residuals', xlab='food', cex.lab=2, col='grey')
stripchart(weight.residuals ~ food, data = treats, add=T, vertical=T, method='jitter', ylim=c(-35, 55), pch=19)

boxplot(meanLifespanHrs ~ food, data = treats, outline=T, las=1, ylab = 'mean lifespan', xlab='food', cex.lab=2, col='grey')
stripchart(meanLifespanHrs ~ food, data = treats, add=T, vertical=T, method='jitter', ylim=c(-35, 55), pch=19)

## is there an effect of food or food:line after block-centering?

tapply(treats$weight.residuals, list(treats$trial), scale) -> x
ddply(treats, c("trial"), transform, scaled.wt.resid = scale(weight.residuals)) -> treats

hist(treats$scaled.wt.resid, col=8)
summary(aov(weight.residuals ~ food + trial + line + food:line, data = treats))
summary(aov(scaled.wt.resid ~ food + trial + line + food:line, data = treats))
## it appears that, block scaling eliminates any main effect of food, however there is still a singificant interaction between line and food - some lines have different survivals on the different foods.


summary(lmer(weight.residuals ~ food + (1 | trial), data = treats)) #testing the effect of food as a fixed effect, and trial as a random effect
summary(lmer(weight.residuals ~ food + (1 | trial) + (1 | line), data = treats))  #testing the effect of food as a fixed effect, and trial and genotype are a random effects
summary(lmer(meanLifespanHrs ~ food + (1 | trial) + (1 | line), data = treats))  #testing the effect of food as a fixed effect, and trial and genotype are a random effects

########################################################################



########################################################################
######### repeatability #############:
head(treats)
table(treats$line, treats$trial) 

## find the average weight.resid for each line by trial - doing this because some lines have more than one value per block (e.g. in trial two there were two groups, purple and blue, that were measured).
j <- treats %>% dplyr::group_by(line, trial) %>% dplyr::summarise(per.ls = mean(weight.residuals))
data.frame(j) -> k
reshape(k, idvar = "line", timevar = "trial", direction = "wide") -> m
rownames(m) <- m$line
m[ ,-c(1) ] -> m

ICC(m, missing=FALSE, alpha=.05) # missing should be FALSE, otherwise (from the white paper: 'if TRUE, remove missing data – work on complete cases only') and df.z has no complete cases!
ICC(m, missing=TRUE, alpha=.05)

# the ICC is an anova, which relies on the data (lifespans) being normal.   Are they?
n <- gather(m) # make a column of all lifespans
hist(n$value, breaks=15, col='grey') # pretty much


## the same analysis, but on the wt.residuals scaled by block - this may help us see how block effects affect ICC/repeatability
x <- treats %>% dplyr::group_by(line, trial) %>% dplyr::summarise(per.ls = mean(scaled.wt.resid))
data.frame(x) -> y
reshape(y, idvar = "line", timevar = "trial", direction = "wide") -> z

# how many measurements per line:
rownames(z) <- z$line
z[ ,-c(1) ] -> z
apply(z, 1, function(x)length(unique(x)))-1 -> replicates

# eliminate lines measured only once:
names(replicates [replicates > 1] ) -> met.the.cuttoff
z[met.the.cuttoff, ]-> zz

dim(zz) # 15 blocks ('judges') ruling on the scaled lifespans (wt.residuals) of 43 lines
ICC(zz, missing=FALSE, alpha=.05) # missing should be FALSE, otherwise (from the white paper: 'if TRUE, remove missing data – work on complete cases only') and df.z has no complete cases!
ICC(zz, missing=TRUE, alpha=.05) # there is no difference between missing T and missing F....?
ICC(zz, missing=TRUE, alpha=.05) -> t
t$summary # shows the anova model used for the ICC


names(replicates [replicates == 2] ) -> measured.twice
z[measured.twice, ]-> dups
head(dups)
t(dups) -> dups2
dups2[sapply(dups2, function(x) !(is.na(x)))] -> dups3 ## a vector where every 2 measurements are from the same line
dup4 <- as.data.frame(matrix(dups3, nrow=39, ncol=2, byrow=T))
library(smatr)
plot(dup4$V1 ~ dup4$V2, pch=19, las=1, main='lines measured twice\n(n=39, major axis R^2 = 0.30)', xlab='scaled wt.residual.1', ylab='scaled.wt.residual.2')
abline(sma(dup4$V1 ~ dup4$V2), lty=2)

perm.results <- numeric()
for(i in 1:10000) {
dup4$perm <- sample(dup4$V2, replace=F)
unlist(sma(dup4$V1 ~ dup4$perm)$r2) -> perm.results[i]
}
hist(perm.results, breaks=20, col='grey', xlab='r2 from permutation', main='permutation', cex.main=2)

unlist(sma(dup4$V1 ~ dup4$V2)$r2) -> actual
table(perm.results > actual) # what proportion of permutation r2 values are greater than the measured value?


########################################################################
## The analysis above was done on lifespans measured in blocks with at least 7 lines, however we ultimately depend on the accuracy of block-scaled values to estimate repeatability.  I should examine the effect of a higher cut-off (e.g at least 11 lines per block):

table(treats$trial)
#exclude trials with fewer than 7 lines from this analysis
treats[treats$trial != '13', ] -> treats
treats[treats$trial != '14', ] -> treats
treats[treats$trial != '16', ] -> treats
#exclude trials with fewer than 11 lines from this analysis
treats[treats$trial != '10', ] -> treats
treats[treats$trial != '12', ] -> treats
treats[treats$trial != '15', ] -> treats
treats[treats$trial != '20', ] -> treats
treats[treats$trial != '3', ] -> treats

table(treats$trial)

## the same analysis, but on the wt.residuals scaled by block - this may help us see how block effects affect ICC/repeatability
x2 <- treats %>% dplyr::group_by(line, trial) %>% dplyr::summarise(per.ls = mean(scaled.wt.resid))
data.frame(x2) -> y2
reshape(y2, idvar = "line", timevar = "trial", direction = "wide") -> z2

# how many measurements per line:
rownames(z2) <- z2$line
z2[ ,-c(1) ] -> z2
apply(z2, 1, function(x)length(unique(x)))-1 -> replicates2
hist(replicates2)

## subset down to lines measured twice
names(replicates2 [replicates2 == 2] ) -> measured.twice2
z2[measured.twice2, ]-> dups2.1
head(dups2.1)
t(dups2.1) -> dups2.2
dups2.2[sapply(dups2.2, function(x) !(is.na(x)))] -> dups2.3 ## a vector where every 2 measurements are from the same line
dup2.4 <- as.data.frame(matrix(dups2.3, nrow=24, ncol=2, byrow=T))
library(smatr)
plot(dup2.4$V1 ~ dup2.4$V2, pch=19, las=1, main='lines measured twice\n(n=24, major axis R^2 = 0.28)', xlab='scaled wt.residual.1', ylab='scaled.wt.residual.2')
sma(dup2.4$V1 ~ dup2.4$V2)$r2
abline(sma(dup2.4$V1 ~ dup2.4$V2), lty=2)

perm.results2 <- numeric()
for(i in 1:10000) {
  dup2.4$perm <- sample(dup2.4$V2, replace=F)
  unlist(sma(dup2.4$V1 ~ dup2.4$perm)$r2) -> perm.results2[i]
}
hist(perm.results2, breaks=20, col='grey', xlab='r2 from permutation', main='permutation', cex.main=2)

unlist(sma(dup2.4$V1 ~ dup2.4$V2)$r2) -> actual2
table(perm.results2 > actual2) # what proportion of permutation r2 values are greater than the measured value?




