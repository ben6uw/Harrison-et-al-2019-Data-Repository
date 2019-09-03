#Analyzing the master data frame made by William Gordon that contians sumary stats for the peroxide Trials: 1 through 20 (no #7 and #19) and 'Extremes1'.  This is a total of 19 trials:
library(data.table)
library(splitstackshape)
library(survival) 
library(dplyr)
library(tidyr)
library(plyr)
library(lme4)
library(arm) # for display of modeling results
library(lmerTest) # this modifies the summery() function (adds a Pvalue, which may be incorrectly estimated).
library(psych) # allows ICC()
library(smatr) # for major-axis regression
library(readxl)

read.csv("/Users/ben/Google_Drive/Documents/Peroxide/master.peroxide.data.csv", row.names=NULL) -> master
str(master)
table(master$treatment)
levels(master$treatment) <- c('control', 'control', 'H2O2') 
table(master$treatment)
subset(master, treatment == 'H2O2') -> temp
table(temp$treatment)
table(temp$trial)
head(temp)
master$trials <- factor(master$trial, levels = c('1','2','3','4','5','7','Extremes1','8','9','10','11','12','13','14','15','16','17','18','20'))

levels(master$trials)
table(master$trials, master$treatment)

#I would like to test the effect of the switch from Debora Xi (DX) food to William Gordon (WG) Food.  I think that trials 13 through 20 were done with WG food, and all of the rest, including 'extremes' were done with DX food.

master$food <- as.factor(ifelse(master$trials == '1'| master$trials =='2'| master$trials =='3'| master$trials =='4'| master$trials =='5'| master$trials =='6'| master$trials =='Extremes'| master$trials =='8' | master$trials =='9'| master$trials =='10' | master$trials =='11' | master$trials =='12' , 'DX', 'WG'))

str(master)
names(master)

master[master$treatment == 'H2O2', ]-> treats
levels(treats$trial) = c('1','2','3','4','5','6','Extremes1','8','9','10','11','12','13','14','15','16','17','18','20')
str(treats)

par(mar=c(6,6,6,6))
par(mfrow=c(1,1))
boxplot(meanLifespanHrs ~ trials, data = treats, col = 4+c(1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2), main = 'Peroxide Lifespan by Block\ncolored by food recipe (as far as I can tell)\n trials are arranged from earliest to latest', xlab = 'trial', ylab = 'lifespan (h)', cex.lab=2, cex.main=2, outline=F)
stripchart(treats$meanLifespanHrs ~ treats$trials, vertical = TRUE, method = "jitter", add = TRUE, pch = 19, col = "black", cex=0.5)

# the df. master has data from all blocks except 'metabolomics' and 'trial 21'.  As far as I can tell, it also has the weight data for each line in EACH BLOCK!  I will do a lm for lifespan ~ weight + block:
names(master)
master[master$treatment=='H2O2', ] -> treats

plot(meanLifespanHrs ~ meanWeightThisTrial, data=treats)
stripchart(meanLifespanHrs ~ meanWeightAllTrials, data=treats, vertical=T, pch=19, cex=0.5)
fit <- lm(meanLifespanHrs ~ meanWeightThisTrial + line, data=treats)
summary(fit)
anova(fit) -> anov
anov$`Sum Sq` -> vars
vars[1]/sum(vars)
vars[2]/sum(vars)

head(treats)

anova(lm(meanLifespanHrs ~ meanWeightThisTrial, data=treats)) 

anova(lm(meanLifespanHrs ~ meanWeightThisTrial + line, data=treats)) -> mod
mod$`Sum Sq`
mod$`Sum Sq`/sum(mod$`Sum Sq`)

# to correct for weight effect, get residuals from lm(lifespan ~ weight):
treats[!is.na(treats$meanWeightThisTrial), ] -> treats
lm.wt <- lm(meanLifespanHrs ~ meanWeightThisTrial, data=treats)
treats$weight.residuals <- lm.wt$residuals
par(mfrow=c(2,2))
plot(treats$meanLifespanHrs ~ treats$meanWeightThisTrial)
abline(lm(treats$meanLifespanHrs ~ treats$meanWeightThisTrial))
plot(treats$weight.residuals ~ treats$meanWeightThisTrial)
abline(lm(treats$weight.residuals ~ treats$meanWeightThisTrial))
hist(treats$meanLifespanHrs, col='indianred')
hist(treats$weight.residuals, col='indianred')

#now that I have a weight residual from each line in each block in which that line had a weight taken, I will need to redo the work below (which was done earlier) to analyze the effect of line on weight.residuals in each block (instead of the effect of line on lifespan in each block).  I start this new effort after the code below potted the field of histograms (~ line 530)


#Try ANOVA to estimate heritablity (broad sense) across trials: each trial has an estimate of heritability (proportion of total variance explained by line), so get the heritability for all trials and then average those hertitabilities.
# Format for one-way ANOVA in Randomized Block Design ('trial' is the blocking factor):
fit <- lm(meanLifespanHrs ~ line + trial, data=treats)
summary(fit)
af <- anova(fit)
af
afss <- af$"Sum Sq"
print(round((cbind(af,PctExp=afss/sum(afss)*100)), 2))

# This is an over estimate of heritability, because the line means are hiding the within line variation.
# get lifespans for individual flies (AgeH for each death) in each block.
read.table('/Users/ben/Google_Drive/Documents/Peroxide/Peroxide Trials/dLifeFilesNoExpCensor/dLife.files.all.trials.csv', header=T, stringsAsFactors = F, sep=',') -> x

x[x$treatment == 'H2O2', ] -> x
head(x)
unique(x$line) # need to separate line (Ral and RAL !):
length(unique(x$line))
transform(x, temp.line = substr(line, 4, 6)) -> x
length(unique(x$temp.line))
x$line <- paste('Ral', x$temp.line, sep="_")
head(x)
length(unique(x$line))
table(x$line)
x$line.block <- paste(x$line, x$trial)
range(table(x$line.block))
hist(table(x$line.block), breaks=30, col='grey')
data.frame(table(x$line.block)) -> p
head(p)
colnames(p) <- c('line.block', 'total.flies.all.blocks')
merge(x, p, by='line.block') -> y
head(y)
hist(y$total.flies.all.blocks, breaks=30, col='grey') # now have df. that gives the number of flies per line measures across all blocks.  I should exclude lines without at least 30  (?) flies measured...
y[y$total.flies.all.blocks >= 30, ] -> y
hist(y$total.flies.all.blocks, breaks=30, col='grey') 

#### limit heritability to only those lines measured in 15+line blocks
rbind(t18.2, t17.2, t11.2, t8.2, t6.2, t5.2, t4.2, t2.2, t1.2) -> fifteen.plus.line.blocks

data.frame(fifteen.plus.line.blocks) -> x
x[x$treatment == 'H2O2', ] -> x
head(x)
unique(x$line) # need to separate line (Ral and RAL !):
length(unique(x$line))
transform(x, temp.line = substr(line, 4, 6)) -> x
length(unique(x$temp.line))
x$line <- paste('Ral', x$temp.line, sep="_")
head(x)
length(unique(x$line))
table(x$line)
x$line.block <- paste(x$line, x$trial)
range(table(x$line.block))
hist(table(x$line.block), breaks=30, col='grey')
data.frame(table(x$line.block)) -> p
head(p)
colnames(p) <- c('line.block', 'total.flies.all.blocks')
merge(x, p, by='line.block') -> y
head(y)
hist(y$total.flies.all.blocks, breaks=30, col='grey') # now have df. that gives the number of flies per line measures across all blocks.  I should exclude lines without at least 30  (?) flies measured...
y[y$total.flies.all.blocks >= 30, ] -> y
hist(y$total.flies.all.blocks, breaks=30, col='grey') 



############# two approaches to esitmate heritability:
# First approach:
# 3_22_2018 Daniel pointed out that Trudy MacKay uses a different measure of heritability (I think she uses even more than this). In Garlapow 2015 PMC4574202, they calu heritability as H2 = vL / vL + vE, where vL is the among line variance and vE is the within line variance.  

head(y)
str(y)
length(unique(y$line))
length(unique(y$line.block))
hist(y$Deaths)
par(mfrow=c(1,1))
hist(y$AgeH, main='All deaths\nblocks with >15 lines\nlines with > 30 flies', breaks=30, col=1)
hist(scale(y$AgeH), main='Scaled - All deaths\nblocks with >15 lines\nlines with > 30 flies', breaks=30, col=1)



stats <- print(survfit(Surv(AgeH) ~ y$line.block, data=y), print.rmean=T) 
as.data.frame(summary(stats)$table) -> stats

summary(aov(AgeH ~ line * trial, data = y))
y -> temp

mvtab2 <- ddply(temp, .(line, trial), summarise, mean = mean(AgeH), var = var(AgeH))
head(mvtab2)



# calculate mean and variance of lifespans by line and trial:
mvtab2 <- ddply(temp, .(line, trial), summarise, mean = mean(AgeH), var = var(AgeH))
head(mvtab2)

aggregate(mvtab2$var, by=list(mvtab2$trial), mean) -> Ve.by.trial
colnames(Ve.by.trial) <- c('trial', 'Ve')
Ve.by.trial
head(mvtab2)
aggregate(mvtab2$mean, by=list(mvtab2$trial), var) -> Vg.by.trial
Vg.by.trial
colnames(Vg.by.trial) <- c('trial', 'Vg')
merge(Vg.by.trial, Ve.by.trial, by='trial') -> vees
vees
vees$H2 <- vees$Vg/(vees$Vg + vees$Ve)
vees
mean(vees$H2) # heritability = 0.5112559
sd(vees$H2) # sd = 0.08892838
sd(vees$H2)/sqrt(nrow(vees)) # se = 0.0296


# Second approach:

# ANOVA to find prop of var by line within a trial:
fit <- lm(AgeH ~ as.factor(line), data=t1.2)
summary(fit)
af <- anova(fit)
af
afss <- af$"Sum Sq"
print(round((cbind(af,PctExp=afss/sum(afss)*100)), 2))
afss/sum(afss) -> percent.var
percent.var[1]

table(all.blocks$trial) # 19 trials for Heritability analysis (those trials where lines were selected 'at random' from the DGRP, e.g. not the trials called 'extremes' or 'metabolomics')
class(all.blocks)
as.data.frame(all.blocks) -> all.blocks.2

x = rep(0,9) # the number of trials with at least 18 lines (N=18 to 35 in this group of blocks)
for (i in c(1, 2, 4, 5, 6, 8, 11, 17, 18) )
{
  df <- subset(fifteen.plus.line.blocks, trial==i)
  ano <- anova(lm(AgeH ~ as.factor(line), data=df))
  afss <- ano$"Sum Sq"
  afss/sum(afss) -> vars
  vars[1] -> var.due.to.line
  x[i] = var.due.to.line
}
x
x[ x > 0.01] -> x # the loop above calculated heritability for the two trials missing from the series 1:21, so this code gets rid of those results
x[!is.na(x)] -> x
x
mean(x) # mean heritiability
sd(x)/sqrt(length(x)) # standard error for heritability
hist(x, breaks=10)
# heritability = 49.2 +/- 2.3%

shapiro.pvalue = rep(0,19) # the number of trials
shapiro.statistic = rep(0,19) # the number of trials
for (i in c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 18, 20, 21) )
{
  df <- subset(all.blocks.2, trial==i)
  norm <- shapiro.test(df$AgeH)
  norm$p.value -> shap.pees
  norm$statistic -> shap.stats
  shapiro.pvalue[i] =  shap.pees
  shapiro.statistic[i] =  shap.stats
}
shapiro.pvalue
shapiro.pvalue[ shapiro.pvalue > 6.387752e-30] -> shapiro.pvalue
shapiro.pvalue
shapiro.statistic
shapiro.statistic[ shapiro.statistic > 0.001] -> shapiro.statistic
shapiro.statistic

# so, AgeH ('time of death') is not normal across any if the whole blocks, BUT!... it really only needs to be normal for each line wihtin the block. REGARDLESS: as long as the residuals of the lm are homoscedastic, then its A-O-K.

lm(AgeH ~ as.factor(line), data=t2.2) -> fit
shapiro.test(fit$residuals)
hist(fit$residuals)


hist(all.blocks.2$AgeH, main='Age at death, ALL FLIES, ALL BLOCKS', col='lightblue')


par(mfrow=c(4,5))
for (i in c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 18, 20, 21) )
{
  df <- subset(all.blocks.2, trial==i)
  print(hist(df$AgeH), breaks=10, main=i)
}

par(mfrow=c(4,5))
for (i in c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 18, 20, 21) )
{
  df <- subset(all.blocks.2, trial==i)
  lm(AgeH ~ as.factor(line), data=df) -> fit
  print(hist(fit$residuals), breaks=10)
}

# Analyze the effect of line on weight.residuals:

#combine line and trial to give line.block, which I can use to merge master$weight.resids into df.all.blocks.2
unique(all.blocks.2$line)
transform(all.blocks.2, Ral = substr(line, 1, 3), numba = substr(line, 4, 6)) -> all.blocks.3
unique(all.blocks.3$numba)
all.blocks.3$line <- paste('Ral', all.blocks.3$numba, sep='')
length(unique(all.blocks.3$line)) # looks like data for 192 lines is in all.batches.3
length(unique(treats$line)) # data for only 165 lines is in treats[any line that did not have a weight was not included]

all.blocks.3$line.block <- paste(all.blocks.3$line, all.blocks.3$trial)
unique(all.blocks.3$line.block)
treats$line.block <- paste(treats$line, treats$trial)
unique(treats$line.block)

merge(all.blocks.3, treats, by = 'line.block') -> all.blocks.4

unique(all.blocks.4$trial.x)
str(all.blocks.4)

# count the number of lines per trial:
as.numeric(unique(all.blocks.4$trials)) -> trials
trials[!is.na(trials)] -> trials 
head(all.blocks.4)
par(mfrow=c(1,1))
stripchart(meanLifespanHrs ~ trials, data=all.blocks.4, method='jitter', vertical=T, pch=1)

ddply(df, .(color), mutate, count = length(unique(type)))

data.frame(all.blocks.4[ ,c(1)]) -> blocks
head(blocks)
names(blocks) <- 'line.block'
separate(blocks, line.block, into=c('line', 'block'), sep=' ', remove=T) -> blocks
head(blocks)
table(blocks$block)
unique(blocks) -> blocks

# trial 19 should be excluded too:
trials[trials!=19] -> trials 
trials
z = rep(0,14) # the number of trials
for (i in trials)
{
  df <- subset(all.blocks.4, trials==i)
  ano <- anova(lm(AgeH ~ as.factor(line.x), data=df))
  afss <- ano$"Sum Sq"
  afss/sum(afss) -> vars
  vars[1] -> var.due.to.line
  z[i] = var.due.to.line
}
z
z[ z > 0.01] -> z # the loop above calculated heritability for trials that didn't exist (these were 0.000), so this code gets rid of those results
z
mean(z)
sd(z)/sqrt(14) # se of the mean heritability 



