mycol <- c("#8B8378", "#EE9A00") 
mycol2 <- c(mycol,"#483D8B")
mycol4 = c("#EEC900", "#CD6600", "#87CEFA", "#104E8B")
mycol3 <- mycol4

### Blue dye assay with and without peroxide treatment (3/1/2018):
read.csv("../3-1-2018.blue.feeding.assay.with.peroxide.treatment.csv", header=T) -> blue
colnames(blue) [6] <- 'well'
as.character(blue$well) -> blue$well
blue[ ,7:47] / blue$number.of.flies -> blue [ ,48:88] # corrects the Abs values for the diffenet numbers of flies per sample


simp <- blue[ , c(1,2,3,4,40,81)] 
simp <- simp[simp$line != 'blank', ]
simp <- simp[simp$time.without.dye == 0 | simp$time.without.dye == 4, ] 
simp$line.time <- as.factor(paste(simp$line, simp$time.without.dye))
simp$line.time.treatment <- as.factor(paste(simp$line, simp$time.without.dye, simp$treatment))
simp$line.time.treatment = factor(simp$line.time.treatment,levels(simp$line.time.treatment)[c(2,4,1,3,6,8,5,7,10,12,9,11)]) # order factor levels in a better way for plotting

plot(X630.1 ~ line.time.treatment, data=simp, col=mycol4[1+c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4)], cex.axis=0.5, vertical=T, pch=c(rep(c(19,1), 8)), cex=2)
stripchart(X630.1 ~ line.time.treatment, data=simp, cex.axis=0.5, vertical=T, pch=c(rep(c(19,1), 8)), cex=1, add=T)

par(mar=c(8,10,4,4))
simp$line.treatment <- as.factor(paste(simp$line, simp$treatment))
interaction.plot(simp$time.without.dye, simp$line.treatment, simp$X630.1, lwd=2, cex.lab=2, xlab='', cex.axis=2, las=1, ylab='', type='b', pch=rep(c(19, 1), 6), cex=2, col=mycol4[c(2,2,3,3,4,4)], fixed=T, lty=c(1,2,1,2,1,2))
title(ylab="Abs 630nm", line=6, cex.lab=2, xlab="time (min)")


#I would like to normalize to the Abs 630nm of the water samples at time=0.  Since the samples aren't paired, I will need to do this for the average of thewater samples at t=0:
head(simp)
aggregate(simp[, 5:6], list(simp$line.time.treatment), mean) -> temp
colnames(temp) <- c('line.time.treatment', 'mean.630', 'mean.630.1')
temp$mean.zero <- temp$mean.630.1[c(1,1,1,1,5,5,5,5,9,9,9,9)] 
merge(simp, temp, by='line.time.treatment') -> simp
simp$norm.abs <- simp$X630.1/simp$mean.zero

interaction.plot(simp$time.without.dye, simp$line.treatment, simp$norm.abs, lwd=2, cex.lab=2, xlab='', cex.axis=2, las=1, ylab='', type='b', pch=rep(c(19, 1), 6), cex=2, col=mycol4[c(2,2,3,3,4,4)], fixed=T, lty=1)
title(ylab="Normalized Abs 630nm", line=6, cex.lab=2, xlab="time (min)")

plot(norm.abs ~ line.time.treatment, data=simp, col=mycol4[1+c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4)], cex.axis=0.5, vertical=T, pch=c(rep(c(19,1), 8)), cex=2)
stripchart(norm.abs ~ line.time.treatment, data=simp, cex.axis=0.5, vertical=T, pch=c(rep(c(19,1), 8)), cex=1, add=T)

stripchart(norm.abs ~ line.time.treatment, data=simp, cex.axis=0.5, vertical=T, pch=rep(c(19,19,1,1),3), cex=2, col=mycol4[1+c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4)])

i = 'RAL_320'
simp[simp$line == i, ] -> temp
temp$line <- as.factor(temp$line)
temp$line <- factor(temp$line)
temp$line.time.treatment <- factor(temp$line.time.treatment)

stripchart(temp$X630.1 ~ temp$line.time.treatment, cex.axis=0.5, vertical=T, pch=c(15,19,0,1), cex=2, col=mycol4[3], lwd=2, at=1:4, xlim=c(0.5,4.5), las=1, cex.axis=1.5, cex.lab=2, ylab='', xaxt='n')
title(ylab="Abs 630nm", line=4.5, cex.lab=2)
points(temp$mean.630.1 ~ temp$line.time.treatment, pch="-", cex=2) 

i = 'RAL_239'
simp[simp$line == i, ] -> temp
temp$line <- as.factor(temp$line)
temp$line <- factor(temp$line)
temp$line.time.treatment <- factor(temp$line.time.treatment)

stripchart(temp$X630.1 ~ temp$line.time.treatment, cex.axis=0.5, vertical=T, pch=c(15,19,0,1), cex=2, col=mycol4[2], lwd=2, at=1:4, xlim=c(0.5,4.5), las=1, cex.axis=1.5, cex.lab=2, ylab='', xaxt='n')
title(ylab="Abs 630nm", line=4.5, cex.lab=2)
points(temp$mean.630.1 ~ temp$line.time.treatment, pch="-", cex=2) 


i = 'RAL_627'
simp[simp$line == i, ] -> temp
temp$line <- as.factor(temp$line)
temp$line <- factor(temp$line)
temp$line.time.treatment <- factor(temp$line.time.treatment)

temp
water0 <- temp$X630.1[temp$treatment == 'water' & temp$time.without.dye == 0]
perox0 <- temp$X630.1[temp$treatment == 'peroxide' & temp$time.without.dye == 0]
water4 <- temp$X630.1[temp$treatment == 'water' & temp$time.without.dye == 4]
perox4 <- temp$X630.1[temp$treatment == 'peroxide' & temp$time.without.dye == 4]

wilcox.test(water0, perox0)$p.value
wilcox.test(water4, perox4)$p.value
t.test(water0, perox0)
t.test(water4, perox4)

stripchart(temp$X630.1 ~ temp$line.time.treatment, cex.axis=0.5, vertical=T, pch=c(15,19,0,1), cex=2, col=mycol4[4], lwd=2, at=1:4, xlim=c(0.5,4.5), las=1, cex.axis=1.5, cex.lab=2, ylab='', xaxt='n')
title(ylab="Abs 630nm", line=4.5, cex.lab=2)
points(temp$mean.630.1 ~ temp$line.time.treatment, pch="-", cex=2) 

# aggregate means and SD for bar-plotting
temp2 <- aggregate(temp$X630.1,
                   by = list(temp$line.time.treatment),
                   FUN = function(x) c(mean = mean(x), sd = sd(x)))

data.frame(temp2$x) -> temp2
temp2[c(1,3,2,4), ] -> temp2
max(temp2$mean) + max(temp2$sd) -> tippytop

plot <- barplot(height = temp2$mean, col=c(mycol4[4], 0, mycol4[4], 0), lwd=1, las=1, space=c(0.5,0,0.5,0), border=T, ylim = c(0, 1.1*tippytop), ann=F, axes=F)
axis(2,at=seq(0,1.1*tippytop,0.05), line=1.2, cex.axis=1.5)
axis(1,at=seq(1, 4), line=1.2, cex.axis=1.5)
title(ylab="Abs 630nm", line=4.5, cex.lab=2)
title(xlab="time prior to dye (h)", line=5, cex.lab=2)
segments(plot, temp2$mean, y1 = temp2$mean + temp2$sd, lwd = 1.5)
arrows(plot, temp2$mean, y1 = temp2$mean + temp2$sd, lwd = 1.5, angle=90, code=2, length=0.05)

i = 'RAL_239'
simp[simp$line == i, ] -> temp
temp$line <- as.factor(temp$line)
temp$line <- factor(temp$line)
temp$line.time.treatment <- factor(temp$line.time.treatment)

temp
water0 <- temp$X630.1[temp$treatment == 'water' & temp$time.without.dye == 0]
perox0 <- temp$X630.1[temp$treatment == 'peroxide' & temp$time.without.dye == 0]
water4 <- temp$X630.1[temp$treatment == 'water' & temp$time.without.dye == 4]
perox4 <- temp$X630.1[temp$treatment == 'peroxide' & temp$time.without.dye == 4]

wilcox.test(water0, perox0)$p.value
wilcox.test(water4, perox4)$p.value
t.test(water0, perox0)
t.test(water4, perox4)

par(mfrow=c(1,1))
stripchart(temp$X630.1 ~ temp$line.time.treatment, cex.axis=0.5, vertical=T, pch=c(15,19,0,1), cex=2, col=mycol4[2], lwd=2, at=1:4, xlim=c(0.5,4.5), las=1, cex.axis=1.5, cex.lab=2, ylab='', xaxt='n')
title(ylab="Abs 630nm", line=4.5, cex.lab=2)
points(temp$mean.630.1 ~ temp$line.time.treatment, pch="-", cex=2) 

# aggregate means and SD for bar-plotting
temp2 <- aggregate(temp$X630.1,
                   by = list(temp$line.time.treatment),
                   FUN = function(x) c(mean = mean(x), sd = sd(x)))

data.frame(temp2$x) -> temp2
temp2[c(1,3,2,4), ] -> temp2
max(temp2$mean) + max(temp2$sd) -> tippytop

plot <- barplot(height = temp2$mean, col=c(mycol4[2], 0, mycol4[2], 0), lwd=1, las=1, space=c(0.5,0,0.5,0), border=T, ylim = c(0, 1.1*tippytop), ann=F, axes=F)
axis(2,at=seq(0,1.1*tippytop,0.05), line=1.2, cex.axis=1.5)
axis(1,at=seq(1, 4), line=1.2, cex.axis=1.5)
title(ylab="Abs 630nm", line=4.5, cex.lab=2)
title(xlab="time prior to dye (h)", line=5, cex.lab=2)
segments(plot, temp2$mean, y1 = temp2$mean + temp2$sd, lwd = 1.5)
arrows(plot, temp2$mean, y1 = temp2$mean + temp2$sd, lwd = 1.5, angle=90, code=2, length=0.05)


i = 'RAL_320'
simp[simp$line == i, ] -> temp
temp$line <- as.factor(temp$line)
temp$line <- factor(temp$line)
temp$line.time.treatment <- factor(temp$line.time.treatment)

stripchart(temp$X630.1 ~ temp$line.time.treatment, cex.axis=0.5, vertical=T, pch=c(15,19,0,1), cex=2, col=mycol4[3], lwd=2, at=1:4, xlim=c(0.5,4.5), las=1, cex.axis=1.5, cex.lab=2, ylab='', xaxt='n')
title(ylab="Abs 630nm", line=4.5, cex.lab=2)
points(temp$mean.630.1 ~ temp$line.time.treatment, pch="-", cex=2) 

water0 <- temp$X630.1[temp$treatment == 'water' & temp$time.without.dye == 0]
perox0 <- temp$X630.1[temp$treatment == 'peroxide' & temp$time.without.dye == 0]
water4 <- temp$X630.1[temp$treatment == 'water' & temp$time.without.dye == 4]
perox4 <- temp$X630.1[temp$treatment == 'peroxide' & temp$time.without.dye == 4]

wilcox.test(water0, perox0)$p.value
wilcox.test(water4, perox4)$p.value

# aggregate means and SD for bar-plotting
temp2 <- aggregate(temp$X630.1,
                   by = list(temp$line.time.treatment),
                   FUN = function(x) c(mean = mean(x), sd = sd(x)))

data.frame(temp2$x) -> temp2
temp2[c(1,3,2,4), ] -> temp2
max(temp2$mean) + max(temp2$sd) -> tippytop

plot <- barplot(height = temp2$mean, col=c(mycol4[3], 0, mycol4[3], 0), lwd=1, las=1, space=c(0.5,0,0.5,0), border=T, ylim = c(0, 1.1*tippytop), ann=F, axes=F)
axis(2,at=seq(0,1.1*tippytop,0.05), line=1.2, cex.axis=1.5)
axis(1,at=seq(1, 4), line=1.2, cex.axis=1.5)
title(ylab="Abs 630nm", line=4.5, cex.lab=2)
title(xlab="time prior to dye (h)", line=5, cex.lab=2)
segments(plot, temp2$mean, y1 = temp2$mean + temp2$sd, lwd = 1.5)
arrows(plot, temp2$mean, y1 = temp2$mean + temp2$sd, lwd = 1.5, angle=90, code=2, length=0.05)

# ALL Lines together:
# aggregate means and SD for bar-plotting
simp2 <- aggregate(simp$X630.1,
                   by = list(simp$line.time.treatment),
                   FUN = function(x) c(mean = mean(x), sd = sd(x)))

data.frame(simp2$x) -> simp2
simp2[c(1,3,2,4, 4+c(1,3,2,4), 8+c(1,3,2,4)), ] -> simp3
max(simp3$mean) + max(simp3$sd) -> tippytop

plot <- barplot(height = simp3$mean, col=c(mycol4[2], 0, mycol4[2], 0, mycol4[3], 0, mycol4[3], 0, mycol4[4], 0, mycol4[4], 0 ), lwd=1, las=1, space=c(0.5,0,0.5,0), border=T, ylim = c(0, 1.1*tippytop), ann=F, axes=F)
axis(2,at=seq(0,1.1*tippytop,0.05), line=1.2, cex.axis=1.5)
axis(1,at=c(1, 4), line=1.2, cex.axis=1.5)
title(ylab="Abs 630nm", line=4.5, cex.lab=2)
title(xlab="time prior to dye (h)", line=5, cex.lab=2)
segments(plot, simp3$mean, y1 = simp3$mean + simp3$sd, lwd = 1.5)
arrows(plot, simp3$mean, y1 = simp3$mean + simp3$sd, lwd = 1.5, angle=90, code=2, length=0.05)
legend(8, 0.5, c("Ral_239", "Ral_320", "Ral_627", 'peroxide'), pch=c(15, 15, 15, 0), col = c(mycol4[c(2,3,4)], 1), cex=1.5, bty='n')

simp
simp[simp$time.without.dye == 0, ] -> temp
temp2 <- aggregate(temp$X630.1, 
          by= list(temp$line.treatment), 
          FUN = function(x) c('mean' = mean(x), sd =sd(x)))

data.frame(temp2$x) -> temp2
temp2[c(2,1,4,3,6,5), ] -> temp2


plot <- barplot(height = temp2$mean, col=c(mycol4[4], 0), lwd=1, las=1, space=c(0, 0, 0.5, 0, 0.5, 0), border=T, ylim = c(0, 0.40), ann=F, axes=F)
axis(2,at=seq(0,1.1*tippytop,0.05), line=1.2, cex.axis=1.5, las=1)
title(ylab = "Abs 630nm", line=6, cex.lab=2)
title(xlab = c('Ral-239  Ral-320  Ral-627'), cex.lab=2, line=2)

segments(plot, temp2$mean, y1 = temp2$mean + temp2$sd, lwd = 1.5)
arrows(plot, temp2$mean, y1 = temp2$mean + temp2$sd, lwd = 1.5, angle=90, code=2, length=0.05)
legend(1.5, 0.4, c('control', 'peroxide', 'no flies'), pch=c(0, 15, 15), col = c(1, mycol4[4], 'grey'), cex=2, bty='n')

