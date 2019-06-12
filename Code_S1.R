rm(list=ls())

####################################################################################
#### 1. PLOTS MAPS FROM FIGURE 1
######################################################################################
library(raster)
library(maps)
library(maptools)
library(mapdata)
library(rgdal)
library(sp)
library(reshape2)
library(mgcv)
library(marmap)

# Load Data
dat <- read.csv(file='Dataset_S1.csv')

# Create the grid
gt<-(GridTopology(c(-60.25, 30.25), c(0.5,0.5), c(200, 120))) 
grt<-SpatialGrid(gt, proj4string=CRS("+init=epsg:4326"))
spix <- as(grt, "SpatialPixels")
spol <- as(spix, "SpatialPolygons")
rnames<-sapply(slot(spol, "polygons"), function(x) slot(x, "ID"))
LOCUNI<-as.data.frame(seq(1,length(spix)))
rownames(LOCUNI)<-rnames
europegrid <-SpatialPolygonsDataFrame(spol, LOCUNI)

### BIOMASS map
dat$biomass <- dat$biomass/1e6
ii <- cut(log(dat$biomass), breaks = seq(log(min(dat$biomass)), log(max(dat$biomass)), len = 10), include.lowest = TRUE)
color <- colorRampPalette(c("blue", "skyblue2", "lightgoldenrodyellow","lightsalmon","indianred4"))(9)[ii]
colfunc <- colorRampPalette(c("blue", "skyblue2", "lightgoldenrodyellow","lightsalmon","indianred4"))
color2 <- (colfunc(50))		  
windows(150,110)
par(fig=c(0.2,1,0,0.9), mar = c(2,4,2,2))
plot(c(-15,15), c(42,62), type='n',xlab="", ylab="", cex.axis=2, axes=TRUE)
points(dat$long, dat$lat,cex=2.3,pch=16, col=color)
map('worldHires', add=T, fill=T, col="black", xlim=c(-15,15), ylim=c(42,62))
par(fig=c(0,0.2,0,0.9), mar = c(2,0,0.5,2), new=TRUE)
legend_image <- grDevices::as.raster(matrix(rev(color2), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '',cex.main=0.5)
text(x=1.5, y = seq(0,1,l=3), labels = c(0.005,0.5,2.5), cex=1.2)
rasterImage(legend_image, 0, 0, 1,1)

### EVENNESS map
ii <- cut(log(dat$evesimpson), breaks = seq(min(log(dat$evesimpson)),max(log(dat$evesimpson)), len = 10), include.lowest = TRUE)
color <- colorRampPalette(c("blue", "skyblue2", "lightgoldenrodyellow","lightsalmon","indianred4"))(9)[ii]
colfunc <- colorRampPalette(c("blue", "skyblue2", "lightgoldenrodyellow","lightsalmon","indianred4"))
color2 <- (colfunc(50))		  
windows(150,110)
par(fig=c(0.2,1,0,0.9), mar = c(2,4,2,2))
plot(c(-15,15), c(42,62), type='n',xlab="", ylab="", cex.axis=2, axes=TRUE)
points(dat$long, dat$lat,cex=2.1,pch=16, col=color)
map('worldHires', add=T, fill=T, col="black", xlim=c(-15,15), ylim=c(42,62))
par(fig=c(0,0.2,0,0.9), mar = c(2,0,0.5,2), new=TRUE)
legend_image <- grDevices::as.raster(matrix(rev(color2), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '',cex.main=0.5)
text(x=1.5, y = seq(0,1,l=3), labels = c(0.03,0.2,0.35), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)

### SPECIES RICHNESS map
ii <- cut(dat$sr, breaks = seq(7, max(dat$sr), len = 10), include.lowest = TRUE)
color <- colorRampPalette(c("blue", "skyblue2", "lightgoldenrodyellow","lightsalmon","indianred4"))(9)[ii]
colfunc <- colorRampPalette(c("blue", "skyblue2", "lightgoldenrodyellow","lightsalmon","indianred4"))
color2 <- (colfunc(50))		  
windows(150,110)
par(fig=c(0.2,1,0,0.9), mar = c(2,4,2,2))
plot(c(-15,15), c(42,62), type='n',xlab="", ylab="", cex.axis=1.25, axes=TRUE)
points(dat$long, dat$lat,cex=2.2,pch=16, col=color)
map('worldHires', add=T, fill=T, col="black", xlim=c(-15,15), ylim=c(42,62))
par(fig=c(0,0.2,0,0.9), mar = c(2,0,0.5,2), new=TRUE)
legend_image <- grDevices::as.raster(matrix(rev(color2), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '',cex.main=0.5)
text(x=1.5, y = seq(0,1,l=3), labels = c(7,25,44), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)

### JACCARD map
ii <- cut(dat$jac, breaks = seq(0.5, 0.9, len = 10), include.lowest = TRUE)
color <- colorRampPalette(c("blue", "skyblue2", "lightgoldenrodyellow","lightsalmon","indianred4"))(9)[ii]
colfunc <- colorRampPalette(c("blue", "skyblue2", "lightgoldenrodyellow","lightsalmon","indianred4"))
color2 <- (colfunc(50))		  
windows(150,110)
par(fig=c(0.2,1,0,0.9), mar = c(2,4,2,2))
plot(c(-15,15), c(42,62), type='n',xlab="", ylab="", cex.axis=1.25, axes=TRUE)
points(dat$long, dat$lat,cex=2.2,pch=16, col=color)
map('worldHires', add=T, fill=T, col="black", xlim=c(-15,15), ylim=c(42,62))
par(fig=c(0,0.2,0,0.9), mar = c(2,0,0.5,2), new=TRUE)
legend_image <- grDevices::as.raster(matrix(rev(color2), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '',cex.main=0.5)
text(x=1.5, y = seq(0,1,l=3), labels = c(0.5,0.7,0.9), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)


######################################################################################
#### 2. PLOTS SCATTER PLOTS AND MODELS FROM FIGURE 2B,C,D,E
######################################################################################

dat <- read.csv(file='Dataset_S1.csv')

# Log values which are overdispersed
dat2 <- dat
dat2$biomass <- dat2$biomass/1e6

### Biomass and evenness
lm0 <- lm(log(dat2$biomass) ~ log(dat2$evesimpson))
windows(150,80)
par(mar=c(5,7,4,4), mfrow=c(2,4))
plot(log(dat2$biomass) ~ log(dat2$evesimpson), yaxt='n', xaxt='n', ylab='Biomass', xlab='Evenness', xlim=c(-3.32,-1.3), cex.lab=3, pch=20, cex=2, col='grey')
axis(side=2, at=sapply(c(0.005,0.05,0.5,2.5), FUN=function(x) log(x)),labels=c(0.005,0.05,0.5,2.5), cex.axis=2.5)
axis(side=1, at=sapply(c(0.01,0.1,0.2,0.3,0,4), FUN=function(x) log(x)),labels=c(0.01,0.1,0.2,0.3,0,4), cex.axis=2.5)
abline(lm0, col='red', lwd=5)

### Biomass and seasonality
lm0 <- lm(log(biomass) ~ sst.sea, data=dat2)
plot(log(dat2$biomass) ~ dat2$sst.sea, yaxt='n', xaxt='n', ylab='Biomass', xlab='SST.sea', cex.lab=3, pch=20, cex=2, col='grey')
axis(side=2, at=sapply(c(0.005,0.05,0.5,2.5), FUN=function(x) log(x)),labels=c(0.005,0.05,0.5,2.5), cex.axis=2.5)
axis(side=1, at=c(0.2,0.3,0.4), labels=c(0.2,0.3,0.4), cex.axis=2.5)
abline(lm0, col='red', lwd=5)

### Biomass and Jaccard
dat2$jac.c <- (dat2$jac-mean(dat2$jac))^2
dat2$jac2 <- dat2$jac^2
lm0 <- lm(log(biomass) ~ jac.c, data=dat2)
dat2<-dat2[order(dat2$jac),]
lm1 <- lm(log(biomass) ~ dat2$jac + dat2$jac2, data=dat2)
pred1 <- predict(lm1,list(var1=dat2$jac,var2=dat2$jac2))
plot(log(dat2$biomass) ~ dat2$jac, yaxt='n',  xaxt='n', ylab='Biomass', xlab='Dissimilarity', cex.lab=3, pch=20, cex=2, col='grey')
axis(side=2, at=sapply(c(0.005,0.05,0.5,2.5), FUN=function(x) log(x)),labels=c(0.005,0.05,0.5,2.5), cex.axis=2.5)
axis(side=1, at=c(0.6, 0.75, 0.9),labels=c(0.6, 0.75, 0.9), cex.axis=2.5)
lines(dat2$jac, pred1, col='black', lwd=5)

### Jaccard and Seasonality
lm0 <- lm(dat2$jac ~ dat2$sst.sea)
plot(dat2$jac ~ dat2$sst.sea, yaxt='n', xaxt='n', ylab='Dissimilarity', xlab='SST.sea', cex.lab=3, pch=20, cex=2, col='grey')
axis(side=2, at=c(0.6,0.75,0.9),labels=c(0.6,0.75,0.9), cex.axis=2.5)
axis(side=1, at=c(0.2,0.3,0.4),labels=c(0.2,0.3,0.4), cex.axis=2.5)
abline(lm0, col='black', lwd=5)	


#####################################################################################################
#### 3A. PIECEWISE SEM MODEL WITH BIOMASS, EVENNESS, RICHNESS AND JACCARD DISSIMILARITY WITH NEW DATA
#####################################################################################################

library(rgdal)
library(sp)
library(reshape2)
library(mgcv)
library(corrplot)
library(PerformanceAnalytics)
library(piecewiseSEM)
library(nlme)
library(gstat)
library(spdep)
library(ggplot2)
library(MuMIn)

dat <- read.csv(file='Dataset_S1.csv')

# Log-transformation
dat$biomass <- log(dat$biomass)
dat$chl <- log(dat$chl)
dat$depth <- log(dat$depth)
dat$evesimpson <- log(dat$evesimpson)

# Squared term
dat$jac.c <- (dat$jac-mean(dat$jac))^2

# Separate relationships - corresponds to results1
gls.b <- gls(biomass ~ evesimpson + jac.c + sst  + chl + pielou.sed + sst.sea + sr + jac, data=dat, correlation=corExp(form= ~ long+lat))
gls.j <- gls(jac ~ depth + sst.sea + chl, data=dat, correlation=corGaus(form= ~ long+lat))
gls.e <- gls(evesimpson ~ depth + pielou.sed, data=dat, correlation=corRatio(form= ~long+lat))
gls.sr <- gls(sr ~ sst.sea + pielou.sed + chl, data=dat, correlation=corExp(form= ~ long+lat))

# SEM with piecewiseSEM version 2.0
dat$effort <- NULL
modelList <- psem(
  gls.sr,
  gls.j,
  gls.e,
  gls.b,
  jac %~~% sr,
  evesimpson %~~% sr,
  jac %~~% evesimpson,
  jac.c %~~% jac,
  jac.c %~~% evesimpson,
  jac.c %~~% sr,
  dat
)

results <- summary(modelList)

rsq <- rsquared(modelList)
results$IC$AIC
results$dTable
results$Cstat
results$coefficients


#####################################################################################################
#### 3A. PIECEWISE SEM MODEL WITH FISHING EFFORT TESTED
#####################################################################################################

dat <- read.csv(file='Dataset_S1.csv')

# Remove missing values
dat <- subset(dat, !is.na(dat$effort))

# Log-transformation
dat$biomass <- log(dat$biomass)
dat$chl <- log(dat$chl)
dat$depth <- log(dat$depth)
dat$evesimpson <- log(dat$evesimpson)
dat$effort <- log(dat$effort)

# Squared term
dat$jac.c <- (dat$jac-mean(dat$jac))^2

# Separate models
gls.b <- gls(biomass ~ evesimpson + jac.c + sst + sst.sea + chl + pielou.sed + effort + sr + jac, data=dat, correlation=corExp(form= ~ long+lat))
gls.j <- gls(jac ~ depth + sst.sea + chl, data=dat, correlation=corGaus(form= ~ long+lat))
gls.e <- gls(evesimpson ~ depth + pielou.sed, data=dat, correlation=corRatio(form= ~long+lat))
gls.sr <- gls(sr ~ sst.sea + pielou.sed + chl, data=dat, correlation=corExp(form= ~ long+lat))

# SEM with piecewiseSEM 2.0
modelList <- psem(
  gls.sr,
  gls.j,
  gls.e,
  gls.b,
  jac %~~% sr,
  evesimpson %~~% jac,
  evesimpson %~~% sr,
  jac.c %~~% jac,
  jac.c %~~% evesimpson,
  jac.c %~~% sr,
  effort %~~% evesimpson,
  effort %~~% sr,
  effort %~~% jac,
  dat
)

results <- summary(modelList)


######################################################################################
#### 4. MODEL DIAGNOSTICS
######################################################################################

### gls models included in the SEM
gls.b <- gls(biomass ~ evesimpson + jac.c + sst + sst.sea + chl + pielou.sed, data=dat, correlation=corExp(form= ~ long+lat))
gls.j <- gls(jac ~ depth + sst.sea + chl, data=dat, correlation=corGaus(form= ~ long+lat))
gls.e <- gls(evesimpson ~ chl + depth + pielou.sed, data=dat, correlation=corRatio(form= ~long+lat))
gls.sr <- gls(sr ~ sst.sea + pielou.sed + chl, data=dat, correlation=corExp(form= ~ long+lat))

### Visual model diagnostics
windows(100,110)
par(mfrow=c(4,3), mar=c(4,5,1,2))
# RICHNESS RESIDUALS
reg <- gls.sr
rawresiduals <- resid(reg, 'pearson')
hist(resid(reg), xlim = c(-20,20), breaks=50, main='', xlab='Residuals', cex.lab=1.5, cex.axis=1.5)
plot( fitted(reg),resid(reg),
      col = "black", xlab = "Fitted class", ylab = "Residuals", cex.lab=1.5, cex.axis=1.5)
abline(h=0, lty='dashed')
qqnorm(resid(reg), main='' ,cex.lab=1.5, cex.axis=1.5)
qqline(resid(reg))
# EVENNESS RESIDUALS
reg <- gls.e
rawresiduals <- resid(reg, 'pearson')
hist(resid(reg), xlim = c(-2,2), breaks=50, main='', xlab='Residuals', cex.lab=1.5, cex.axis=1.5)
plot( fitted(reg),resid(reg),
      col = "black", xlab = "Fitted class", ylab = "Residuals", cex.lab=1.5, cex.axis=1.5)
abline(h=0, lty='dashed')
qqnorm(resid(reg), main='' ,cex.lab=1.5, cex.axis=1.5)
qqline(resid(reg))
# BETA-DIVERSITY RESIDUALS
reg <- gls.j
rawresiduals <- resid(reg, 'pearson')
hist(resid(reg), xlim = c(-0.2,0.2), breaks=50, main='', xlab='Residuals', cex.lab=1.5, cex.axis=1.5)
plot( fitted(reg),resid(reg),
      col = "black", xlab = "Fitted class", ylab = "Residuals", cex.lab=1.5, cex.axis=1.5)
abline(h=0, lty='dashed')
qqnorm(resid(reg), main='' ,cex.lab=1.5, cex.axis=1.5)
qqline(resid(reg))
# BIOMASS RESIDUALS
reg <- gls.b
rawresiduals <- resid(reg, 'pearson')
hist(resid(reg), xlim = c(-2,2), breaks=50, main='', xlab='Residuals', cex.lab=1.5, cex.axis=1.5)
plot( fitted(reg),resid(reg),
      col = "black", xlab = "Fitted class", ylab = "Residuals", cex.lab=1.5, cex.axis=1.5)
abline(h=0, lty='dashed')
qqnorm(resid(reg), main='' ,cex.lab=1.5, cex.axis=1.5)
qqline(resid(reg))




######################################################################################
#### 5. PLOT FROM FIGURE 3A,B,C
######################################################################################

# FIGURE 3A
dat <- read.csv(file='Dataset_S1.csv')
dat$biomass <- dat$biomass/1e6
ii <- cut(dat$jac, breaks = seq(min(dat$jac), max(dat$jac), len = 10), include.lowest = TRUE)
color <- colorRampPalette(c("blue", "skyblue2", "lightgoldenrodyellow","lightsalmon","indianred4"))(9)[ii]
windows(80,80)
par(mar=c(6,6,4,4))
plot(log(dat$biomass) ~ log(dat$evesimpson),col=color, yaxt='n', xaxt='n', ylab='Biomass', xlab='Evenness', xlim=c(-3.32,-1.3), cex.lab=3, pch=20, cex=2)
axis(side=2, at=sapply(c(0.005,0.05,0.5,2.5), FUN=function(x) log(x)),labels=c(0.005,0.05,0.5,2.5), cex.axis=2.5)
axis(side=1, at=sapply(c(0.01,0.05,0.1,0.2,0.3), FUN=function(x) log(x)),labels=c(0.01,0.05,0.1,0.2,0.3), cex.axis=2.5)

# FIGURE3B,C
lm0 <- lm(log(dat$biomass) ~ log(dat$evesimpson))
res <- resid(lm0)

windows(100,80)
par(mar=c(4,6,4,4))
plot(res ~ dat$feeding.mode_benthivorous, ylab='Residuals',xlab='CWM Benthivorous', pch=20, col='grey', cex=1.75, cex.lab=2.5, cex.axis=1.5, xaxt='n', yaxt='n')
abline(lm(res ~ dat$feeding.mode_benthivorous), lwd=8)
axis(side=1, at=c(0.2,0.4,0.6), cex.axis=2.5)
axis(side=2, at=c(-2,0,2), cex.axis=2.5)

windows(100,80)
par(mar=c(4,6,4,4))
plot(res ~ dat$tl, ylab='Residuals', xlab='CWM TL', pch=20, col='grey',xlim=c(3.5,4.0), cex=1.75, cex.lab=2.5, cex.axis=2, xaxt='n', yaxt='n')
abline(lm(res ~ dat$tl), lwd=8)
axis(side=1, at=c(3.5, 3.75, 4.0), cex.axis=2.5)
axis(side=2, at=c(-2,0,2), cex.axis=2.5)
