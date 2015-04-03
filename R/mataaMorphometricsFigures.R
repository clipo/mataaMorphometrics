## Figures for Mata'a Morphometrics paper
setwd("/Users/clipo/mataaMorphometrics/R") ## SET THE WORKING DIRECTORY
library(knitr)
library(httr)
library(Rmisc)
library(lattice)
library(plyr)
library(plotrix)
library(rgl)
#set_config( config( ssl.verifypeer = 0L ) )
opts_chunk$set(fig.align='left') 
options(repos = c(CRAN = "http://cran.rstudio.com"))
set_config( config( ssl.verifypeer = 0L ) )
library(devtools)
install_github("vbonhomme/Momocs")
library(Momocs)

# some constants we use in the file
data.path <- "../Data/xyfinaldatawithids.txt" 
scalingFactor = 75 ## this is the # of pixels per cm. 
xCenter = 12 # x coordinate of centroid (if its fixed, otherwise can be aligned based on average centroid, etc.)
yCenter = 11 # y coordinate of centroid

# We read the whole file
allPAST <- read.table(data.path, header = TRUE, sep="\t")
PAST <-allPAST[which(allPAST$Island=="Rapa_Nui"),]
nonRapaNui <-PAST[which(allPAST$Island != "Rapa_Nui"),]

last.meta <- which(colnames(PAST)=="X1") - 1
fac <- PAST[, 1:last.meta]
allfac<-allPAST[,1:last.meta]
nonRapaNuiFac <-nonRapaNui[,1:last.meta]

# xy will contain coordinates only
xy <- as.matrix(PAST[, -c(1:last.meta)])
allxy<- as.matrix(allPAST[,-c(1:last.meta)])
nonRapaNuiXY <- as.matrix(nonRapaNui[,-c(1:last.meta)])

# scale the coordinates so that pixels = cm
scaledXY <- xy/scalingFactor
allscaledXY <-allxy/scalingFactor
nonRapaNuiScaledXY <- nonRapaNuiXY/scalingFactor

# a short loop to reorder thing and store them in a list
coo <- list()
allcoo <-list()
nonRapaNuicoo <-list()
for (i in 1:nrow(scaledXY)){ 
  coo[[i]] <- cbind(scaledXY[i, seq(1, ncol(scaledXY), 2)], scaledXY[i, seq(2, ncol(scaledXY), 2)])
}

for (i in 1:nrow(allscaledXY)){ 
  allcoo[[i]] <- cbind(allscaledXY[i, seq(1, ncol(allscaledXY), 2)], allscaledXY[i, seq(2, ncol(allscaledXY), 2)])
}

for (i in 1:nrow(nonRapaNuiScaledXY)){ 
  nonRapaNuicoo[[i]] <- cbind(nonRapaNuiScaledXY[i, seq(1, ncol(nonRapaNuiScaledXY), 2)], nonRapaNuiScaledXY[i, seq(2, ncol(nonRapaNuiScaledXY), 2)])
}
# we renames the components of the list (ie the shapes)
names(coo) <- fac[, 1]
#names(allcoo) <-allfac[,1]

# now we create the Out object (formerly Coo,
# but Coo is now a super-class now in order to handle outlines, open outlines, and
# landmarks)
RapaNui <- Out(coo, fac=fac)

## include the non Rapa Nui examples as well
allMataa <-Out(allcoo, fac=allfac)

## Just the non-Rapa Nui examples 
allMataaNonRapaNui <-Out(nonRapaNuicoo, fac=allfac)

## Descriptive information for the paper
numberOfMataa <- nrow(scaledXY)
numberOfNonRapaNuiMataa <- nrow(nonRapaNuiScaledXY)

# Create vectors of the lengths and widths for descr. stats. 
## note that in the images 1 cm = 75 pixels
lengths <- 0
widths <- 0
for (i in 1:nrow(scaledXY)){ 
  lengths[i] <- (max(scaledXY[i, seq(2, ncol(scaledXY), 2)])-min(scaledXY[i, seq(2, ncol(scaledXY), 2)]))
  widths[i] <- (max(scaledXY[i, seq(1, ncol(scaledXY), 2)])-min(scaledXY[i, seq(1, ncol(scaledXY), 2)]))
}

# descriptive values
meanLength <-mean(lengths)
sdLength <-sd(lengths)
meanWidth <- mean(widths)
sdWidth <-sd(widths)

# Determine maximum width and length
xMaxLength<-0
xMaxWidth<-0
for (i in 1:nrow(scaledXY)){ 
  xMaxLength[i]<-(max(scaledXY[i, seq(2, ncol(scaledXY), 2)])-min(scaledXY[i, seq(2, ncol(scaledXY), 2)]) )
  xMaxWidth[i]<-(max(scaledXY[i, seq(1, ncol(scaledXY), 2)])-min(scaledXY[i, seq(1, ncol(scaledXY), 2)]) )
}

# Create some descriptions of the numbers of mataa in different categories
englertMataa <- length(which(fac[3] == "Englert"))
bishopMataa <- length(which(fac[3]=="Bishop"))
unknownMataa <- length(which(fac[3]=="Unknown"))
numberOfSites <-length(fac[4])-1

# the dmp function to create the polar plots
dmp <- function(q, xlab="", ylab="", title="", th0=pi/2,
                cols, palette=col.summer, leg, ...){
  op <- par(mar=c(5.1, 0, 4.1, 0))
  q <- q/max(q)
  nr <- nrow(q)
  nc <- ncol(q)
  th <- seq(0, 2*pi, length=nc) + th0
  xi <- t(apply(q, 1, function(r) r*cos(th)))
  yi <- t(apply(q, 1, function(r) r*sin(th)))
  if (missing(cols)) cols <- palette(nr)
  plot(NA, asp=1, xlim=c(-1, 1), ylim=c(-1, 1),
       xlab=xlab, ylab=ylab, main=title,
       xaxs="i", las=1, ann=FALSE, axes=FALSE, frame=FALSE)#, ...)
  lines(cos(th), sin(th), col="grey80")
  for (i in seq(0, 1, 0.25)){
    lines(i*cos(th), i*sin(th), col="grey80", lty=2)}
  th.grid <- seq(0, 2*pi, length=9)
  segments(0, 0, cos(th.grid), sin(th.grid), col=c("grey80"), lty=2)
  for (i in 1:nr) lines(xi[i, ], yi[i, ], col=cols[i])
  if (missing(leg)) leg <- rownames(q)
  legend("topright", lwd=1, col=cols, legend=leg, cex=3/4, bty="n")}

## plot the confidence intervals
ci.plot <- function(x){
  x <- lapply(x$coo,coo.interpolate,360)
  x <- matrix(unlist(x), nrow=length(x), byrow=TRUE)
  xCI <- apply(x, 2, CI,ci = 0.99)
  ## now add some and substract to differentiate between the 3 lines.. this is just arbitrary spacing.
  l<-xCI[1,]-1
  m<-xCI[2,]
  u<-xCI[3,]+1
  xxCI=rbind(l,m,u)
  dmp(xxCI)
}

quartiles.plot<-function(x){
  x <-lapply(x$coo,coo.interpolate,360)
  x <- matrix(unlist(x), nrow=length(x), byrow=TRUE)
  quartilesCI <- apply(x, 2, quantile, probs=seq(0, 1, 0.25))
  dmp(quartilesCI)
}

#Table 1
byCollection <- table(fac$Collection,fac$Site)
bySource <- table(fac$Source,fac$Site)
knitr:::kable(byCollection, digits=0)

#Table 2:  *Mata'a* included in analyses by obsidian source.
knitr:::kable(bySource,digits=0)

##Figure 1. Pacific island locations mentioned in the text.
##Figure 2. Examples of mata’a from Rapa Nui. These mata’a are part of the collections at the P. Sebastian Englert Museum, Hanga Roa, Isla de Pascua. 
##Figure 3. Locations of mata’a collections and obsidian sources on Rapa Nui, Chile.
##Figure 4. Measurement process used to generate outline coordinates for each mata’a in the study. (A) First, we took a scaled digital photo is taken of the object. We also ensure that all images are resized so that they are equivalent in scale. (B) Second, we isolated the outline of mata’a from the background using TPS Dig software (Rohlf 2014). (C) Along the perimeter of the artifact, we placed 400 semilandmarks and recorded their x- and y-coordinates. We then use these coordinate data in the morphometric analyses.

#Figure 5. *Mata'a* length and width.
par(mfrow=c(1, 2))
hist(lengths, main="Length (cm)" ,  xlab="length (cm)", breaks = 8) 
hist(widths, main="Width (cm)", xlab="width (cm)", breaks = 8) 
par(mfrow=c(1, 1))

#Figure 6. Superimposed *mata'a* outlines from Rapa Nui.  For comparison, all *mata'a* are aligned at the center point of the haft where it meets the blade.
stack(RapaNui)

RapaF <- efourier(RapaNui, 12, norm=TRUE)  #norm makes sure they are normalized by total length of perimeter

##Figure 7. For all positive integers, the sum of a cosine curve and a sine curve defines an ellipse in the plane. Elliptic Fourier analysis is based on an harmonic sum of such ellipses. Five harmonics are here shown at four locations on the original outline of a *mata'a*. As the number of harmonics is increased the reconstruction better approximates the original shape outline.
Ptolemy(RapaNui[1], nb.h = 5, t = seq(0, 2*pi, pi/2), legend = TRUE)


#Figure 8.  First two principal components (PC1 and PC2 are on the x- and y-axis, respectively) for the Rapa Nui *mata'a*.
plot(PCA(RapaF),  ellipses=TRUE) 

#Figure 9.  First two principal components of *mata'a* grouped by site location. 
plot(PCA(RapaF), "Site",  ellipses=TRUE, conf_ellipses=0.90, conf_ellipsesax = c(0.2, 0.4, 0.75, 0.8))
# regular EFT with normalized coefficients
manovaSite<-MANOVA_PW(PCA(RapaF), "Site")
knitr:::kable(manovaSite$summary, digits=4)

#Figure 10.  First two principal components of *mata'a* grouped by obsidian source. 
plot(PCA(RapaF), "Source",  ellipses=TRUE,  conf_ellipses=0.90, conf_ellipsesax = c(0.2, 0.4, 0.75, 0.8)) # regular EFT with normalized coefficients
manovaObsidian<-MANOVA_PW(PCA(RapaF), "Source")
# p-values for table in figure
knitr:::kable(manovaObsidian$summary, digits=4)

##Figure 11:  Factorial maps depicting the two principal compoents (PC1 and PC2 are the x- and y-axis, respectively) of morphological variation for stemmed lithic shaped objects from Rapa Nui, New Britain, New Zealand, Chatham and Pitcairn Islands. The shapes are reconstructed from the factorial map using the first two component axes.
allF <- efourier(allMataa, 13, norm=TRUE)  ## note that this is scaled

allP <- PCA(allF)
manovaIsland<-MANOVA_PW(filter(allP, !Island %in% c("New_Zealand", "Pitcairn")), "Island")
plot(PCA(allF), "Island", ellipses=TRUE,  conf_ellipses=0.90, conf_ellipsesax = c(0.2, 0.4, 0.75, 0.8))
knitr:::kable(manovaIsland$summary, digits=4)

## p-values for smaller island comparisons
New_Zealand_scores <- filter(allP, Island=="New_Zealand")$x[, 1]
New_Britain_scores <- filter(allP, Island=="New_Britain")$x[, 1]
Chatham_scores <- filter(allP, Island=="Chatham")$x[, 1]
Pitcairn_scores <- filter(allP, Island=="Pitcairn")$x[, 1]
Rapa_Nui_scores <- filter(allP, Island=="Rapa_Nui")$x[, 1]

wilcox.test(New_Zealand_scores, New_Britain_scores)
wilcox.test(New_Zealand_scores, Chatham_scores)
wilcox.test(New_Zealand_scores, Pitcairn_scores)
wilcox.test(New_Zealand_scores, Rapa_Nui_scores)

wilcox.test(Chatham_scores, New_Britain_scores)
wilcox.test(Chatham_scores, New_Zealand_scores)
wilcox.test(Chatham_scores, Pitcairn_scores)
wilcox.test(Chatham_scores, Rapa_Nui_scores)

########################## SUPPLEMENTAL FIGURES AND TABLES ###########################

###Figure S1. Sample size and *mata'a* parameter estimation. The sample size (N=423) of 
###mata'a* appears to be sufficient to estimate variability in the basic shapes.
###

par(mfrow=c(1, 2))
lengthFrame <- data.frame(number=rep(NA,(numberOfMataa*2)), mean=NA, sd=NA,LL=NA, UL=NA)
widthFrame <- data.frame(number=rep(NA,(numberOfMataa*2)), mean=NA, sd=NA,LL=NA, UL=NA)
lengthMean = NULL
widthMean = NULL
lengthSD = NULL
numberOfSamples=100
for (i in 1:(2*numberOfMataa) ) {
  for (n in 1:numberOfSamples) {
    subsampleL <- sample(lengths,size=i,replace=TRUE)
    subsampleW <- sample(widths,size=i,replace=TRUE)
    meanL <- mean(subsampleL)
    sdL <- sd(subsampleL)
    meanW <-var(subsampleW)
    lengthMean[n] = meanL
    lengthSD[n] = sdL
    widthMean[n] = meanW
  }
  lengthFrame[i, ] <-c( i, mean(lengthMean), sd(lengthMean),mean(lengthMean)-2*sd(lengthMean), mean(lengthMean)+2*sd(lengthMean))
  widthFrame[i, ] <-c( i, mean(widthMean), sd(widthMean),mean(widthMean)-2*sd(widthMean), mean(widthMean)+2*sd(widthMean))
}
#plotCI(lengthFrame$number,lengthFrame$mean, ui=lengthFrame$UL, li=lengthFrame$LL )
plot(lengthFrame$number, lengthFrame$mean,  type = "l", ylim=c(4,12), xlab="Sample Size (N)", ylab="Length (cm)")
polygon(c(lengthFrame$number,rev(lengthFrame$number)),c(lengthFrame$LL,rev(lengthFrame$UL)),col = "grey75", border = FALSE)
lines(lengthFrame$number,lengthFrame$mean, col="red",lty=2)
abline(v=numberOfMataa, lty="dotted")

#plotCI(widthFrame$number,widthFrame$mean, ui=widthFrame$UL, li=widthFrame$LL )
plot(widthFrame$number, widthFrame$mean,  type = "l",ylim=c(0,12), xlab="Sample Size (N)", ylab="Width (cm)")
polygon(c(widthFrame$number,rev(widthFrame$number)),c(widthFrame$LL,rev(widthFrame$UL)),col = "grey75", border = FALSE)
lines(widthFrame$number,widthFrame$mean, col="red",lty=2)
abline(v=numberOfMataa, lty="dotted")

par(mfrow=c(1, 1))

## Figure S2. Mata’a included in the current analyses. 
##The 5 colors indicate the collection locations on Rapa Nui 
#(Blue=Ahu Tautiri, Green=Orito, Yellow/Green=Orongo, Orange=Rano Kao, 
#Red=Location only known to the level of the island, Yellow=Parcela).
Momocs:::panel(RapaNui, fac = "Site")

###Figure S3. *Mata'a* reconstructed from different numbers of 
### harmonics. Twelve harmonics provide a satisfactory reconstruction.

## Note: Im not sure what function this should be in 0.99 
hcontrib(RapaNui, method = "efourier", id = 16, title="*Mata'a* reconstructions with harmonics", harm.range = 1:36,  palette = col.sari, plot.method = "panel")

###Figure S4.  Cumulated harmonic Fourier power calculated from Rapa Nui
###*mata'a*. The 12 first harmonics gather nearly 100% of the harmonic power. 
###Maxima, minimaRapa and medians are also plotted.
harm_pow(RapaNui, method = "efourier", title="Fourier coefficients power spectrum")

##########################EXTRA 3D STUFF ##########################


library(rgl)
## playing with 3D
bp <- PCA(efourier(RapaNui, 12, norm=TRUE))
x <- bp$x[, 1]
y <- bp$x[, 2]
z <- bp$x[, 3]
f <- bp$fac$Site
cg <- col.qual(8)
cp <- cg[f]
level=0.5
open3d()
plot3d(x, y, z, col=cp, xlab = "PC1", ylab="PC2", zlab="PC3", size=3, aspect=TRUE)
for (i in seq(along=levels(f))){
  m <- cbind(x, y, z)[f==levels(f)[i],]
  cm <- apply(m, 2, mean)
  #points3d(cm[1], cm[2], cm[3], col=cg[i], size=5)
  apply(m, 1, function(x)
    segments3d(c(cm[1], x[1]),
               c(cm[2], x[2]),
               c(cm[3], x[3]),
               col=cg[i], alpha=0.5))
  
  shade3d(ellipse3d(cov(m), col=cg[i], level=level,  alpha=0.1))
  wire3d(ellipse3d(cov(m), col=cg[i], level=level, alpha=0.05))}

xx <- spin3d(axis = c(1, 1, 1), rpm = 10)
play3d(xx, duration = 6)


