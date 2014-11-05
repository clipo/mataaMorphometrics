## Figures for Mata'a Morphometrics paper
setwd("/Volumes/Macintosh HD/Users/carllipo/mataaMorphometrics/R") ## SET THE WORKING DIRECTORY
library(knitr)
library(httr)
library(Rmisc)
library(lattice)
library(plyr)
library(plotrix)
#set_config( config( ssl.verifypeer = 0L ) )
opts_chunk$set(fig.align='left') 
options(repos = c(CRAN = "http://cran.rstudio.com"))
set_config( config( ssl.verifypeer = 0L ) )
library(devtools)
RCurl::getBinaryURL(url='https://github.com/vbonhomme/Momocs.zip', cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl"), followlocation=1L)
library(Momocs)

# some constants we use in the file
data.path <- "../Data/xyfinaldatawithids.txt" 
scalingFactor = 75 ## this is the # of pixels per cm. 
xCenter = 12 # x coordinate of centroid (if its fixed, otherwise can be aligned based on average centroid, etc.)
yCenter = 11 # y coordinate of centroid

# We read the whole file
allPAST <- read.table(data.path, header = TRUE, sep="\t")
PAST <-allPAST[which(allPAST$Island=="Rapa_Nui"),]
nonRapaNui <-PAST[which(allPAST$Island=="New_Britain"),]

last.meta <- which(colnames(PAST)=="X1") - 1
fac <- PAST[, 1:last.meta]
allfac<-allPAST[,1:last.meta]

# xy will contain coordinates only
xy <- as.matrix(PAST[, -c(1:last.meta)])
allxy<- as.matrix(allPAST[,-c(1:last.meta)])
scaledXY <- xy/scalingFactor
allscaledXY <-allxy/scalingFactor

# a short loop to reorder thing and store them in a list
coo <- list()
allcoo <-list()
for (i in 1:nrow(scaledXY)){ 
  coo[[i]] <- cbind(scaledXY[i, seq(1, ncol(scaledXY), 2)], scaledXY[i, seq(2, ncol(scaledXY), 2)])
}

for (i in 1:nrow(allscaledXY)){ 
  allcoo[[i]] <- cbind(allscaledXY[i, seq(1, ncol(allscaledXY), 2)], allscaledXY[i, seq(2, ncol(allscaledXY), 2)])
}
# we renames the components of the list (ie the shapes)
names(coo) <- fac[, 1]
#names(allcoo) <-allfac[,1]

# now we create the Out object (formerly Coo,
# but Coo is now a super-class now in order to handle outlines, open outlines, and
# landmarks)
RapaNui <- Out(coo, fac=fac)

## include the non Rapa Nui examples
allMataa <-Out(allcoo, fac=allfac)

## Descriptive information for the paper
numberOfMataa <- nrow(scaledXY)

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


# the dmp function
dmp <- function(q, xlab="", ylab="", title="", th0=pi/2,
                cols, palette=col.summer, leg, ...){
  op <- par(mar=c(5.1, 0, 4.1, 0))
  q <- q/max(q)
  nr <- nrow(q)
  nc <- ncol(q)
  #th <- seq(0, 2*pi, length=nc+1)[-(nc+1)] + th0
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
  #x <- lapply(x$coo, function(x) apply(x, 1, ed, pt2=c(-xCenter, -yCenter)))
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
  #x$coo <- lapply(coo.trans$coo, x$coo, - 625, - 500)
  #x <- lapply(x$coo, function(x) apply(x, 1, ed, pt2=c(-xCenter, -yCenter)))
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

##Figure 2. Locations of *mata'a* collections from Rapa Nui, Chile.
Momocs:::panel(RapaNui, fac = "Site")

##Figure 4. Measurement process used to generate semilandmark data for each *mata'a* image.

#Figure 5. *Mata'a* length and width.
par(mfrow=c(1, 2))
hist(lengths, main="Length (cm)" ,  xlab="length (cm)", breaks = 8) 
hist(widths, main="Width (cm)", xlab="width (cm)", breaks = 8) 
par(mfrow=c(1, 1))

#Figure 6. Superimposed *mata'a* outlines from Rapa Nui.  For comparison, all *mata'a* are aligned at the center point of the haft where it meets the blade.
stack(RapaNui)

RapaF <- eFourier(RapaNui, 12, norm=TRUE)  #norm makes sure they are scaled

##Figure 8. For all positive integers, the sum of a cosine curve and a sine curve defines an ellipse in the plane. Elliptic Fourier analysis is based on an harmonic sum of such ellipses. Five harmonics are here shown at four locations on the original outline of a *mata'a*. As the number of harmonics is increased the reconstruction better approximates the original shape outline.
Ptolemy(RapaNui[1], nb.h = 5, t = seq(0, 2*pi, pi/2), legend = TRUE)
##Figure 9. *Mata'a* reconstructed from different numbers of harmonics. Twelve harmonics provide a satisfactory reconstruction.
hqual(RapaNui, method = "eFourier", id = 16, title="*Mata'a* reconstructions with harmonics", harm.range = 1:36,  palette = col.sari, plot.method = "panel")

##Figure 10.  Cumulated harmonic Fourier power calculated from Rapa Nui *mata'a*. The 12 first harmonics gather nearly 100% of the harmonic power. Maxima, minimaRapa and medians are also plotted.
hpow(RapaNui, method = "efourier", title="Fourier coefficients power spectrum")

#Figure 11.  First two principal components (PC1 and PC2 are on the x- and y-axis, respectively) for the Rapa Nui *mata'a*.
plot(PCA(RapaF),  ellipses=TRUE) 

#Figure 12.  First two principal components of *mata'a* grouped by site location. 
plot(PCA(RapaF), "Site",  ellipses=TRUE) # regular EFT with normalized coefficients

#Figure 13.  First two principal components of *mata'a* grouped by obsidian source. 
plot(PCA(RapaF), "Source",  ellipses=TRUE) # regular EFT with normalized coefficients

allF <- eFourier(allMataa, 13, norm=TRUE)  ## note that this is scaled

##Figure 14:  Factorial maps depicting the two principal compoents (PC1 and PC2 are the x- and y-axis, respectively) of morphological variation for stemmed lithic shaped objects from Rapa Nui, New Britain, New Zealand, Chatham and Pitcairn Islands. The shapes are reconstructed from the factorial map using the first two component axes.
plot(PCA(allF), "Island", ellipses=TRUE)

##Figure S1. Sample size and *mata'a* parameter estimation. The sample size (N=417) of *mata'a* appears to be sufficient to estimate variability in the basic shapes.
par(mfrow=c(1, 2))
lengthFrame <- data.frame(number=rep(NA,(numberOfMataa*2)), mean=NA, sd=NA,LL=NA, UL=NA)
widthFrame <- data.frame(number=rep(NA,(numberOfMataa*2)), mean=NA, sd=NA,LL=NA, UL=NA)
lengthMean = NULL
widthMean = NULL
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
#plot(lengthMean, ylab="Lengths (cm)",xlab="Sample Size (N)")

#plot(widthMean, ylab="Widths (cm)", xlab="Sample Size (N)")
#abline(v=numberOfMataa, lty="dotted")
par(mfrow=c(1, 1))

##########################

PAST <- subset(allPAST, Island="Rapa_Nui" & Site != "Unknown")
last.meta <- which(colnames(PAST)=="X1") - 1
fac <- PAST[, 1:last.meta]


# xy will contain coordinates only
xy <- as.matrix(PAST[, -c(1:last.meta)])
scaledXY <- xy/scalingFactor

# a short loop to reorder thing and store them in a list
coo <- list()
for (i in 1:nrow(scaledXY)){ 
  coo[[i]] <- cbind(scaledXY[i, seq(1, ncol(scaledXY), 2)], scaledXY[i, seq(2, ncol(scaledXY), 2)])
}

# we renames the components of the list (ie the shapes)
names(coo) <- fac[, 1]

# now we create the Out object (formerly Coo,
# but Coo is now a super-class now in order to handle outlines, open outlines, and
# landmarks)
RapaNui <- Out(coo, fac=fac)


## Descriptive information for the paper
numberOfMataa <- nrow(scaledXY)
library(rgl)
## playing with 3D
bp <- PCA(eFourier(RapaNui, 12, norm=TRUE))
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


