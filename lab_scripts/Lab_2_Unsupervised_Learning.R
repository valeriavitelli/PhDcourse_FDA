
rm(list=objects())
setwd("C:/Users/valeriv/UiO Dropbox/Valeria Vitelli/Teaching/FDA_Palermo/R scripts")
set.seed(25032024)

library(refund)
library(mgcv)
library(fields)
library(fda)
library(funFEM)
library(fdacluster)

###############################################
###               PART 1 -- AIM:            ###
### FUNCTIONAL PRINCIPAL COMPONENT ANALYSIS ###
###############################################

### Illustration of fPCA
### DATA: diffusion tensor imaging (DTI) 

data(DTI); attach(DTI)
names(DTI)

DTI.complete <- subset(DTI, complete.cases(DTI))
DTI.baseline <- subset(DTI.complete, visit == 1 & case == 1)
p <- dim(DTI.baseline$cca)[2]
tract <- 1:p
n <- length(unique(DTI.baseline$ID))
dim(DTI.baseline$cca) # sanity check: dimensions

# Recall from Lab 1:
# we focus on Fractional Anisotropy (FA) along the corpus callosum (CCA) tract
# collected from multiple sclerosis (MS) patients without any missing values

# plot the sample and 3 randomly selected curves
x11()
matplot(tract, t(DTI.baseline$cca), 
        type='l', lty=1, col="light grey",
        main = "Diffusion Tensor Imaging : CCA",
        xlab="tract", ylab="Fractional anisotropy (FA)")
sel.crv <- sample(1:n, size = 3, replace = FALSE)
matlines(tract, t(DTI.baseline$cca[sel.crv,]), 
         type='l', lty=1, lwd=2, col = rainbow(3))

## STEP 1.1: Smooth each curve and take pointwise average
#########################################################
smooth.curves <- array(0, dim=c(n, p))
 
# alternative approach to smoothing: gam function in the mgcv package
# various basis functions, selection methods for smoothing parameter, flexible
# look into one case closely: for more info see ?mgcv::gam -> smooth.terms
j = 1
fit <- gam(DTI.baseline$cca[j,] ~ s(tract, k = 10, bs = 'cr'), method = "REML") # 'cr' = Cubic regression splines
plot(tract, DTI.baseline$cca[j,], ylab="Fractional anisotropy (FA)",
     main='Smoothing performance, one curve')
lines(tract, fit$fitted)

# seems ok! let's do this for all curves
for(j in 1:n){
  fit <- gam(DTI.baseline$cca[j,] ~ s(tract, k = 10, bs = 'cr'), method = "REML")
  smooth.curves[j,] <- fit$fitted
}

# compare data and smoothed from the 3 randomly selected curves we looked at before
matplot(tract, t(DTI.baseline$cca[sel.crv,]), 
        type='p', pch=1, lwd=1, col = rainbow(3),
        ylab="Fractional anisotropy (FA)",
        main='Smoothing performance, 3 random curves')
matlines(tract, t(smooth.curves[sel.crv,]), 
         type='l', lty=1, lwd=1, col = rainbow(3))

# now compute pointwise average
mean.hat <- colMeans(smooth.curves)
# and plot it
matplot(tract, t(smooth.curves), 
        type='l', lty=1, col="light grey",
        main = "DTI data: smooth curves and smooth mean function",
        xlab="tract", ylab="Fractional anisotropy (FA)")
lines(tract, mean.hat, col='blue', lwd=2)


# Note that this is slightly different from taking the average of un-smoothed 
# raw functional data! Some noise has been removed (pre-smoothing approach...)
matplot(tract, t(DTI.baseline$cca), 
        type='l', lty=1, col="light grey",
        main = "DTI data: original curves and point-wise raw mean function",
        xlab="tract", ylab="Fractional anisotropy (FA)")
lines(tract, colMeans(DTI.baseline$cca), col='blue', lwd=2)


## STEP 1.2: fPCA -- Spectral decomposition of the estimated covariance
#######################################################################

# point-wise covariance function
smooth.cov <- cov(smooth.curves)
image.plot(tract, tract, smooth.cov, 
           main='Smooth covariance of FA')

# Estimated eigenfunctions and eigenvalues of the covariance function
# can be obtained from the spectral decomposition of the estimated
# point-wise covariance function

svd.result0 <- eigen(smooth.cov, symmetric = TRUE)
names(svd.result0)
evectors <- svd.result0$vectors[,svd.result0$values > 0]
evalues <- svd.result0$values[svd.result0$values > 0]
head(colSums(evectors^2)) # returns unitary vectors 

# eigen function returns unitary vectors -> scale them by \sqrt{p}
efns0 <- evectors*sqrt(p)
evals0 <- evalues/p
pve <- cumsum(evals0)/sum(evals0)
npc <- sum(pve < 0.95) + 1 # select the PCs that give at least 95% explained variance

# truncated estimated eigen components
efns <- efns0[,1:npc]
evals <- evals0[1:npc]

# scree plot
x11()
plot(1:20, pve[1:20], pch = 16, 
     ylab="percentage of variance explained", xlab="number of PCs",
     main="DTI data: scree plot")
abline(h = 0.95, lty=2, col='red')


## STEP 1.3 -- Interpretation: the first 6 PCs explain together > 95% variance
##############################################################################

# let's look at them
x11()
matplot(tract, efns[,1:6], col=rainbow(6), 
        type='l', lty=1, lwd=2,
        ylab="eigenfunctions", xlab="tract",
        main="First 6 eigenfunctions")
legend('topleft',col = rainbow(6), legend = paste('eigenfunction ',1:6, sep = ''), 
       bty = 'n', lwd = 2, cex = .6)

## Useful visualization (remember the lecture!):
## visualize the effect of each PC on the mean! 
## HOW? For the i-th component, plot the mean \mu(t) as a line, 
## and then '+' for \mu(t) + 2*\sqrt{\lambda_i} \phi_i(t)
##          '-' for \mu(t) - 2*\sqrt{\lambda_i} \phi_i(t)
k.pc <- 1
effect <- efns[, k.pc]*2*sqrt(evals[k.pc])
mat <- cbind(mean.hat - effect,  mean.hat + effect)

x11()
par(mfrow=c(2,1))
plot(tract, efns[,k.pc], lty=1, lwd=2, type='l', ylim=c(-2,2),
     main = paste0("fPC",k.pc), ylab="", xlab="tract" )
abline(h = 0, lty=3)

matplot(tract, mat, type='p', col=c(2,4), pch = c("-", "+"),
        ylab="", xlab="tract", 
        main = paste0("fPC",k.pc, " (",round(pve[k.pc]*100,1),"%)"))
lines(tract, mean.hat, lty=1, lwd=1)
# (repeat for all fPCs and try to interpret)

## Estimation and Visualization of scores & fitted curves (using the first 6 fPCs)
demeaned <- DTI.baseline$cca - t(matrix(rep(mean.hat, n),
                                        nrow=length(mean.hat)))
scores <- matrix(NA, nrow=n, ncol=npc)
fitted <- array(NA, dim(DTI.baseline$cca))
for(i in 1:n){
  scores[i,] <- colMeans(matrix(rep(demeaned[i,], npc), nrow=p) * efns)
  fitted[i,] <- mean.hat + scores[i,]%*%t(efns)
}

matplot(tract, t(DTI.baseline$cca[sel.crv,]), pch = "o", cex = 0.5,
        ylab="", xlab="tract", main="Fitted curves via fPCA", col = rainbow(3))
matlines(tract, t(fitted[sel.crv,]), type='l', lwd=2, lty=1, col = rainbow(3))


##################################
###       PART 2 -- AIM:       ###
### CLUSTERING FUNCTIONAL DATA ###
##################################

# Recall from Lab 1:
# we focus on the Berkeley Growth Study data
# heights of 54 girls and 39 boys measured at a set of 31 ages

?growth
names(growth)

data <- t(cbind(growth$hgtm,growth$hgtf))
n <- dim(data)[1]
x <- growth$age
gender <- c(rep(1, dim(growth$hgtm)[2]), rep(2, dim(growth$hgtf)[2]))

# plot the data
x11()
matplot(x, t(data), 
        type='l', lty=1, col=c('purple','forestgreen')[gender],
        main = "Berkeley Growth Study data",
        xlab = 'Age (yrs)', ylab = 'Height (cm)')
#for(i in 1:n)points(x, data[i,])
legend('topleft', bty = 'n', lwd=2, 
       col=c('purple','forestgreen'), legend = c('males', 'females'))

# plot with point-wise sample mean
matplot(x, t(data), type='l', lty=1, col='lightgrey',
        main = "Berkeley Growth Study data",
        xlab = 'Age (yrs)', ylab = 'Height (cm)')
sampleMean <- colMeans(data)
lines(x, sampleMean, lty=2, lwd=2, col='red')

## data are clearly monotone!!
# is there a way to exploit this knowledge in smoothing?
## monotone splines (see R&S chapter 6) ##
# -> define any basis to start with
# -> compute the smoother as the integral of exp(basis expansion)
# -> already implemented in the function "smooth.monotone"

## STEP 2.1 -- Monotone Smoothing Spline (for smoothing Growth data)
####################################################################

# set the initial necessary variables
x_iniz <- x
p <- 31 # number of age measurements
p.fine <- 110 # fine grid on which to evaluate the data -> a bit larger! (border effects)

# I will need to estimate the derivatives! 
# (both for the analysis and for the definition of the functional similarity)
deriv0 <- deriv1 <- deriv2 <- NULL

## create a B-spline basis
# -> knots every 2 years except at border
# -> better computational stability; border effects, stop at 16 years
age <- c(min(x):3,(2:7)*2+1,16)
age
agefine <- seq(min(x),16,length=p.fine)
fbobj <- create.bspline.basis(range(x),norder=4,breaks=age)
WfdParobj <- fdPar(fbobj)

for(i in 1:n){ # takes approx 10 minutes to run -> saved results!
  # monotone smoothing
  result <- smooth.monotone(x,data[i,],WfdParobj)
  wfd <- result$Wfdobj
  beta <- result$beta
  hgtfine <- beta[1]+beta[2]*eval.monfd(agefine,wfd)
  deriv0 <- rbind(deriv0,t(hgtfine))
  # evaluate first derivative -> growth velocity!
  hgtfine1 <- beta[2]*eval.monfd(agefine,wfd,1)
  deriv1 <- rbind(deriv1,t(hgtfine1))
  # evaluate second derivative -> growth acceleration!
  hgtfine2 <- beta[2]*eval.monfd(agefine,wfd,2)
  deriv2 <- rbind(deriv2,t(hgtfine2))
}

# remove border points (derivatives are much more noisy at borders!)
agefine <- agefine[-(1:10)]
deriv0 <- deriv0[,-(1:10)]
deriv1 <- deriv1[,-(1:10)]
deriv2 <- deriv2[,-(1:10)]
p.fine <- length(agefine)

# save the smoothed data so that we can use it later!
save(agefine,deriv0,deriv1,deriv2,file='Lab_2_growth_smoothing.Rdata')

## plot the smoothed data (same as lecture)
par(mfrow=c(1,3))
matplot(agefine, t(deriv0), type = 'l', col=c('purple','forestgreen')[gender],
        lty=1, xlab = 'Age (yrs)', ylab = 'Height (cm)')
legend('topleft', bty = 'n', lwd=2, 
       col=c('purple','forestgreen'), legend = c('males', 'females'))
matplot(agefine, t(deriv1), type = 'l', col=c('purple','forestgreen')[gender],
        lty=1, xlab = 'Age (yrs)', ylab = 'Velocity (cm/yr)')
matplot(agefine, t(deriv2), type = 'l', col=c('purple','forestgreen')[gender],
        lty=1, xlab = 'Age (yrs)', ylab = 'Acceleration (cm/yr^2)')
abline(h=0, lty=2)

# evaluate the estimates for some selected data
par(mfrow=c(1,2))
i <- 20 # 10, 15, 20
dt <- (x[-1]+x[-length(x)])/2
df <- diff(data[i,])/diff(x)
plot(dt, df, pch=16,xlab='age',ylab='growth velocity',ylim = range(c(df,deriv1[i,])),
     main=paste('Growth data: empirical 1st derivative for curve ', i,sep=''))
lines(agefine,deriv1[i,],col=2, lwd=2)
dt2 <- (dt[-1]+dt[-length(dt)])/2
df2 <- diff(df)/diff(dt)
plot(dt2, df2, pch=16,xlab='age',ylab='growth acceleration', ylim = range(c(df2,deriv2[i,])),
     main=paste('Growth data: empirical 2nd derivative for curve ', i,sep=''))
lines(agefine,deriv2[i,],col=2, lwd=2)

# basis seems to capture the data oscillations!
# this check is quite crucial as we fixed the number of
# basis functions (remember from lecture: nbasis = order + nr. internal knots)

# Exercise: try to use CV to optimize the number of basis functions for some
#           selected observations, and see what happens

## STEP 2.2 -- Functional Clustering via k-means
################################################

#install.packages("fdacluster")
?fdacluster

## Example: functional clustering of Growth velocities
?fdakmeans

# start with functional k-means 
growth.kmean <- fdakmeans(agefine, deriv1, n_clusters = 2, check_total_dissimilarity = FALSE,
                  centroid_type = "mean", warping_class = "none", metric = "l2")

x11()
plot(growth.kmean, type = "amplitude")

### how to decide on the number of groups?
### several cluster quality criteria are automatically stored in the result
# distances_to_center : vector of length n, distance of each curve to the center of its cluster
# silhouettes: vector of length n, silhouette values of each observation
# amplitude_variation: numeric value, fraction of total variation explained by amplitude variability

hist(growth.kmean$distances_to_center, xlab = expression(d(x[i], mu^c)),
     main = 'Histogram of distances to cluster centers')

hist(growth.kmean$silhouettes, xlab = expression(s[i]),
     main = 'Curves silhouette values')

which.L <- 2:7 # choose reasonable max number of clusters depending on n
growth.kmean.L <- NULL
for(L in which.L){# takes approx 10 minutes to run -> saved results!
  
  print('---------------------')
  print(paste('start on L = ', L,sep=''))
  out <- fdakmeans(agefine, deriv1, n_clusters = L, check_total_dissimilarity = FALSE,
                            centroid_type = "mean", warping_class = "none", metric = "l2")
  growth.kmean.L <- c(growth.kmean.L, list(out))

}
save(which.L, growth.kmean.L, file= 'Lab_2_functional_kmeans.Rdata')

# postprocess results
av.dist <- silh <- matrix(NA, n, length(which.L))
for(l in 1:length(growth.kmean.L)){
  av.dist[,l] <- growth.kmean.L[[l]]$distances_to_center
  silh[,l] <- growth.kmean.L[[l]]$silhouettes
  #ampl <- c(ampl, growth.kmean.L[[l]]$amplitude_variation)
}

## plot of the distance of each curve to the center of its cluster
## as a function of L (elbow plot)
x11()
#pdf('./images/Lect2_fClust_growth_kmean_avdist.pdf', width = 9, height = 9)
plot(which.L, which.L, pch='', ylim = range(av.dist), xlim = range(which.L)+c(-1,1),
     xlab = 'number of clusters', ylab = 'distance',
     main = 'Boxplots of distances to cluster centers')
for(l in 1:length(which.L))boxplot(av.dist[,l], at = which.L[l], add = TRUE)
#dev.off()

## plot of the silhouette value of each curve as a function of L (elbow plot)
x11()
#pdf('./images/Lect2_fClust_growth_kmean_silh.pdf', width = 9, height = 9)
plot(which.L, which.L, pch='', ylim = range(silh), xlim = range(which.L)+c(-1,1),
     xlab = 'number of clusters', ylab = 'silhouette',
     main = 'Boxplots of silhouette values')
for(l in 1:length(which.L))boxplot(silh[,l], at = which.L[l], add = TRUE)
#dev.off()

## STEP 2.3 -- Hierarchical Functional Clustering
#################################################
?fdahclust

growth.hclust <- fdahclust(agefine, deriv1, n_clusters = 2, check_total_dissimilarity = FALSE,
                          centroid_type = "mean", warping_class = "none", 
                          metric = "l2", linkage_criterion = "complete")
# number of clusters has to be specified, because the function needs it for alignment!
# (weakness of the coding.. it should return the dendrogram in case of no alignment)

x11()
plot(growth.hclust, type = "amplitude")


## STEP 2.4 -- Alternative Functional Clustering methods
########################################################


## funFEM 
#install.packages("funFEM")
?funFEM

# Clustering the Canadian Temperature data (they are included in the package)
basis <- create.bspline.basis(c(0, 365), nbasis=21, norder=4)
fdobj <- smooth.basis(day.5, CanadianWeather$dailyAv[,,"Temperature.C"],basis,
                      fdnames=list("Day", "Station", "Deg C"))$fd
res = funFEM(fdobj,K=2)

# Visualization of the partition and the group means
par(mfrow=c(1,2))
plot(fdobj); lines(fdobj,col=res$cls,lwd=2,lty=1)
fdmeans = fdobj; fdmeans$coefs = t(res$prms$my)
plot(fdmeans); lines(fdmeans,col=1:max(res$cls),lwd=2)

# Visualization in the discriminative subspace (projected scores)
par(mfrow=c(1,1))
plot(t(fdobj$coefs) %*% res$U,col=res$cls,pch=19,main="Discriminative space")


####################################
###       PART 3 -- AIM:         ###
### FUNCTIONAL DATA REGISTRATION ###
####################################


## STEP 3.1 -- Landmark registration 
####################################

# landmark registration can be performed with the function
# 'landmarkreg' in the package fda. Was including errors!
# -> my debugged version
source('Lab_2_my_landmarkreg.R')

## Example: Female Growth Velocities from the Berkeley Growth study
# to use the function 'landmarkreg', we need the smoothed functional
# objects: let us do the smoothing only on the first 10 girls
attach(growth)
female.sel <- hgtf[,1:10]
n.sel <- 10

# monotone smoothing with order 6 spline basis,
# knots at ages of observations, roughness penalty third derivatives
# and a smoothing parameter 1/sqrt{10}
wbasis <- create.bspline.basis(c(1,18), 35, 6, age)
growfdPar <- fdPar(wbasis, 3, 10^(-0.5))
growthMon <- smooth.monotone(age, female.sel, growfdPar)
Wfd <- growthMon$Wfd
betaf <- growthMon$beta
hgtfhatfd <- growthMon$yhatfd

## for landmark registration, we first need to define the landmarks!
## Good choices? We can reason in this way
# 1. single good landmark t_i is the age for girl i at which her
#    acceleration curve crossed 0 with a negative slope during the
#    pubertal growth spurt
# 2. we can also define t_0 as a fixed time specified for the middle
#    of the average pubertal growth spurt
# Then we specify warping functions h_i by fitting a smooth function
# to the three points (1,1); (t_0, t_i); and (18;18)

# first find the acceleration (2nd derivatives)
accelfdUN <- deriv.fd(hgtfhatfd, 2)

# then find the point t_i for each curve:
# the R function locator() can be applied to plots to select
# a particular abscissa point;
# in this case, to a functional data object 'accfd' that contains
# estimated acceleration curves to select the age
# of the center of the pubertal growth spurt
PGSctr <- rep(0,10)
newagefine <- seq(1,18,len=101)
for (icase in 1:n.sel) {
  accveci <- predict(accelfdUN[icase], newagefine)
  pos.spurt <- which(diff(accveci > 0) == -1)
  if(length(pos.spurt)>1)pos.spurt <- pos.spurt[length(pos.spurt)]
  PGSctr[icase] <- newagefine[pos.spurt+1]
}
PGSctrmean <- mean(PGSctr)
# the function fitting the three landmark points does not need to be
# *very flexible* -> order 3 spline basis functions 
wbasisLM <- create.bspline.basis(c(1,18), 4, 3,#norder=3, nbasis = 4,
                                 breaks=c(1,PGSctrmean,18))
WfdLM <- fd(matrix(0,4,1),wbasisLM)
WfdParLM <- fdPar(WfdLM,1,1e-12)

## actual landmark registration can now be performed:
regListLM <- mylandmarkreg(accelfdUN, cbind(PGSctr), PGSctrmean,
                         c(1,18), WfdPar = WfdParLM)
accelfdLM <- regListLM$regfd
accelmeanfdLM <- mean(accelfdLM)
warpfdLM <- regListLM$warpfd
WfdLM <- regListLM$Wfd

# plot the final result (misalinged and landmark aligned curves)
x11()
#pdf('./images/Lect2_Registr_landmark.pdf', width = 12, height =7)
par(mfrow=c(2,1))
newagefine <- seq(2,17.1,len=101)# remove border effects again
plot(newagefine,newagefine,pch="", ylim=c(-6,4),
     xlab="Year", ylab="Height Acceleration (cm/yr^2)",
     main='Misaligned Accelerations')
for(icase in 1:n.sel){
  accveci <- predict(accelfdUN[icase], newagefine)
  lines(newagefine,accveci) #, col=rainbow(12)[icase+1]
}
abline(h=0, lty=2)
abline(v=PGSctrmean, lty=2, lwd=2)
#abline(v=PGSctr, lty=2, col = rainbow(12)[2:11])
plot(newagefine,newagefine,pch="", ylim=c(-6,4),
     xlab="Year", ylab="Height Acceleration (cm/yr^2)",
     main='Landmark Aligned Accelerations')
for (icase in 1:n.sel) {
  accveci <- predict(accelfdLM[icase], newagefine)
  lines(newagefine,accveci)
}
abline(h=0, lty=2)
abline(v=PGSctrmean, lty=2, lwd=2)
#dev.off()

x11()
#pdf('./images/Lect2_Registr_landmark_w.pdf', width = 6, height = 6)
plot(warpfdLM, lty=1, xlab = 't', ylab = 'h(t)',
     main = 'Estimated warping functions')
abline(a=0, b=1, lwd=2)
#dev.off()


## STEP 3.2 -- Continuous registration 
######################################
## Back to the function 'fdacluster',
## which handles clustering and alignment
?fdacluster

## Example: functional clustering of Growth velocities
?fdakmeans

# start with functional k-means 
growth.kmean <- fdakmeans(agefine, deriv1, n_clusters = 2, check_total_dissimilarity = FALSE,
                          centroid_type = "mean", warping_class = "none", metric = "l2")

x11()
plot(growth.kmean, type = "amplitude")

