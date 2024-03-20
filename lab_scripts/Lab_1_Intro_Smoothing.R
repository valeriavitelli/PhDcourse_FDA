
rm(list=objects())
setwd("C:/Users/valeriv/UiO Dropbox/Valeria Vitelli/Teaching/FDA_Palermo/R scripts")
set.seed(24032024)

library(fda)
library(refund)
library(fields)
library(fdaoutlier)

###########################################
###           PART 1 -- AIM:            ###
### EXPLORE SEVERAL FUNCTIONAL DATASETS ###
###########################################

# 3 examples of functional datasets

#######################################
### Dataset 1: Berkeley Growth Data ###
#######################################

## BACKGROUND:
# heights of 54 girls and 39 boys
# measured at a set of 31 ages
# Berkeley Growth Study
# used mostly in Lecture 2! (we will see smoothing then)

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


################################################
### Dataset 2: Diffusion Tensor Imaging Data ###
################################################

## BACKGROUND:
# Tractography Application: Multiple sclerosis (MS) is an demyalinating 
# autoimmune disease associated with brain lesions.
# DTI is an MRI technique which measures proxies of demyalination
# by quantifying the water diffusivity in the brain.
# Fractional anisotropy (FA) is one measure of water diffusion in the brain.

## DATA:
# Here we consider FA tract profiles for the corpus callosum (CCA) measured
# without missing values for 66 MS subjects at their first visit.

data(DTI); attach(DTI)
names(DTI)

head(ID)
head(case)
head(visit)

DTI.complete <- subset(DTI, complete.cases(DTI))
DTI.baseline <- subset(DTI.complete, visit == 1 & case == 1)
n <- length(unique(DTI.baseline$ID))
dim(DTI.baseline$cca)
tract <- 1:dim(DTI.baseline$cca)[2]

# plot the data
#pdf('./images/Lect1_data_MS.pdf', width = 8, height = 4.5)
matplot(tract, t(DTI.baseline$cca), 
        type='l', lty=1, col=rainbow(n),
        main = "Diffusion Tensor Imaging : CCA",
        xlab="tract", ylab="Fractional anisotropy (FA)")
#dev.off()

# point-wise sample mean
matplot(tract, t(DTI.baseline$cca), type='l', lty=1, col='lightgrey',
        main = "Diffusion Tensor Imaging : CCA",
        xlab="tract", ylab="Fractional anisotropy (FA)")
sampleMean <- colMeans(DTI.baseline$cca)
lines(tract, sampleMean, lty=2, lwd=2, col='red')

# empirical covariance
sampleCov <- cov(DTI.baseline$cca)
image.plot(tract, tract, sampleCov, main='sample covariance of FA')

# empirical correlation
sampleCor <- cor(DTI.baseline$cca)
image.plot(tract, tract, sampleCor, main='sample correlation of FA')

## plot relations among variables:
# the DTI dataset also includes the patients’ Paced Auditory Serial Addition Test
# (PASAT) scores, a measure of cognitive functions. 
# Below we plot FA measurements again, but while also visualizing the relationship
# to the PASAT scores: by color-coding the PASAT scores (red = low, blue = high)
clrs <- rev(colorRampPalette(c("blue", "green", "yellow", "red"))(40))    
colfct <- as.numeric(cut(DTI.baseline$pasat, 40))
matplot(tract, t(DTI.baseline$cca), 
        type='l', lty=1, col = clrs[colfct],
        main = "Diffusion Tensor Imaging : CCA",
        xlab="tract", ylab="Fractional anisotropy (FA)")

# same plot, but now 3D
## 3D rainbow plot
par(mar=c(1,1,0,0), cex.axis=1, cex.lab=1)
clrs <- rev(colorRampPalette(c("blue", "green", "yellow", "red"))(40))    
proj = persp(x = tract, y = seq(min(DTI.baseline$pasat), max(DTI.baseline$pasat), l=length(DTI.baseline$pasat)),  z=t(DTI.baseline$cca),
             xlab="tract", ylab="PASAT", zlab="FA", col=NA, border=NA,
             ticktype = "detailed", axes=TRUE, theta=30, phi=30)
o <- rev(order(DTI.baseline$pasat))
for(i in o){
  lines(trans3d(x = tract, y=rep(DTI.baseline$pasat[i], ncol(DTI.baseline$cca)),  z=DTI.baseline$cca[i,], pmat=proj), col=clrs[colfct[i]])
}


########################################
### Dataset 3: Canadian Weather Data ###
########################################

data("CanadianWeather")       # load data
attach(CanadianWeather)       # attach data
?CanadianWeather
names(CanadianWeather)

# plot the average monthly temperature across locations -> already prepared
month <- 1:12
dim(monthlyTemp)
n <- ncol(monthlyTemp)
n
matplot(month, monthlyTemp,
        type='l', lty=1, col = rainbow(n),
        xlab="months", ylab="temperature", 
        main="monthly temperatures")

# plot the daily average temperature across locations -> select inside the array 'dailyAv'
day <- 1:365  # define functional argument
dim(dailyAv)

matplot(day, dailyAv[,,1], # see help
        type='l', lty=1, col = rainbow(n),
        xlab="days", ylab="temperature", 
        main="daily temperatures")

##################################################
## FUNCTIONAL BAND DEPTH AND FUNCTIONAL BOXPLOT ##
##################################################

?band_depth
weatData_depth <- band_depth(t(dailyAv[,,1])) # J = 2
weatData_depth

?fbplot
#pdf('./images/Lect1_boxplot_weatCanad_mod.pdf', width = 8, height = 4.5)
fbplot(dailyAv[,,1], method = 'MBD', ylim = range(dailyAv[,,1]),
       xlab = 'days', ylab = 'temperature', main = 'Functional Boxplot for Canadian Weather data')
dev.off()


#####################################
###         PART 2 -- AIM:        ###
### SMOOTHING OF FUNCTIONAL DATA  ###
#####################################

## STEP 2.1: Construct different basis functions in R
#####################################################

## 1st approach: the fda-package
## allows constructing several bases, returns functional data objects

library(fda)
?"fda-package"

## monomial basis: requires the domain and the basis number
# Example: monomial basis K=8 basis functions defined on the interval [0,1]
K <- 8
bbasis_obj <- create.monomial.basis(rangeval=c(0,1), nbasis = K)

# this returns a functional object that needs to be evaluated
x <- seq(0,1,length.out=100)
bbasisevals <- eval.basis(x, bbasis_obj)
dim(bbasisevals) # m x K matrix, m = number of evaluations, K = basis dimension
matplot(x, bbasisevals, type='l', lty=1, col=rainbow(6),
        xlab="x", ylab="basis functions", 
        main=paste('monomial basis with K = ', K, sep=''))

## Fourier basis:
# need to define the domain, the period of oscillation, and the number of basis functions
fbasis_obj <- create.fourier.basis(rangeval=c(0,1), 
                                   nbasis=65, period = 1)
fbasisevals <- eval.basis(x, fbasis_obj)
matplot(x, fbasisevals[, 1:5], type='l', lty=1, col=rainbow(5),
        xlab="x", ylab="basis functions", 
        main="first 5 Fourier basis functions")

## B-spline basis:
# need to specify the domain, the number basis functions, and the order
K <- 12; Bsp.order <- 4
bsbasis_obj <- create.bspline.basis(rangeval=c(0,1),
                                    nbasis=K, norder=Bsp.order)
bsbasisevals <- eval.basis(x, bsbasis_obj)
matplot(x, bsbasisevals, type='l', lty=1, col=rainbow(K+4),
        xlab="x", ylab="basis functions", 
        main=paste('First ',K, ' elemens of the B-spline basis of order ', Bsp.order, sep = '') )

## Exercise: try to vary K and/or Bsp.order and see the effect

## STEP 2.2: Create functional data
###################################

## Functional Data (FD) Object in the package fda
# Need a basis of functions and a vector of coefficients  
# The function 'fd' is needed to build functional data objects
  
nb <- 10
norder <- 4
coef <- rnorm(nb)
bsbasis_obj <- create.bspline.basis(rangeval=c(0,1),
                                    nbasis=nb, norder=norder)
fd_obj <- fd(coef, bsbasis_obj)
x <- seq(0,1,length.out = 100)
fd_eval <- eval.fd(x, fd_obj)
plot(x, fd_eval, type='l', ylab="", xlab="x", 
     main=paste('Random function generated using ',nb,' B-splines of order ',norder, sep='') )

# Exercise: try running the code several times; try to vary the different choices

## Sample of n random functions: need a matrix of random coefficients (dimension K x n)
## and the order of B-spline basis of choice (here, cubic)
n <- 20 ; nb <- 10; norder <- 4
coefs <- matrix(rnorm(n*nb),nb,n)
dim(coefs)   # [1]  10 20

bsbasis_obj <- create.bspline.basis(rangeval=c(0,1),
                                    nbasis=nb, norder=norder)
fd_obj <- fd(coef, bsbasis_obj)
x <- seq(0,1,length.out = 100)
fd_eval <- eval.fd(x, fd_obj)
dim(fd_eval)   # [1]  100 20

matplot(x, fd_eval, 
        type='l', lty = 1, col=rainbow(n),
        ylab='f(x)', xlab='x', 
        main='Sample of random functions via a B-spline basis')

## STEP 2.3: Linear regression on basis functions (OLS)
#######################################################

data("CanadianWeather")

## Let’s focus on log-transformed average daily precipitation
## of Victoria (one observed curve)

y.lprec <- CanadianWeather$dailyAv[,,3]
l <- which(CanadianWeather$place == "Victoria")
y <- y.lprec[,l]
day <- 1:365
plot(day, y, 
     type='o', pch = 16, cex = 0.5, col='forestgreen',
     xlab="day", ylab="log-precipitation",
     main="Log of daily average precipitation in Victoria")


## Periodic data -> Data smoothing via a Fourier basis

# define the domain, period, and the number of basis
rangval <- range(day)
period <- 365 
nbasis <- 5 

# create the basis
fbasis <- create.fourier.basis(rangval, nbasis=nbasis, period=period)  
bvals <- eval.basis(day, fbasis)
Xbasis <- bvals

# smoothing via LS: fits linear regression on basis functions
lm.fit <- lm(y ~ 0 + Xbasis)   
y.fit <- lm.fit$fitted.values
coef <- lm.fit$coefficient

# Additional information one can derive from smoothing:
# Second Derivative of the fit evaluated at day
yfitfd <- fd(coef,fbasis)  # obtain the functional data object
yfit2D <- eval.fd(day, yfitfd, 2) # evaluate the 2nd derivative

x11()
par(mfrow=c(2,1))
plot(day, y, type="n",lwd=4, col="black",
     xlab="day", ylab="log-precipitation", 
     main=paste(nbasis, "Fourier basis functions"), cex=1)
points(day, y, pch=1, cex=.5, col="blue", lwd=1)
lines(day, lm.fit$fitted.values, lwd=1, col="red")


plot(day, yfit2D, type="l",lwd=2, col="black", 
     xlab="day", ylab="Second derivative of log-precipitation", 
     main=paste("L2 Norm of second derivative = ", 
                round(mean(yfit2D^2),2)))

# Exercise: try to experiment by changing the number of basis functions, 
#           and examine the impact


## STEP 2.4: How to choose the number of basis functions?
#########################################################

## Cross-validation: useful to choose the number of basis functions 

# choose the same setting as before
attach(CanadianWeather)
y.precip <- dailyAv[,,2]
l <- which(place=='Victoria') 
t.day <- 1:365  
y <- y.precip[,l]
p <- length(y)

# choose which values for the number of basis functions to try
which.K <-  2*c(2:11)+1 
which.K

# calculate all Fourier bases that are needed
t.day <- 1:365
fbasis <- create.fourier.basis(rangeval = c(1, 365), nbasis=max(which.K), period=365)
bvals <- eval.basis(t.day,fbasis)


# compute CV score for each t_j and each K
CVfit <-  matrix(0,  nrow=p, ncol=length(which.K))
for(j in 1:p){
  
  Y.star <- y[-j] #L-O-O cross validation
  
  # fit using Fourier basis and K basis functions
  index <- 0
  for (K in which.K){
    index <- index+1
    Xbasis <- bvals[, 1:K]
    Xbasis.j <- Xbasis[-j, ]  # remove t_j
    lm.fit <-  lm(Y.star~0+Xbasis.j)
    Xbasis.coeff <- lm.fit$coefficients
    y.fit <- Xbasis%*%Xbasis.coeff
    CVfit[j,index] <- (y[j] - y.fit[j])^2 # compute CV score
  }
}
CV_L2 <- apply(CVfit, 2, sum) # integral of the CV score over T

# plot the CV score and choose K
x11()
plot(which.K, CV_L2, type='n',
     xlab='Number of basis functions', 
     ylab='Total cross-validation error')
points(which.K, CV_L2, type='b', col=2, lwd=2)
title(paste0("K = ", which.K[which(CV_L2==min(CV_L2))], " with the smallest CV score!"))

myK <- c(13,21) # choose 2 values to compare
x11()
par(mfrow = c(1,2))
for(k in myK){
  Xbasis <- bvals[, 1:k]
  lm.fit <-  lm(y~0+Xbasis)
  Xbasis.coeff <- lm.fit$coefficients
  y.fit <- Xbasis%*%Xbasis.coeff
  plot(day, y, type="n",lwd=4, col="black",
       xlab="day", ylab="log-precipitation", 
       main=paste('Smoothing with ',k,' Fourier basis functions'), cex=1)
  points(day, y, pch=1, cex=.5, col="blue", lwd=1)
  lines(day, lm.fit$fitted.values, lwd=1, col="red")
}


# Exercise: Try to change location (Vancouver!) and see what happens
#           do you get a good optimal behavior of the CV score?


## STEP 2.5: Penalized smoothing & selection of smoothing parameter
###################################################################

## Smoothing with Roughness Penalty:
# define the smoothness of a function and control it through the penalty parameter

# Fix the penalty parameter first: 
# Focus on Victoria precipitations
# use a B-spline basis for smoothing with LOTS of basis elements (one for each day!)

y.precip=dailyAv[,,2]
l <-  which(place=="Victoria") 
t.day <- 1:365  
y <- y.precip[,l]

nbasis <- 365
ybasis  <- create.bspline.basis(rangeval = c(1,365), nbasis = nbasis, norder=4)
bvals <- eval.basis(t.day, ybasis)
Xbasis  <- bvals
lm.fit <- lm(y ~ 0 + Xbasis)   
y.fit <-  lm.fit$fitted.values

x11()
plot(t.day, y, type="n",lwd=4, col="black",
     xlab="day", ylab="precipitation", 
     main=paste('LS smoothing with ', nbasis, ' B-spline functions'), cex=1)
points(t.day, y, pch=1, cex=.5, col="blue", lwd=1)
lines(t.day, lm.fit$fitted.values, lwd=1, col="red")
# Without roughness penalty the OLS fit connects each of the measurements


# We now use the penalty parameter to control roughness of the fit
# and see the effect on smoothness
lambda <- 10^4

## some useful functions and arguments..
# int2Lfd(m): defines the m-th order derivative penalty term
# fdPar(): sets functional parameters; 2nd order derivative penalty & smoothing parameter
# smooth.basis(): smoothing with roughness penalty & smoothing parameter as in 'tD2fdPar' 

ybasis  <- create.bspline.basis(rangeval = c(1,365), nbasis = 365, norder=4)
tD2fdPar <- fdPar(ybasis, Lfdobj=int2Lfd(2), lambda=lambda)
tyfd <- smooth.basis(t.day,y,tD2fdPar) 
names(tyfd)
# fd:   functional data object containing the smoothed data
# df:   degrees of freedom
# gcv:  value of the generalized cross-validation criterion 
# beta: regression coefficients (i.e., basis coefficient)
# SSE:  error sums of squares
# penmat: penalty matrix
# y2cMap: matrix mapping the data to the coefficients (Phi^T Phi + R)^(-1) \Phi^T

x11()
plot(t.day, y, type="n", ylim=range(y), 
     ylab="Precipitation", xlab="day", 
     main= paste("Victoria data (lambda =", round(lambda,2), ")", sep=""))
points(t.day, y, pch=1, cex=.5, col="blue", lwd=1)
lines(tyfd$fd,col="red",lwd=4)

## Exercise 1: Try to change the value of lambda to very small (ex. 0.0001)
##             or much larger (ex. 1e8) values, and see what happens
## Exercise 2: Try to change location (Vancouver!) and see what happens


