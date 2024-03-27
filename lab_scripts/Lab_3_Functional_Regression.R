rm(list=objects())
setwd("C:/Users/valeriv/UiO Dropbox/Valeria Vitelli/Teaching/FDA_Palermo/R scripts")
set.seed(27032024)

library(fda)
library(refund)

#####################################
###       PART 1 -- AIM:          ###
### SCALAR-ON-FUNCTION REGRESSION ###
#####################################

### Example 1.1: Canadian Weather Data
## Question of interest:
## can the daily temperature along the year in each weather station
## predict the log10 of the total annual precipitation in that station?

?CanadianWeather
# recall the structure of the Canadian Weather Data:
# 'dailyAv' = a three dimensional array c(365, 35, 3) summarizing
#             data collected at 35 different weather stations in Canada
#             on the following measures:
# [,,1] = [,,'Temperature.C']: average daily temperature for each day of the year
# [,,2] = [,, 'Precipitation.mm']: average daily rainfall for each day of the year
#         (rounded to 0.1 mm)
# [,,3] = [,, 'log10precip']: base 10 logarithm of Precipitation.mm 
#         (after first replacing 27 zeros by 0.05 mm (see Ramsay and Silverman 2006, p. 248))

## let us fit a S-o-F regression model on the Canadian Weather data 
## we use the function 'fRegress'

#  set up log10 of annual precipitation for 35 weather stations
annualprec <- log10(apply(CanadianWeather$dailyAv[,,"Precipitation.mm"], 2,sum))
# The simplest 'fRegress' call is singular with more bases
# than observations, so we use only 25 basis functions
# (still a quite large basis for the 35 stations)
smallbasis  <- create.fourier.basis(c(0, 365), 25)
# The covariate is the temperature curve for each station.
tempfd <- smooth.basis(day.5, CanadianWeather$dailyAv[,,"Temperature.C"], smallbasis)$fd

precip.Temp1 <- fRegress(annualprec ~ tempfd)
#  the output is a list with class name fRegress, display names
names(precip.Temp1)

#  the vector of fits to the data is object  precip.Temp1$yfdPar,
#  but since the dependent variable is a vector, so is the fit
annualprec.fit1 <- precip.Temp1$yhatfdobj
#  plot the data and the fit
plot(annualprec.fit1, annualprec, type="p", pch="o", xlab = 'fitted', ylab = 'observed')
lines(annualprec.fit1, annualprec.fit1, lty=2)

#  print root mean squared error
RMSE <- sqrt(mean((annualprec-annualprec.fit1)^2))
print(paste("RMSE =",RMSE))

##  plot the estimated regression coefficient (functional)
plot(precip.Temp1$betaestlist[[2]], xlab = 't', ylab = expression(beta(t)))
# this is NOT AT ALL interpretable, we are clearly overfitting the data


## alternative way for carrying out *the exact same analysis*
precip.Temp.mdl1 <- fRegress(annualprec ~ tempfd, method="model")
precip.Temp.m <- do.call('fRegress', precip.Temp.mdl1)
all.equal(precip.Temp.m, precip.Temp1)

# now we can set up a smaller basis for beta2 than for temperature so that we
#  get a more parsimonious fit to the data

nbetabasis2 <- 15  #  a bit less, and we add some roughness penalization
betabasis2  <- create.fourier.basis(c(0, 365), nbetabasis2)
betafd2     <- fd(rep(0, nbetabasis2), betabasis2)

# add penalization
betafdPar2  <- fdPar(betafd2, lambda=10)

# replace the regress coefficient function with this fdPar object
precip.Temp.mdl2 <- precip.Temp.mdl1
precip.Temp.mdl2[['betalist']][['tempfd']] <- betafdPar2

# Now do re-fit the data
precip.Temp2 <- do.call('fRegress', precip.Temp.mdl2)
# Compare the two fits:
#  degrees of freedom
precip.Temp1[['df']] # 26
precip.Temp2[['df']]  # a bit more than 15
#  root-mean-squared errors:
RMSE1 <- sqrt(mean(with(precip.Temp1, (yhatfdobj-yfdobj)^2)))
RMSE1
RMSE2 <- sqrt(mean(with(precip.Temp2, (yhatfdobj-yfdobj)^2)))
RMSE2
# error is slightly increased

#  display further results for the more parsimonious model
annualprec.fit2 <- precip.Temp2$yhatfdobj
plot(annualprec.fit2, annualprec, type="p", pch="o", 
     xlab = 'fitted', ylab = 'observed')
lines(annualprec.fit2, annualprec.fit2, lty=2)

#  plot the estimated regression function
plot(precip.Temp2$betaestlist[[2]])
#  still a bit overfitted, but at least
#  now we see that it is primarily the temperatures in the
#  early winter that provide the fit to log precipitation by temperature


#####################################
###       PART 2 -- AIM:          ###
### FUNCTION-ON-SCALAR REGRESSION ###
#####################################

## Example 2.1: Canadian Wheather Data (again)
## But now the question of interest is:
## What is the association between total yearly precipitation
## and daily temperature curve?

## We can use a F-o-S regression model to answer this question
## Temp_i(t)=\beta_0(t)+\beta_1(t)*'TotalPreci' + \epsilon_i(t)

# First define response and predictor
data("CanadianWeather")
attach(CanadianWeather)
day <- 1:365

# therefore: Y is the daily temperature (for each station)
Y <- t(as.matrix(dailyAv[,,1])) 
# X is the total annual precipitation in each station
# (mean normalized, so that \beta(t) can be interpreted)
X <- as.vector(colSums(dailyAv[,,2])) - mean(as.vector(colSums(dailyAv[,,2])))

myDat <- data.frame(X = X)
myDat$Y <- Y

dim(Y);length(X) # sanity check
## [1]  35 365
## [1] 35

### ALTERNATIVE TO THE FDA PACKAGE:
## we fit a F-o-S regression model using the pffr function
## (The pffr function in the 'refund' package can fit any
## functional linear model with functional response)

fit <- pffr(Y ~ X, data = myDat)
yhat <- predict(fit, newdata = myDat)

Rsq_t <- 1-colSums((Y - yhat)^2) / colSums((Y - colMeans(Y))^2)

x11()
plot(day, Rsq_t, type='l', ylab = expression(R^2),
     main = 'Time dependent goodness-of-fit measure')
mean(Rsq_t)
## [1] 0.7296885

# plot the original and fitted curves
x11()
matplot(day, t(Y), type='l', lty = 1, 
        col='light grey', ylab="fitted", xlab="day")
matlines(day, t(yhat), type='l', lty = 1)

# same plot but now dependence on precipitation is explicit:
# yellow = dry stations
# blue   = wet stations
# which means that the colder places are more dry (makes sense!)
clrs <- rev(colorRampPalette(c("blue", "green", "yellow"))(20))
colfct <- as.numeric(cut(X, 20))
x11()
matplot(day, t(Y), type='l', lty = 1, 
        col='light grey', ylab="fitted", xlab="day")
matlines(day, t(yhat), type='l', lty = 1, col = clrs[colfct])



## Compute the estimated coefficient functions
# (in the F-o-S regression, two functional coefficients
#  must explain all variability in the response)
coef <- coef(fit)
beta0.hat <- coef$smterms$`Intercept(yindex)`$coef
beta1.hat <- coef$smterms$`X(yindex)`$coef

# plot the original curves with the estimated intercept function
# (\beta_0(t))
matplot(day, t(Y), type='l', lty = 1, col='light grey', 
        ylab=expression(paste(beta[0](t))), xlab="day")
lines(beta0.hat$yindex.vec, beta0.hat$value, type='l', lwd=2)
# clearly the intercept reflects a general behavior of temperature
# along the year in Canada (colder in winter  and warmer in summer)

plot(beta1.hat$yindex.vec, beta1.hat$value, type='l', lwd=2, 
     ylab=expression(paste(beta[1](t))),
     xlab="day")
# the behavior of \beta_1(t) is instead the opposite! why?
# this tells us that wet places are warmer in winter and colder in summer
# (stations along the coast, where temp is more stable along the year)
# while dry places are warmer in summer and colder in winter
# (inland areas, where summer can be *very* hot and winter *very* cold)
# Mean-normalizing the precipitation here helps interpretation *a lot*!!


## Other functions that can fit F-o-S regression:

# -> 'bayes_fosr' and 'fosr' in the refund package
# fosr takes as inputs a matrix or a fd object from the fda package
# it selects optimal smoothing parameters using various methods
# (GCV, REML, ML..)
?fosr

# -> 'fRegress' in the fda package 
# fRegress is a very flexible and useful function, however,
# the selection of smoothing parameters is not implemented inside
# (one needs to do it in advance - as we have seen for the S-o-F model)
?fRegress

## bayes_fosr uses Bayesian estimation 
## plot_shiny takes the output of bayes_fosr for interactive plots.
#fit <- bayes_fosr(Y ~ X, data = myDat)
#plot_shiny(fit)

## Visualization using refund.shiny package (document)
## For installation of refund.shiny see:
## https://cran.r-project.org/web/packages/refund.shiny/refund.shiny.pdf

#######################################
###       PART 3 -- AIM:            ###
### FUNCTION-ON-FUNCTION REGRESSION ###
#######################################


## Function-on-function (FoF) regression - Concurrent model
## The concurrent functional regression model

## Temp_i(t) = \beta_0(t) + \beta_1(t)*Precip(t) + \epsilon_i(t)

## relates daily average temperature at the current time point t
## to daily average precipitation at the same time point t

day <- 1:365
# temperature data
Y <- t(as.matrix(dailyAv[,,1]))
#fPCA basis expansion with very high PVE for the precipitation function:
fit <- fpca.sc(t(as.matrix(dailyAv[,,2])), pve=0.99) 
X <- fit$Yhat 

myDat <- list()
myDat$X <- X
myDat$Y <- Y

# use again pffr as easier
fit <- pffr(Y ~ X, data = myDat)
yhat <- predict(fit, newdata = myDat)
Rsq_t <- 1-colSums((Y - yhat)^2) / colSums((Y - colMeans(Y))^2)
mean(Rsq_t)

matplot(day, t(Y), type='l', lty = 1, 
        col='light grey', ylab="fitted", xlab="day")
matlines(day, t(yhat), type='l', lty=1)

## This model explains about 74%
## of total variability. The estimated coefficient functions are

coef <- coef(fit)
beta0.hat <- coef$smterms$`Intercept(yindex)`$coef
beta1.hat <- coef$smterms$`X(yindex)`$coef

# intercept function
matplot(day, t(Y), type='l', lty = 1, 
        col='light grey',
        ylab=expression(paste(beta[0](t))), xlab="day")
lines(beta0.hat$yindex.vec, beta0.hat$value, type='l', lwd=2, lty=1)

plot(beta1.hat$yindex.vec, beta1.hat$value, type='l', lwd=2, lty=1,
     ylab=expression(paste(beta[1](t))),
     xlab="day")

## The 'fRegress' function in the fda package also fits concurrent
## functional regression but selecting smoothing parameters is not
## implemented in the function (same problem as before)


### EXTRA EXERCISE
##  Analyze the DTI data using the regression models that we learned today:
## 1. What is the association between FA and PASAT scores of MS patients at their first visit? (scalar-on-function regression)
## 2. The dataset also includes rcst - FA profiles collected from the right cortisopinal tract.
## How do these measurements relate to FA measurements along CCA?

