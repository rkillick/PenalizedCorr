library(tinytest)


set.seed(123)
y <- arima.sim(n=100, model=list(ar=0.5))
#also checked
# y <- arima.sim(n=100, model=list(ma=0.5))
y1 <- c("a", "b", "c", "d", "e", "f")
y2 <- c("3", "1", "5", "3", "8", "6")

#################################################
## testing that PenalizedCorr can only be logical
#################################################
### acf ###
expect_error(acf(y,penalized =NULL), "penalized must be logical")
expect_error(acf(y,penalized = 3), "penalized must be logical")
expect_error(acf(y,penalized = "afkj"), "penalized must be logical")

### pacf ###
expect_error(pacf(y,penalized =NULL), "penalized must be logical")
expect_error(pacf(y,penalized = 3), "penalized must be logical")
expect_error(pacf(y,penalized = "adkj"), "penalized must be logical")

### ar ###
expect_error(ar(y,penalized =NULL), "penalized must be logical")
expect_error(ar(y,penalized = 3), "penalized must be logical")
expect_error(ar(y,penalized = "sfdlkj"), "penalized must be logical")

### auto.acf ###
expect_error(auto.acf(y,penalized =NULL), "penalized must be logical")
expect_error(auto.acf(y,penalized = 3), "penalized must be logical")
expect_error(auto.acf(y,penalized = "sfdkjl"), "penalized must be logical")

### invertpacf ###
expect_error(PenalizedCorr:::invertpacf(y,penalized =NULL), "penalized must be logical")
expect_error(PenalizedCorr:::invertpacf(y,penalized = 3), "penalized must be logical")
expect_error(PenalizedCorr:::invertpacf(y,penalized = "sdklj"), "penalized must be logical")

### DLpencoaf ###
expect_error(PenalizedCorr:::DLpencoef(y,penalized =NULL), "penalized must be logical")
expect_error(PenalizedCorr:::DLpencoef(y,penalized = 3), "penalized must be logical")
expect_error(PenalizedCorr:::DLpencoef(y,penalized = "fdlkj"), "penalized must be logical")


### ar.penyw ###
expect_error(PenalizedCorr:::ar.penyw(y,penalized =NULL), "penalized must be logical")
expect_error(PenalizedCorr:::ar.penyw(y,penalized = 3), "penalized must be logical")
expect_error(PenalizedCorr:::ar.penyw(y,penalized = "sdfklj"), "penalized must be logical")



###############################################
## testing that data series can only be numeric
###############################################
expect_error(acf(y1), "'x' must be numeric")
expect_error(acf(y2), "'x' must be numeric")

expect_error(pacf(y1), "'x' must be numeric")
expect_error(pacf(y2), "'x' must be numeric")

expect_error(auto.acf(y1), "'x' must be numeric")
expect_error(auto.acf(y2), "'x' must be numeric")

expect_error(PenalizedCorr:::invertpacf(y1), "'x' must be numeric")
expect_error(PenalizedCorr:::invertpacf(y2), "'x' must be numeric")

expect_error(PenalizedCorr:::DLpencoef(y1), "'x' must be numeric")
expect_error(PenalizedCorr:::DLpencoef(y2), "'x' must be numeric")

expect_error(ar(y1), "'x' must be numeric")
expect_error(ar(y2), "'x' must be numeric")

expect_error(PenalizedCorr:::ar.penyw(y1), "'x' must be numeric")
expect_error(PenalizedCorr:::ar.penyw(y2), "'x' must be numeric")


##################################
## testing lh, numeric + dimension
##################################
lh=c(0.4,0.2)
expect_error(acf(y, lh = lh))
expect_error(pacf(y, lh = lh))
expect_error(auto.acf(y, lh = lh))
expect_error(PenalizedCorr:::invertpacf(y, lh = lh))
expect_error(PenalizedCorr:::DLpencoef(y, lh = lh))
expect_error(ar(y, lh = lh))
expect_error(PenalizedCorr:::ar.penyw(y, lh = lh))

lh="sfdjkl"
expect_error(acf(y, lh = lh),"lh must be numeric")
expect_error(pacf(y, lh = lh),"lh must be numeric")
expect_error(auto.acf(y, lh = lh),"lh must be numeric")
expect_error(PenalizedCorr:::invertpacf(y, lh = lh),"lh must be numeric")
expect_error(PenalizedCorr:::DLpencoef(y, lh = lh),"lh must be numeric")
expect_error(ar(y, lh = lh),"lh must be numeric")
expect_error(PenalizedCorr:::ar.penyw(y, lh = lh),"lh must be numeric")

lh=TRUE
expect_error(acf(y, lh = lh),"lh must be numeric")
expect_error(pacf(y, lh = lh),"lh must be numeric")
expect_error(auto.acf(y, lh = lh),"lh must be numeric")
expect_error(PenalizedCorr:::invertpacf(y, lh = lh),"lh must be numeric")
expect_error(PenalizedCorr:::DLpencoef(y, lh = lh),"lh must be numeric")
expect_error(ar(y, lh = lh),"lh must be numeric")
expect_error(PenalizedCorr:::ar.penyw(y, lh = lh),"lh must be numeric")


########################################
## testing lambda, numeric and dimension
########################################
lambda=c(0.4,0.2)
expect_error(acf(y, lambda = lambda))
expect_error(pacf(y, lambda = lambda))
expect_error(PenalizedCorr:::DLpencoef(y, lambda = lambda))
expect_error(ar(y, lambda = lambda))
expect_error(PenalizedCorr:::ar.penyw(y, lambda = lambda))

lambda="sfdjkl"
expect_error(acf(y, lambda = lambda),"lambda must be numeric")
expect_error(pacf(y, lambda = lambda),"lambda must be numeric")
expect_error(PenalizedCorr:::DLpencoef(y, lambda = lambda),"lambda must be numeric")
expect_error(ar(y, lambda = lambda),"lambda must be numeric")
expect_error(PenalizedCorr:::ar.penyw(y, lambda = lambda),"lambda must be numeric")

lambda=TRUE
expect_error(acf(y, lambda = lambda),"lambda must be numeric")
expect_error(pacf(y, lambda = lambda),"lambda must be numeric")
expect_error(PenalizedCorr:::DLpencoef(y, lambda = lambda),"lambda must be numeric")
expect_error(ar(y, lambda = lambda),"lambda must be numeric")
expect_error(PenalizedCorr:::ar.penyw(y, lambda = lambda),"lambda must be numeric")


###########################################
## testing that lambda can only be positive
###########################################
expect_error(acf(y, lambda = -1), "lambda must be positive")
expect_error(pacf(y, lambda = -1), "lambda must be positive")
expect_error(PenalizedCorr:::DLpencoef(y, lambda = -1),"lambda must be positive")
expect_error(ar(y, lambda = -1),"lambda must be positive")
expect_error(PenalizedCorr:::ar.penyw(y, lambda = -1),"lambda must be positive")


########################################
## testing target, numeric and dimension
########################################
target=c(0.4,0.2)
expect_error(acf(y, target = target))
expect_error(pacf(y, target = target))
expect_error(PenalizedCorr:::DLpencoef(y, target = target))
expect_error(ar(y, target = target))
expect_error(PenalizedCorr:::ar.penyw(y, target = target))

target="sfdjkl"
expect_error(acf(y, target = target),"target must be numeric")
expect_error(pacf(y, target = target),"target must be numeric")
expect_error(PenalizedCorr:::DLpencoef(y, target = target),"target must be numeric")
expect_error(ar(y, target = target),"target must be numeric")
expect_error(PenalizedCorr:::ar.penyw(y, target = target),"target must be numeric")

target=TRUE
expect_error(acf(y, target = target),"target must be numeric")
expect_error(pacf(y, target = target),"target must be numeric")
expect_error(PenalizedCorr:::DLpencoef(y, target = target),"target must be numeric")
expect_error(ar(y, target = target),"target must be numeric")
expect_error(PenalizedCorr:::ar.penyw(y, target = target),"target must be numeric")


###################################################
## testing that target can only be between 1 and -1
###################################################
expect_error(acf(y, target = -2), "target must be between 1 and -1")
expect_error(acf(y, target = 2), "target must be between 1 and -1")
expect_error(pacf(y, target = -2), "target must be between 1 and -1")
expect_error(pacf(y, target = 2), "target must be between 1 and -1")
expect_error(PenalizedCorr:::DLpencoef(y, target = 2), "target must be between 1 and -1")
expect_error(PenalizedCorr:::DLpencoef(y, target = -2), "target must be between 1 and -1")
expect_error(ar(y, target = 2), "target must be between 1 and -1")
expect_error(ar(y, target = -2), "target must be between 1 and -1")
expect_error(PenalizedCorr:::ar.penyw(y, target = 2), "target must be between 1 and -1")
expect_error(PenalizedCorr:::ar.penyw(y, target = -2), "target must be between 1 and -1")


