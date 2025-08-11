library(tinytest)


set.seed(123)
y <- arima.sim(n=100, model=list(ar=0.5))
#also checked
# y <- arima.sim(n=100, model=list(ma=0.5))
y1 <- c("a", "b", "c", "d", "e", "f")
y2 <- c("a", "1", "c", "3", "e", "6")
acf(y,penalized =TRUE)

## testing that PenalizedCorr can only be logical
### acf ###
expect_error(acf(y,penalized =TRUE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(acf(y,penalized =FALSE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(acf(y,penalized =NULL), "penalized must be logical")
#Passed as errors found as expected
expect_error(acf(y,penalized = 3), "penalized must be logical")
#Passed as errors found as expected

### pacf ###
expect_error(pacf(y,penalized =TRUE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(pacf(y,penalized =FALSE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(pacf(y,penalized =NULL), "penalized must be logical")
#Passed as errors found as expected
expect_error(pacf(y,penalized = 3), "penalized must be logical")
#Passed as errors found as expected

### ar ###
expect_error(ar(y,penalized =TRUE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(ar(y,penalized =FALSE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(ar(y,penalized =NULL), "penalized must be logical")
#Passed as errors found as expected
expect_error(ar(y,penalized = 3), "penalized must be logical")
#Passed as errors found as expected

### auto.acf ###
expect_error(auto.acf(y,penalized =TRUE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(auto.acf(y,penalized =FALSE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(auto.acf(y,penalized =NULL), "penalized must be logical")
#Passed as errors found as expected
expect_error(auto.acf(y,penalized = 3), "penalized must be logical")
#Passed as errors found as expected

### invertpacf ###
expect_error(invertpacf(y,penalized =TRUE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(invertpacf(y,penalized =FALSE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(invertpacf(y,penalized =NULL), "penalized must be logical")
#Passed as errors found as expected
expect_error(invertpacf(y,penalized = 3), "penalized must be logical")
#Passed as errors found as expected

### DLpencoaf ###
expect_error(DLpencoef(y,penalized =TRUE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(DLpencoef(y,penalized =FALSE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(DLpencoef(y,penalized =NULL), "penalized must be logical")
#Passed as errors found as expected
expect_error(DLpencoef(y,penalized = 3), "penalized must be logical")
#Passed as errors found as expected


### ar.penyw ###
expect_error(ar.penyw(y,penalized =TRUE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(ar.penyw(y,penalized =FALSE), "penalized must be logical")
#Failed as no errors found as expected
expect_error(ar.penyw(y,penalized =NULL), "penalized must be logical")
#Passed as errors found as expected
expect_error(ar.penyw(y,penalized = 3), "penalized must be logical")
#Passed as errors found as expected


## testing that data series can only be numeric
#### acf ####
expect_error(acf(y), "'x' must be numeric")
#Failed as no error found as expected
expect_error(acf(y1), "'x' must be numeric")
#Passed as error found as expected
expect_error(acf(y2), "'x' must be numeric")
#Passed as error found as expected

#### pacf ####
expect_error(pacf(y), "'x' must be numeric")
#Failed as no error found as expected
expect_error(pacf(y1), "'x' must be numeric")
#Passed as error found as expected
expect_error(pacf(y2), "'x' must be numeric")
#Passed as error found as expected

#### auto.acf ####
expect_error(auto.acf(y), "'x' must be numeric")
#Failed as no error found as expected
expect_error(auto.acf(y1), "'x' must be numeric")
#Passed as error found as expected
expect_error(auto.acf(y2), "'x' must be numeric")
#Passed as error found as expected

#### invertpacf ####
expect_error(invertpacf(y), "'x' must be numeric")
#Failed as no error found as expected
expect_error(invertpacf(y1), "'x' must be numeric")
#Passed as error found as expected
expect_error(invertpacf(y2), "'x' must be numeric")
#Passed as error found as expected

#### DLpencoef ####
expect_error(DLpencoef(y), "'x' must be numeric")
#Failed as no error found as expected
expect_error(DLpencoef(y1), "'x' must be numeric")
#Passed as error found as expected
expect_error(DLpencoef(y2), "'x' must be numeric")
#Passed as error found as expected

#### ar ####
expect_error(ar(y), "'x' must be numeric")
#Failed as no error found as expected
expect_error(ar(y1), "'x' must be numeric")
#Passed as error found as expected
expect_error(ar(y2), "'x' must be numeric")
#Passed as error found as expected

#### ar.penyw ####
expect_error(ar.penyw(y), "'x' must be numeric")
#Failed as no error found as expected
expect_error(ar.penyw(y1), "'x' must be numeric")
#Passed as error found as expected
expect_error(ar.penyw(y2), "'x' must be numeric")
#Passed as error found as expected


## lh
expect_error(acf(y, lh = NULL), "lh must either be NULL, length 1, length order.max, ncol(x), or a matrix with dimension order.max x nser.")
expect_error(pacf(y, lh = NULL), "lh must either be NULL, length 1, length order.max, ncol(x), or a matrix with dimension order.max x nser.")
expect_error(auto.acf(y, lh = NULL), "lh must either be NULL, length 1, length order.max, ncol(x), or a matrix with dimension order.max x nser.")
expect_error(invertpacf(y, lh = NULL), "lh must either be NULL, length 1, length order.max, ncol(x), or a matrix with dimension order.max x nser.")
expect_error(DLpencoef(y, lh = NULL), "lh must either be NULL, length 1, length order.max, ncol(x), or a matrix with dimension order.max x nser.")
expect_error(ar(y, lh = NULL), "lh must either be NULL, length 1, length order.max, ncol(x), or a matrix with dimension order.max x nser.")
expect_error(ar.penyw(y, lh = NULL), "lh must either be NULL, length 1, length order.max, ncol(x), or a matrix with dimension order.max x nser.")


## testing that lambda can only be numeric
### acf ###
expect_error(acf(y, lambda = 0), "lambda must be numeric")
#Failed as no errors found as expected
expect_error(acf(y, lambda = TRUE), "lambda must be numeric")
#Passed as errors found as expected
expect_error(acf(y, lambda = NULL), "lambda must be numeric")
#Failed as no errors found as expected

### pacf ###
expect_error(pacf(y, lambda = 0), "lambda must be numeric")
#Failed as no errors found as expected
######failed for na not found so not correct yet
expect_error(pacf(y, lambda = TRUE), "lambda must be numeric")
#Passed as errors found as expected
expect_error(pacf(y, lambda = NULL), "lambda must be numeric")
#Failed as no errors found as expected



## testing that lambda can only be positive
### acf ###
expect_error(acf(y, lambda = 0), "lambda must be positive")
#Failed as no errors found as expected
expect_error(acf(y, lambda = -1), "lambda must be positive")
#Passed as errors found as expected
expect_error(acf(y, lambda = NULL), "lambda must be positive")
#Failed as no errors found as expected

### pacf ###
expect_error(pacf(y, lambda = 0), "lambda must be positive")
#Failed as no errors found as expected
expect_error(pacf(y, lambda = -1), "lambda must be positive")
#Passed as errors found as expected
expect_error(pacf(y, lambda = NULL), "lambda must be positive")
#Failed as no errors found as expected



## testing that lambda can only be a matrix
### acf ###
######## na not found so all fail which is wrong
expect_error(acf(y, lambda = 0), "lambda must be a matrix")
#Failed as no errors found as expected
expect_error(acf(y, lambda = NULL), "lambda must be a matrix")
#Failed as no errors found as expected

### pacf ###
expect_error(pacf(y, lambda = 0), "lambda must be a matrix")
#Failed as no errors found as expected
expect_error(pacf(y, lambda = NULL), "lambda must be a matrix")
#Failed as no errors found as expected



## testing that lambda must either be NULL, length 1, length lag.max, ncol(x), or a matrix with dimension lag.max x nser.
### acf ###
expect_error(acf(y, lambda = 0), "lambda must either be NULL, length 1, length lag.max, ncol(x), or a matrix with dimension lag.max x nser.")
#Failed as no errors found as expected
expect_error(acf(y, lambda = ncol(y)), "target must either be NULL, length 1, length lag.max, ncol(x), or a matrix with dimension lag.max x nser.")
#Failed as no errors found as expected
expect_error(acf(y, lambda = NULL), "lambda must either be NULL, length 1, length lag.max, ncol(x), or a matrix with dimension lag.max x nser.")
#Failed as no errors found as expected

### pacf ###
expect_error(pacf(y, lambda = 0), "lambda must either be NULL, length 1, length lag.max, ncol(x), or a matrix with dimension lag.max x nser.")
#Failed as no errors found as expected
expect_error(acf(y, lambda = ncol(y)), "target must either be NULL, length 1, length lag.max, ncol(x), or a matrix with dimension lag.max x nser.")
#Failed as no errors found as expected
expect_error(pacf(y, lambda = NULL), "lambda must either be NULL, length 1, length lag.max, ncol(x), or a matrix with dimension lag.max x nser.")
#Failed as no errors found as expected



## testing that target can only be numeric
### acf ###
expect_error(acf(y, target = 0), "target must be numeric")
#Failed as no errors found as expected
expect_error(acf(y, target = TRUE), "target must be numeric")
#Passed as errors found as expected
expect_error(acf(y, target = NULL), "target must be numeric")
#Failed as no errors found as expected

### pacf ###
expect_error(pacf(y, target = 0), "target must be numeric")
#Failed as no errors found as expected
expect_error(pacf(y, target = TRUE), "target must be numeric")
#Passed as errors found as expected
expect_error(pacf(y, target = NULL), "target must be numeric")
#Failed as no errors found as expected



## testing that target can only be between 1 and -1
### acf ###
expect_error(acf(y, target = 0), "target must be between 1 and -1")
#Failed as no errors found as expected
expect_error(acf(y, target = -2), "target must be between 1 and -1")
#Passed as errors found as expected
expect_error(acf(y, target = 2), "target must be between 1 and -1")
#Passed as errors found as expected
expect_error(acf(y, target = NULL), "target must be between 1 and -1")
#Failed as no errors found as expected

### pacf ###
expect_error(pacf(y, target = 0), "target must be between 1 and -1")
#Failed as no errors found as expected
expect_error(acf(y, target = -2), "target must be between 1 and -1")
#Passed as errors found as expected
expect_error(acf(y, target = 2), "target must be between 1 and -1")
#Passed as errors found as expected
expect_error(pacf(y, target = NULL), "target must be between 1 and -1")
#Failed as no errors found as expected



## testing that target can only be a matrix
### acf ###
expect_error(acf(y, target = 0), "target must be a matrix")
#Failed as no errors found as expected
expect_error(acf(y, target = NULL), "target must be a matrix")
#Failed as no errors found as expected

### pacf ###
expect_error(pacf(y, target = 0), "target must be a matrix")
#Failed as no errors found as expected
expect_error(pacf(y, target = NULL), "target must be a matrix")
#Failed as no errors found as expected



## testing that target must either be NULL, length 1, length lag.max, ncol(x), or a matrix with dimension lag.max x nser.
### acf ###
expect_error(acf(y, target = 0), "target must either be NULL, length 1, length lag.max, ncol(x), or a matrix with dimension lag.max x nser.")
#Failed as no errors found as expected
expect_error(acf(y, target = ncol(y)), "target must either be NULL, length 1, length lag.max, ncol(x), or a matrix with dimension lag.max x nser.")
#Failed as no errors found as expected
expect_error(acf(y, target = NULL), "target must either be NULL, length 1, length lag.max, ncol(x), or a matrix with dimension lag.max x nser.")
#Failed as no errors found as expected

### pacf ###
expect_error(pacf(y, target = 0), "target must either be NULL, length 1, length lag.max, ncol(x), or a matrix with dimension lag.max x nser.")
#Failed as no errors found as expected
expect_error(pacf(y, target = ncol(y)), "target must either be NULL, length 1, length lag.max, ncol(x), or a matrix with dimension lag.max x nser.")
#Failed as no errors found as expected
expect_error(pacf(y, target = NULL), "target must either be NULL, length 1, length lag.max, ncol(x), or a matrix with dimension lag.max x nser.")
#Failed as no errors found as expected