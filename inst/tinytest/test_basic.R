library(tinytest)

set.seed(1234)
data=cbind(arima.sim(model=list(ar=0.7),n=500),arima.sim(model=list(ar=0.3),n=500))

acfdata=acf(data)
expect_equal_to_reference(acfdata,'./expected/acfdata.Rdata')

pacfdata=pacf(data)
expect_equal_to_reference(pacfdata,'./expected/pacfdata.Rdata')

ardata=ar(data)
expect_equal_to_reference(ardata,'./expected/ardata.Rdata')

# test a single lag, one series
set.seed(5923)
datama=arima.sim(model=list(ma=0.99),n=100)
acfdatamalag1=acf(datama,lag.max=1)
expect_equal_to_reference(acfdatamalag1,'./expected/acfdatamalag1.Rdata')

# test a NND
acfdatama=acf(datama)
expect_equal_to_reference(acfdatama,'./expected/acfdatama.Rdata')

# test AR order larger than 3 warning
acfdatamawarn=acf(datama,estimate="invertpacf")
expect_equal_to_reference(acfdatamawarn,'./expected/acfdatamawarn.Rdata')

pacfdatama=pacf(datama)
expect_equal_to_reference(pacfdatama,'./expected/pacfdatama.Rdata')



# Colin's example
set.seed(56334)
x=arima.sim(n=1000,model=list(ar=.9))
original=acf(x,penalized=FALSE)
expect_equal_to_reference(original,'./expected/original.Rdata')

pen=acf(x)
expect_equal_to_reference(pen,'./expected/pen.Rdata')

invert=acf(x,estimate="invertpacf")
expect_equal_to_reference(invert,'./expected/invert.Rdata')
