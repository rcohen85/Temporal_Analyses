# this was turned into an RMD!

library(splines)
library(splines2)
library(mgcv)

# M-splines ---------------------------------------------------------------

### example given in the reference paper by Ramsay (1988)
x <- seq.int(0, 1, 0.01)
knots <- c(0.3, 0.5, 0.6)
msMat <- mSpline(x, knots = knots, degree = 2, intercept = TRUE)
matplot(x, msMat, type = "l", ylab = "y", ylim = c(0,5), main = "m-spline")
abline(v = knots, lty = 2, col = "gray")

# Using real data, but (mostly) the same settings:
knots = quantile(Jday, probs=c(0.333,0.666))
Jday_mat1 = mSpline(Jday,knots=knots,
                     degree = 2, intercept = TRUE)
matplot(Jday, Jday_mat1, type = "l", ylab = "y", main = "m-spline, real JDay data")
abline(v = knots, lty = 2, col = "gray")

# Using real data and my own settings:
knots = quantile(Jday, probs=c(0.333,0.666))
Jday_mat2 = mSpline(Jday,knots=knots,
                     Boundary.knots=c(1,365),
                     periodic=T)
matplot(Jday, Jday_mat2, type = "l", ylab = "y", main = "m-spline, real JDay data, RC settings")
abline(v = knots, lty = 2, col = "gray")

# B-splines ------------------------------------------------------

# Same as above but with a b-spline
x <- seq.int(0, 1, 0.01)
knots <- c(0.3, 0.5, 0.6)
msMat <- bs(x, knots = knots, degree = 2, intercept = TRUE)
matplot(x, msMat, type = "l", ylab = "y", ylim = c(0,5), main = "b-spline")
abline(v = knots, lty = 2, col = "gray")

# No option for cyclic variables in bs()
knots = quantile(Jday, probs=c(0.333,0.666))
Jday_mat3 <- bs(Jday, knots = knots, degree = 2, intercept = TRUE)
matplot(Jday, Jday_mat2, type = "l", ylab = "y", main = "b-spline, real JDay data, non-cyclic, degree = 2")
abline(v = knots, lty = 2, col = "gray")

Jday_mat3 <- bs(Jday, knots = knots, degree = 3, intercept = TRUE)
matplot(Jday, Jday_mat2, type = "l", ylab = "y", main = "b-spline, real JDay data, non-cyclic, degree = 3")
abline(v = knots, lty = 2, col = "gray")

# Cubic Splines -----------------------------------------------------------




# Cyclic Variance-Covariance Splines ----------------------------------------------
# from Benjamins '17 harbour porpoise paper
knots = seq(1,365,length=6)
JdayBasis = gam(thisSite~s(Jday, bs="cc", k=6), fit=FALSE, family=poisson, 
                knots = list(knots))$X[,2:5]
matplot(Jday, JdayBasis, type = "l", ylab = "y", ylim = c(0,5), main = "Variance-Covariance Smooth")
abline(v = knots, lty = 2, col = "gray")

