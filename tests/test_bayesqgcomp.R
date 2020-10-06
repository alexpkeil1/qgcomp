cat("# bayesqgcomp test\n")
#devtools::install_github("alexpkeil1/qgcomp@dev")
library("qgcomp")
data("metals", package="qgcomp")

Xnm <- c(
  'arsenic','barium','cadmium','calcium','chromium','copper',
  'iron','lead','magnesium','manganese','mercury','selenium','silver',
  'sodium','zinc'
)

# continuous outcome, example with perfect collinearity
metals$cadmium2 = metals$cadmium
Xnm2 = c(Xnm, "cadmium2")
res = try(qc.fit <- qgcomp.noboot(y~.,dat=metals[,c(Xnm2, 'y')], family=gaussian()), silent=TRUE)
stopifnot(class(res)=="try-error")
# error

res = try(qc.fit <- qgcomp.noboot(y~.,dat=metals[,c(Xnm2, 'y')], family=gaussian(), bayes=TRUE))
stopifnot(class(res)=="qgcompfit")

# compare to results with colinear exposure removed
res = try(qc.fit2 <- qgcomp.noboot(y~.,dat=metals[,c(Xnm, 'y')], family=gaussian()))
stopifnot(class(res)=="qgcompfit")

#hitting code coverage just to check
dat = qgcomp:::.dgm_quantized(N=100)
dat$y = as.numeric(dat$y>median(dat$y))

qgcomp(y~.,dat=dat, expnms=c("x1", "x2"), family="binomial", rr=FALSE)

res = try(qgcomp(y~x1+x2,dat=dat, expnms=c("x1", "x2"), family=NULL, rr=FALSE), silent=TRUE)
stopifnot(inherits(res, "try-error"))

qgcomp(y~x1+x2 | x1+x2,dat=dat, expnms=c("x1", "x2"), family="poisson")

res = try(qgcomp(y~x1+x2 | x1+x2 +I(x2^2), B=2,dat=dat, expnms=c("x1", "x2"), bayes=TRUE, family="poisson"), silent=TRUE)
stopifnot(inherits(res, "try-error"))


ft = qgcomp(y~x1+x2 + I(x2^2), B=2,dat=dat, expnms=c("x1", "x2"), bayes=TRUE, family="poisson")
pp = msm.predict(ft)
ft = qgcomp(y~x1+x2,dat=dat, expnms=c("x1", "x2"), bayes=TRUE, family="poisson")
res = try(msm.predict(ft), silent = TRUE)
stopifnot(inherits(res, "try-error"))


ft = qgcomp(y~x1+x2 + I(x2^2), B=2, q=NULL,dat=dat, expnms=c("x1", "x2"), bayes=TRUE, family=binomial())
summary(ft)
qc.fit2 <- qgcomp.boot(y~.,expnms = Xnm, dat=metals[,c(Xnm, 'y')],B=2, q=NULL, family=gaussian())


## slight shrinkage with heavily non-linear models that don't converge
## non-Bayesian
#results5 = qgcomp(y~. + .^2 + .^3 + arsenic*cadmium,
#                  expnms=Xnm,
#                  metals[,c(Xnm, 'y')], family=gaussian(), q=10, B=10, 
#                  seed=125, degree=3)
#
#print(results5)
#plot(results5)
#results5$bootsamps # samples are wack
#
## Bayesian (takes a few mins)
#results5bayes = qgcomp.boot(y~. + .^2 + .^3 + arsenic*cadmium,
#                  expnms=Xnm,
#                  metals[,c(Xnm, 'y')], family=gaussian(), q=10, B=10, 
#                  seed=125, degree=3, 
#                  bayes=TRUE)
#
## note p-values are not computed due to negative
## degrees of freedom (can fix by using a z-statistic instead of a t-statistic)
#print(results5bayes)
#plot(results5bayes)
#results5bayes$bootsamps # samples are reasonable


cat("done")
