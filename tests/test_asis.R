cat("# asis test\n")
library("qgcomp")
# does the checknames correctly intuit which are linear/non-linear models for common specifications
set.seed(50)
N=50
dat <- data.frame(time=(tmg <- pmin(.1,rweibull(N, 10, 0.1))), 
                d=1.0*(tmg<0.1), x1=runif(N), x2=runif(N), z=runif(N))
expnms=paste0("x", 1:2)


##### non additive

f0 = survival::Surv(time, d)~ x1 + x2
f1 = survival::Surv(time, d)~ x1*x2
f2 = survival::Surv(time, d)~ x1 + x2 + I(x1*x2)
f3 = survival::Surv(time, d)~x1 + x2 + x1:x2
set.seed(12321)
#res = try(qgcomp(f0, expnms = expnms, data = dat, B=2, MCsize=50))
#stopifnot(class(res)=="try-error")
res = try(obj0 <- qgcomp(f0, expnms = expnms, data = dat))
stopifnot(inherits(res, "qgcompfit"))
set.seed(12321)
res = try(obj1 <- qgcomp(f1, expnms = expnms, data = dat, B=2, MCsize=50))
stopifnot(inherits(res, "qgcompfit"))
set.seed(12321)
res = try(obj2 <- qgcomp(f2, expnms = expnms, data = dat, B=2, MCsize=50))
stopifnot(inherits(res, "qgcompfit"))
set.seed(12321)
res = try(obj3 <- qgcomp(f3, expnms = expnms, data = dat, B=2, MCsize=50))
stopifnot(inherits(res, "qgcompfit"))
stopifnot(all.equal(
  obj0$fit$coefficients, 
  obj1$fit$coefficients, check.names=FALSE
)=="Numeric: lengths (2, 3) differ")

stopifnot(all.equal(
  obj1$fit$coefficients, 
  obj2$fit$coefficients, 
  obj3$fit$coefficients, check.names=FALSE))
stopifnot(all.equal(
  obj1$msmfit$coefficients, 
  obj2$msmfit$coefficients, 
  obj3$msmfit$coefficients, check.names=FALSE))

f0b = d ~ x1 + x2
f4 = d ~ x1*x2
f5 = d ~ x1 + x2 + I(x1*x2) 
f6 = d ~ x1 + x2 + x1:x2
res = try(obj0b <- qgcomp(f0b, expnms = expnms, data = dat, B=2, family=binomial()))
stopifnot(inherits(res,"qgcompfit"))
res = try(obj4 <- qgcomp(f4, expnms = expnms, data = dat, B=2, family=binomial()))
stopifnot(inherits(res,"qgcompfit"))
res = try(obj5 <- qgcomp(f5, expnms = expnms, data = dat, B=2, family=binomial()))
stopifnot(inherits(res,"qgcompfit"))
res = try(obj6 <- qgcomp(f6, expnms = expnms, data = dat, B=2, family=binomial()))
stopifnot(inherits(res,"qgcompfit"))

stopifnot(all.equal(
  obj0b$fit$coefficients, 
  obj4$fit$coefficients, check.names=FALSE
)=="Numeric: lengths (3, 4) differ")

stopifnot(all.equal(
  obj4$fit$coefficients, 
  obj6$fit$coefficients, check.names=FALSE))
stopifnot(all.equal(
  obj4$msmfit$coefficients, 
  obj6$msmfit$coefficients, check.names=FALSE))
stopifnot(all.equal(
  obj4$fit$coefficients, 
  obj5$fit$coefficients, check.names=FALSE))
stopifnot(all.equal(
  obj4$msmfit$coefficients, 
  obj5$msmfit$coefficients, check.names=FALSE))


##### non linear


# behavior changes when considering self interactions!
f0c = survival::Surv(time, d)~ x1 + x1:x1 
f1c = survival::Surv(time, d)~ x1^2
f3c = survival::Surv(time, d)~ x1 + I(x1*x1)
f4c = survival::Surv(time, d)~ x1 + I(x1^2)
stopifnot(all.equal(
  survival::coxph(f0c, data=dat)$coefficients,
  survival::coxph(f1c, data=dat)$coefficients, check.names=FALSE))
stopifnot(all.equal(
  survival::coxph(f3c, data=dat)$coefficients,
  survival::coxph(f4c, data=dat)$coefficients, check.names=FALSE))

stopifnot(all.equal(
  survival::coxph(f0c, data=dat)$coefficients,
  survival::coxph(f3c, data=dat)$coefficients, check.names=FALSE)=="Numeric: lengths (1, 2) differ")



f0 = d~ x1 + x1:x1 
f1 = d~ x1^2
f3 = d~ x1 + I(x1*x1)
f4 = d~ x1 + I(x1^2)
glm(f0, data=dat, family=binomial())$coefficients
glm(f1, data=dat, family=binomial())$coefficients
glm(f3, data=dat, family=binomial())$coefficients
glm(f4, data=dat, family=binomial())$coefficients

stopifnot(all.equal(
  glm(f0, data=dat, family=binomial())$coefficients,
  glm(f1, data=dat, family=binomial())$coefficients, check.names=FALSE))
stopifnot(all.equal(
  glm(f3, data=dat, family=binomial())$coefficients,
  glm(f4, data=dat, family=binomial())$coefficients, check.names=FALSE))

stopifnot(all.equal(
  glm(f0, data=dat, family=binomial())$coefficients,
  glm(f3, data=dat, family=binomial())$coefficients, check.names=FALSE)=="Numeric: lengths (2, 3) differ")


res = try(obj0 <- qgcomp(f0, expnms = expnms, 
                         data = dat, B=2, family=binomial())) # defaults to RR
stopifnot(inherits(res,"qgcompfit"))
#res = try(obj0 <- qgcomp(f0, expnms = expnms, 
#         data = dat, B=2, family=gaussian())) # gives error due to B

res = try(obj1 <- qgcomp(f1, expnms = expnms, data = dat, B=2, family=binomial()))
#res = try(obj1 <- qgcomp(f1, expnms = expnms, 
#         data = dat, B=2, family=gaussian())) # gives error due to B


res = try(obj3 <- qgcomp(f3, expnms = expnms, 
                         data = dat, B=2, family=binomial())) # defaults to RR
stopifnot(inherits(res,"qgcompfit"))
#res = try(obj0 <- qgcomp(f0, expnms = expnms, 
#         data = dat, B=2, family=gaussian())) # gives error due to B

res = try(obj4 <- qgcomp(f4, expnms = expnms, data = dat, B=2, family=binomial()))
#res = try(obj1 <- qgcomp(f1, expnms = expnms, 
#         data = dat, B=2, family=gaussian())) # gives error due to B

stopifnot(inherits(res,"qgcompfit"))
stopifnot(all.equal(
  obj0$fit$coefficients, 
  obj1$fit$coefficients, check.names=FALSE))
stopifnot(all.equal(
  obj3$fit$coefficients, 
  obj4$fit$coefficients, check.names=FALSE))

stopifnot(all.equal(
  obj0$fit$coefficients,
  obj3$fit$coefficients, check.names=FALSE)=="Numeric: lengths (2, 3) differ")

# should result in different MSMs
stopifnot(substr(all.equal(
  obj0$msmfit$coefficients,
  obj3$msmfit$coefficients, check.names=FALSE), 1, 10)=="Mean relat")

cat("done")
