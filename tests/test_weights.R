cat("# weights test\n")
library("qgcomp")
# are results at a given seed numerically stable across versions?
set.seed(50)
N=100
dat <- data.frame(time=(tmg <- pmin(.1,rweibull(N, 10, 0.1))), 
                  d=1.0*(tmg<0.1), x1=runif(N), x2=runif(N), z=runif(N))
dat$wt = runif(N)*40
dat$wt = dat$wt/mean(dat$wt)
dat$wt2 = 1
expnms=paste0("x", 1:2)



##### binomial
set.seed(123123)
f0 = d ~ x1 + x2
obj0a <- qgcomp.noboot(f0, expnms = expnms, data = dat, family=binomial())
obj0b <- qgcomp.noboot(f0, expnms = expnms, data = dat, weight=wt2, family=binomial())
obj0c <- qgcomp.noboot(f0, expnms = expnms, data = dat, weight=wt, family=binomial())

stopifnot(all.equal(
  coef(obj0a),
  coef(obj0b),
  check.names=FALSE, tolerance = 1e-4))

##### survival
#set.seed(123123)
#f1 = survival::Surv(tmg, d) ~ x1 + x2
#obj1a <- qgcomp(f1, expnms = expnms, data = dat)
#obj1b <- qgcomp.cox.boot(f1, expnms = expnms, B=2, data = dat)
#
#stopifnot(all.equal(
#  coef(obj1a),
#  c(0.02461938), 
#  check.names=FALSE, tolerance = 1e-4))
#stopifnot(all.equal(
#  coef(obj1b),
#  c(0.02654175), 
#  check.names=FALSE, tolerance = 1e-4))
#
#
###### zi
#set.seed(123123)
#f2 = d ~ x1 + x2 | x1+x2
#obj2a <- qgcomp(f2, expnms = expnms, data = dat, family=poisson())
#obj2b <- qgcomp.zi.boot(f2, expnms = expnms, data = dat, dist="poisson", B = 2)
#
#stopifnot(all.equal(
#  coef(obj2a),
#  list(count=c( -0.42823824,-0.05519007 ),
#       zero=c(9.057221, -40.356026)),
#  check.names=FALSE, tolerance = 1e-4))
#stopifnot(all.equal(
#  coef(obj2b),
#  list(count=c(  -0.43554489,-0.05279048 ),
#       zero=c(7.369426 , -18.411273)),
#  check.names=FALSE, tolerance = 1e-4))
#
