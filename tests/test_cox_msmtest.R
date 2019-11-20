cat("# cox msm test\n")
library(qgcomp)
# does the simulation version of a MS cox model yield the same result as
# the sum-of-random-variables version when there are no non-exposure covariates
# THIS TEST TAKES A LONG TIME!
#set.seed(50)
#N=500
#dat <- data.frame(time=(tmg <- pmin(.1,rweibull(N, 10, 0.1))), 
#                d=1.0*(tmg<0.1), x1=runif(N), x2=runif(N), z=runif(N))
#expnms=paste0("x", 1:2)
#f = survival::Surv(time, d)~x1 + x2
#f2 = survival::Surv(time, d)~x1 + x2 + z
#f3 = survival::Surv(time, d)~x1 + x2 + I(z-mean(z))
#survival::coxph(f, data = dat)
#survival::coxph(f2, data = dat)
#(obj <- qgcomp.cox.noboot(f, expnms = expnms, data = dat))
#(obj2 <- qgcomp.cox.boot(f, expnms = expnms, data = dat, B=1000, MCiter=20000))
#(obj3 <- qgcomp.cox.boot(f2, expnms = expnms, data = dat, B=1000, MCiter=20000, parallel=TRUE))
#
#stopifnot(all.equal(obj$psi, obj2$psi, tolerance = .01))
#stopifnot(all.equal(obj$var.psi, obj2$var.psi, tolerance = .1))
#
##(r00 <- qgcomp.cox.noboot(f2, expnms = expnms, data = dat))
##system.time(r0 <- qgcomp.cox.boot(f2, expnms = expnms, data = dat, 
##                            B=16, MCiter=20000, parallel=FALSE))
###  Expected time to finish: 1.54 minutes 
###  user  system elapsed 
###  90.438  11.711 102.804 
##system.time(r1 <- qgcomp.cox.boot(f2, expnms = expnms, data = dat, 
##                            B=16, MCiter=20000, parallel=TRUE))
### fairly high overhead vs. mclapply (but works on windows)
### Expected time to finish: 3.75 minutes 
### user  system elapsed 
### 10.922   1.846  51.252 
#
#
#
##(obj4 <- qgcomp.cox.boot(survival::Surv(time, d)~factor(x1) + splines::bs(x2) + z, 
##                         expnms = expnms, data = dat, 
##                         B=10, MCiter=20000, parallel=FALSE, degree=2))
##
##lapply(1:5, mean)
##parallel::mclapply(1:5, mean)
### windows friendly version
##future.apply::future_sapply(1:5, mean)
#
## LATE ENTRY

set.seed(50)
N=50
dat <- data.frame(stop=(tmg <- pmin(.1,rweibull(N, 10, 0.1))), 
                  start = pmax(0, tmg-runif(N)),
                  d=1.0*(tmg<0.1), x1=runif(N), x2=runif(N), z=runif(N))
expnms=paste0("x", 1:2)
f = survival::Surv(start,stop, d)~x1 + x2
(obj <- qgcomp.cox.boot(f, expnms = expnms, data = dat, B=1, MCsize=5000))
plot(obj)


ymat = obj$fit$y
tval = grep("stop|time",colnames(ymat) , value=TRUE)
stop = as.numeric(ymat[,tval])
times = sort(-sort(-unique(stop))[-1])
datatimes = with(dat, sort(unique(stop*d))[-1])

stopifnot(all.equal(times, datatimes))

