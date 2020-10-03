cat("# factors test\n")
library("qgcomp")
library("survival")
library('pscl')
# are results at a given seed numerically stable across versions?
set.seed(50)
N=50
dat <- data.frame(time=(tmg <- pmin(.1,rweibull(N, 10, 0.1))), 
                d=1.0*(tmg<0.1), x1=runif(N), x2=runif(N), z=runif(N), 
                z2 = as.factor(sample(c(1,2,3), size=N, replace=TRUE)))
expnms=paste0("x", 1:2)



##### binomial
set.seed(123123)
f0 = d ~ x1 + x2 + z2
obj0a <- qgcomp.noboot(f0, expnms = expnms, data = dat, family=binomial())
print(obj0a)
pointwisebound.noboot(obj0a)

##### survival
f1 = survival::Surv(time, d)~ x1 + x2 + z2
obj0b <- qgcomp.cox.noboot(f1, expnms = expnms, data = dat)
print(obj0b)

res = try(pointwisebound.noboot(obj0b), silent=TRUE)
stopifnot(class(res)=="try-error")

##### zi
f2 = d ~ x1 + x2 + z2 | x1 + x2 + z2
pp = pscl::zeroinfl(formula = f2, data = dat)
obj0c <- qgcomp.zi.noboot(f2, expnms = expnms, data = dat)
res = try(pointwisebound.noboot(obj0c), silent=TRUE)
stopifnot(class(res)=="try-error")

