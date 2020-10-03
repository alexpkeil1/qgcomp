cat('testing utility functions')


qgcomp:::construction("msg", "test")

qgcomp:::cox()

qgcomp:::zi()


# zero inflated models
set.seed(50)
n=100
dat <- data.frame(y=rbinom(n, 1, 0.5)*rpois(n, 1.2), x1=runif(n), x2=runif(n), z=runif(n))

# poisson count model, mixture in both portions
qgcomp::qgcomp.zi.noboot(f=y ~ z + x1 + x2 | x1 + x2, expnms = c('x1', 'x2'), 
    data=dat, q=2, dist="poisson")

qgcomp::qgcomp.zi.noboot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), 
   data=dat, q=2, dist="negbin") # equivalent


ftz = qgcomp::qgcomp.zi.boot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), 
                            data=dat, q=2, dist="negbin", B=2, parallel=TRUE) # equivalent
fth = qgcomp::qgcomp.hurdle.boot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), 
                            data=dat, q=2, dist="poisson", B=2, parallel=TRUE) # equivalent


res = try(qgcomp:::glance.qgcompfit(ftz), silent=TRUE)
stopifnot(inherits(res,"try-error"))
res = try(qgcomp:::glance.qgcompfit(fth), silent=TRUE)
stopifnot(inherits(res,"try-error"))

ft = qgcomp::qgcomp.noboot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), 
                           data=dat, q=2, family="gaussian") # equivalent
qgcomp:::glance.qgcompfit(ft)
ft = qgcomp::qgcomp.noboot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), 
                           data=dat, q=2, family="poisson") # equivalent

qgcomp:::tidy.qgcompfit(ft)
