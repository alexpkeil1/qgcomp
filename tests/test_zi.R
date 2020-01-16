cat("# zi test\n")
library(qgcomp)
n=100
dat <- data.frame(y=rbinom(n, 1, 0.5)*rpois(n, 1.2), x1=runif(n), x2=runif(n), z=runif(n))

# should not cause error
qgcomp.zi.noboot(f=y ~ z + x1 + x2 | z, expnms = c('x1', 'x2'), data=dat, q=2, dist="negbin")
qgcomp.zi.noboot(f=y ~ z + x1 + x2 | x1 + x2 + z, expnms = c('x1', 'x2'), data=dat, q=2, dist="negbin")

# should cause error
res <- try(qgcomp.zi.noboot(f=y ~ z + x1 + x2 | x1 + z, expnms = c('x1', 'x2'), data=dat, q=2, dist="negbin"))
stopifnot(class(res)=="try-error")

cat("done")
