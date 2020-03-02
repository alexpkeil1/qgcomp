cat("# zi test\n")
library(qgcomp)
set.seed(1231)
n=100
dat <- data.frame(y=rbinom(n, 1, 0.5)*rpois(n, 1.2), x1=runif(n), x2=runif(n), z=runif(n))

# should not cause error
ee = qgcomp.zi.noboot(f=y ~ z + x1 + x2 | z, expnms = c('x1', 'x2'), data=dat, q=2, dist="negbin")
qgcomp.zi.noboot(f=y ~ z + x1 + x2 | x1 + x2 + z, expnms = c('x1', 'x2'), data=dat, q=2, dist="negbin")

qgcomp.zi.boot(f=y ~ z + x1 + x2 | z, expnms = c('x1', 'x2'), data=dat, q=4, B = 10, MCsize = 100, 
               dist="negbin")
#'\donttest{
#' qgcomp.zi.boot(f=y ~ z + x1 + x2 | z, expnms = c('x1', 'x2'), data=dat, q=4, B = 10, MCsize = 1000, 
#'                dist="negbin", msmcontrol = zimsm.fit.control(predmethod="catprobs"))
#' ff = qgcomp.zi.boot(f=y ~ z + x1 + x2 | z, expnms = c('x1', 'x2'), data=dat, q=4, B = 10, MCsize = 100, 
#'                dist="geometric")
#' plot(ff)
#'}              


# should cause error
res <- try(qgcomp.zi.noboot(f=y ~ z + x1 + x2 | x1 + z, expnms = c('x1', 'x2'), data=dat, q=2, dist="negbin"))
stopifnot(class(res)=="try-error")

cat("done")
