cat("# zi test\n")
library(qgcomp)
set.seed(1231)
n=300
dat <- data.frame(y=rbinom(n, 1, 0.5)*rpois(n, 1.2), x1=runif(n), x2=runif(n), z=runif(n))

# should not cause error
qgcomp.zi.noboot(f=y ~ z + x1 + x2 | z, expnms = c('x1', 'x2'), data=dat, q=2, dist="negbin")
qgcomp.zi.noboot(f=y ~ z + x1 + x2 | x1 + x2 + z, expnms = c('x1', 'x2'), data=dat, q=2, dist="negbin")
qgcomp.hurdle.noboot(f=y ~ z + x1 + x2 | x1 + x2 + z, expnms = c('x1', 'x2'), data=dat, q=2, dist="negbin")

qgcomp.zi.boot(f=y ~ z + x1 + x2 | x1 + x2 + z, B=1, expnms = c('x1', 'x2'), msmcontrol =zimsm.fit.control(predmethod="catprobs") , data=dat, q=2, dist="negbin")
qgcomp.hurdle.boot(f=y ~ z + x1 + x2 | x1 + x2 + z, B=1, expnms = c('x1', 'x2'), msmcontrol =zimsm.fit.control(predmethod="components") , data=dat, q=2, dist="negbin")



qgcomp.zi.boot(f=y ~ z + x1 + x2 | z, expnms = c('x1', 'x2'), data=dat, q=4, B = 10, MCsize = 100, 
               dist="poisson")

# qgcomp.hurdle.boot is iffy when mixture is not in zero model
#rr = qgcomp.hurdle.boot(f=y ~ z + x1 + x2 | x1 + x2 + z, expnms = c('x1', 'x2'), data=dat, q=2, B = 10, MCsize = 1000)
#summary(rr)
#'\dontrun{
#' qgcomp.zi.boot(f=y ~ z + x1 + x2 | z, expnms = c('x1', 'x2'), data=dat, q=4, B = 10, MCsize = 1000, 
#'                dist="negbin", msmcontrol = zimsm.fit.control(predmethod="catprobs"))
#' ff = qgcomp.zi.boot(f=y ~ z + x1 + x2 | z, expnms = c('x1', 'x2'), data=dat, q=4, B = 10, MCsize = 100, 
#'                dist="geometric")
#' plot(ff)
#'}              


# should cause error
res <- try(qgcomp.zi.noboot(f=y ~ z + x1 + x2 | x1 + z, expnms = c('x1', 'x2'), data=dat, q=2, dist="negbin"), silent = TRUE)
stopifnot(class(res)=="try-error")

cat("done")
