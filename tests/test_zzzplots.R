cat("# plot test\n")
library(qgcomp)
set.seed(112312)
n=100

#base
# binomial
   dat <- data.frame(y=rbinom(n, 1, 0.5), x1=runif(n), x2=runif(n), z=runif(n))
   ee = qgcomp.noboot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=7, family=binomial())
   plot(ee)
   ff = qgcomp.boot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=7, B=20, family=binomial())
   plot(ff)

# gaussian
   dat <- data.frame(y=rnorm(n), x1=runif(n), x2=runif(n), z=runif(n))
   ee = qgcomp.noboot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=7, family=gaussian())
   plot(ee)
   ff = qgcomp.boot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=7, B=20, family=gaussian())
   plot(ff)
   
# poisson
   dat <- data.frame(y=rpois(n, 1.2), x1=runif(n), x2=runif(n), z=runif(n))
   ee = qgcomp.noboot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=7, family=poisson())
   plot(ee)
   ff = qgcomp.boot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=7, B=20, family=poisson())
   qgcomp::modelbound.boot(ff)
   plot(ff)
      
#cox
   dat <- data.frame(stop=(tmg <- pmin(.1,rweibull(n, 10, 0.1))), 
                     start = pmax(0, tmg-runif(n)),
                     d=1.0*(tmg<0.1), x1=runif(n), x2=runif(n), z=runif(n))
   expnms=paste0("x", 1:2)
   f = survival::Surv(start,stop, d)~x1 + x2
   suppressWarnings(ee <- qgcomp.cox.noboot(f, expnms = expnms, data = dat))
   plot(ee)
   suppressWarnings(ff <- qgcomp.cox.boot(f, expnms = expnms, data = dat, B=12, MCsize=1000))
   plot(ff)
   
# zi
  dat <- data.frame(y=rbinom(n, 1, 0.5)*rpois(n, 1.2), x1=runif(n), x2=runif(n), z=runif(n))
  ee = qgcomp.zi.noboot(f=y ~ z + x1 + x2 | z, expnms = c('x1', 'x2'), data=dat, q=7, dist="negbin")
  plot(ee)
  ff = qgcomp.zi.boot(f=y ~ z + x1 + x2 | z, expnms = c('x1', 'x2'), data=dat, q=7, B=10, MCsize=1000, dist="negbin")
  modelbound.boot(ff, pwonly=TRUE)
  plot(ff)

  
cat("done")
