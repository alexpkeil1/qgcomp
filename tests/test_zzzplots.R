cat("# plot test\n")
library(qgcomp)
set.seed(112312)
n=100

#base
# binomial
   dat <- data.frame(y=rbinom(n, 1, 0.5), x1=runif(n), x2=runif(n), z=runif(n))
   ee = qgcomp.noboot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=7, family=binomial())
   plot(ee)
   plot(ee)
   ff = qgcomp.boot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=7, B=5, family=binomial())
   plot(ff)
   pointwisebound.boot(ff)
   modelbound.boot(ff)
   
   gg = qgcomp.boot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=7, B=5, family=binomial(), rr=TRUE)
   plot(gg)
   pointwisebound.boot(gg)
   modelbound.boot(gg)
   
   
# gaussian
   dat <- data.frame(y=rnorm(n), x1=runif(n), x2=runif(n), z=runif(n))
   ee = qgcomp.noboot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=7, family=gaussian())
   plot(ee)
   pointwisebound.noboot(ee) 
   ff = qgcomp.boot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=7, B=8, family=gaussian())
   plot(ff)
   modelbound.boot(ff)
   
# poisson
   dat <- data.frame(y=rpois(n, 1.2), x1=runif(n), x2=runif(n), z=runif(n))
   ee = qgcomp.noboot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=7, family=poisson())
   plot(ee)
   pointwisebound.noboot(ee) 
   ff = qgcomp.boot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=7, B=5, family=poisson())
   modelbound.boot(ff)
   plot(ff)
      
#cox
   dat <- data.frame(stop=(tmg <- pmin(.1,rweibull(n, 10, 0.1))), 
                     start = pmax(0, tmg-runif(n)),
                     d=1.0*(tmg<0.1), x1=runif(n), x2=runif(n), z=runif(n))
   expnms=paste0("x", 1:2)
   f = survival::Surv(start,stop, d)~x1 + x2
   suppressWarnings(ee <- qgcomp.cox.noboot(f, expnms = expnms, data = dat))
   #pointwisebound.noboot(ee) # not working
   
   plot(ee)
   suppressWarnings(ff <- qgcomp.cox.boot(f, expnms = expnms, data = dat, B=2, MCsize=1000))
   plot(ff)
   #modelbound.boot(ff, pwonly=TRUE) # not working
   
   
# zi
  dat <- data.frame(y=rbinom(n, 1, 0.5)*rpois(n, 1.2), x1=runif(n), x2=runif(n), z=runif(n))
  ee = qgcomp.zi.noboot(f=y ~ z + x1 + x2 | z, expnms = c('x1', 'x2'), data=dat, q=7, dist="negbin")
  plot(ee)
  #pointwisebound.noboot(ee) # not working
  ffz = qgcomp.zi.boot(f=y ~ z + x1 + x2 | z, expnms = c('x1', 'x2'), data=dat, q=7, B=2, MCsize=1000, dist="negbin")
  pointwisebound.boot(ffz)
  modelbound.boot(ffz, pwonly=TRUE)
  plot(ffz)

  
cat("done")
