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
   ff = qgcomp.boot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=7, B=5, family=binomial(), rr=FALSE)
   plot(ff)
   ff$msmfit$family$link
   pointwisebound.boot(ff)
   qgcomp:::pointwisebound.boot_old(ff)
   modelbound.boot(ff)
   qgcomp:::modelbound.boot_old(ff)
   
   gg = qgcomp.boot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=7, B=5, family=binomial(), rr=TRUE)
   plot(gg)
   pointwisebound.boot(gg)
   modelbound.boot(gg)
   
   
# gaussian
   dat <- data.frame(y=rnorm(n), x1=runif(n), x2=runif(n), z=runif(n))
   ee = qgcomp.noboot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=7, family=gaussian())
   plot(ee)
   pointwisebound.noboot(ee) 
   qgcomp:::pointwisebound.noboot_old(ee) 
   ff = qgcomp.boot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=3, B=3, family=gaussian())
   plot(ff, flexfit = FALSE)
   modelbound.boot(ff)
   qgcomp:::modelbound.boot_old(ff)
   
# poisson
   dat <- data.frame(y=rpois(n, 1.2), x1=runif(n), x2=runif(n), z=runif(n))
   ee = qgcomp.noboot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=7, family=poisson())
   plot(ee)
   pointwisebound.noboot(ee) 
   qgcomp:::pointwisebound.noboot_old(ee) 
   ff = qgcomp.boot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=7, B=5, family=poisson())
   modelbound.boot(ff)
   qgcomp:::modelbound.boot_old(ff)
   plot(ff)
      
#cox
   dat <- data.frame(stop=(tmg <- pmin(.1,rweibull(n, 10, 0.1))), 
                     start = pmax(0, tmg-runif(n)),
                     d=1.0*(tmg<0.1), x1=runif(n), x2=runif(n), z=runif(n))
   expnms=paste0("x", 1:2)
   f = survival::Surv(start,stop, d)~x1 + x2
   suppressWarnings(ee <- qgcomp.cox.noboot(f, expnms = expnms, data = dat))
   res = try(pointwisebound.noboot(ee), silent = TRUE) # not working
   stopifnot(class(res)=="try-error")   
   plot(ee)
   suppressWarnings(ff <- qgcomp.cox.boot(f, expnms = expnms, data = dat, B=2, MCsize=500))
   plot(ff)
   res = try(modelbound.boot(ff, pwonly=TRUE), silent = TRUE) # not working
   stopifnot(class(res)=="try-error")   
   curvelist = qgcomp.survcurve.boot(ff)
   stopifnot(inherits(curvelist, "list"))
   
# zi
  dat <- data.frame(y=rbinom(n, 1, 0.5)*rpois(n, 1.2), x1=runif(n), x2=runif(n), z=runif(n))
  ee = qgcomp.zi.noboot(f=y ~ z + x1 + x2 | z, expnms = c('x1', 'x2'), data=dat, q=7, dist="negbin")
  plot(ee)
  res = try(pointwisebound.noboot(ee), silent = TRUE) # not working
  stopifnot(class(res)=="try-error")   
  ffz = qgcomp.zi.boot(f=y ~ z + x1 + x2 | z, expnms = c('x1', 'x2'), data=dat, q=7, B=2, MCsize=500, dist="negbin")
  pointwisebound.boot(ffz)
  modelbound.boot(ffz, pwonly=TRUE)
  plot(ffz)

  
cat("done")
