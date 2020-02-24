cat("# boot chooser test\n")
library("qgcomp")
library("splines")
library("survival")

dgm <- function(N){
 dat <- data.frame(id=1:N) 
 dat <- within(dat, {
     time=(tmg <- pmin(.1,rweibull(N, 10, 0.1)))
     d=1.0*(tmg<0.1)
     u = 0
     x1 = runif(N)*4 + u
     x2 = runif(N)*4 + u
     x3 = runif(N)*4 + u
     x4 = runif(N)*4 + u
     x5 = runif(N)*4 + u
     x6 = runif(N)*4 + u
     y = rnorm(N, x1+x2, 2)
})
 dat[,c("y", "time", "d", paste0("x", 1:6))]
}

set.seed(11232)
dat = dgm(200)
Xnm = c(paste0("x", 1:6))


f = y~x1 + x2 + x3 + I(x3^2)
qgcomp(f, expnms = c("x1", "x2"), data = dat)
qgcomp.noboot(f, expnms = c("x1", "x2"), data = dat)



# should cause errors
f = y~. + .^2
res = try(qgcomp(f, data = dat), silent=TRUE)
stopifnot(class(res)=="try-error")

fs = Surv(time,d)~. + .^2
res = try(qgcomp(fs, data = dat), silent=TRUE)
stopifnot(class(res)=="try-error")


f = y~. 
res = try(fit1 <- qgcomp(f, expnms=Xnm, data = dat, parallel=TRUE), silent=TRUE)
stopifnot(class(res)=="try-error")

fs = Surv(time,d)~.
res = try(fit1 <- qgcomp(fs, expnms=Xnm, data = dat, B=5, MCsize=100), silent=TRUE)
stopifnot(class(res)=="try-error")

fs = Surv(time,d)~.
res = try(fit1 <- qgcomp(fs, expnms=Xnm, family=cox(), data = dat), silent=TRUE)
stopifnot(class(res)=="try-error")


f = y~. + .^2
res = try(fit1 <- qgcomp(f, expnms=Xnm, data = dat), silent=TRUE)
stopifnot(class(res)=="qgcompfit")

fs = Surv(time,d)~. + .^2
res = try(fit1 <- qgcomp(fs, expnms=Xnm, data = dat, B=5, MCsize=100), silent=TRUE)
stopifnot(class(res)=="qgcompfit")

fs = Surv(time,d)~. + .^2
res = try(fit1 <- qgcomp(fs, expnms=Xnm, data = dat, B=5, MCsize=100, parallel=TRUE), silent=TRUE)
stopifnot(class(res)=="qgcompfit")



# splines splines do work
f = y ~ x2 + x3 + x4 + x5 + x6 + splines::ns(x1, df=2)
res = try(qgcomp(f, data = dat), silent=TRUE)
stopifnot(class(res)=="try-error") # should give error that expnms not defined
res = try(fit1 <- qgcomp(f, expnms=Xnm, q=8, data = dat, deg=2), silent=TRUE)
stopifnot(class(res)=="qgcompfit")



# splines splines do work
f = y ~ x2 + x3 + x4 + x5 + x6 + splines::ns(x1, df=2)
res = try(qgcomp(f, data = dat), silent=TRUE)
stopifnot(class(res)=="try-error") # should give error that expnms not defined
res = try(fit1 <- qgcomp(f, expnms=Xnm, q=8, data = dat, deg=2), silent=TRUE)
stopifnot(class(res)=="qgcompfit")

# splines splines + bayes
f = y ~ x2 + x3 + x4 + x5 + x6 + splines::ns(x1, df=2)
res = try(fit2 <- qgcomp(f, expnms=Xnm, q=8, data = dat, deg=2, bayes=TRUE), silent=TRUE)
stopifnot(class(res)=="qgcompfit")


# indicator functions
f = y ~ factor(x1) + x2 + x3 + x4 + x5 + x6
res = try(qgcomp(f, data = dat), silent=TRUE)
stopifnot(class(res)=="try-error") # should give error that expnms not defined
res = try(fit1 <- qgcomp(f, expnms=Xnm, q=8, data = dat, deg=3), silent=TRUE)
stopifnot(class(res)=="qgcompfit")



cat("done")
