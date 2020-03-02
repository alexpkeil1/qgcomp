cat("# poisson test\n")
library("qgcomp")

# grouped survival data
N = 250
t = 4
dat <- data.frame(row.names = 1:(N*t))
dat <- within(dat, {
  group = do.call("c", lapply(1:N, function(x) rep(x, t)))
  age = do.call("c", lapply(1:N, function(x) 1:t))
  u =  do.call("c", lapply(1:N, function(x) rep(runif(1), t)))
  x1 = rnorm(N, u)
  x2 = do.call("c", lapply(1:N, function(x) rep(sample(0:3,1), t)))
  time = runif(N*t, 20, 200)
  logtime = log(time)
  lograte = -5 + x1
  deaths = rpois(N, exp(log(time) + lograte)) 
})

summary(glm(deaths ~ x2, data = dat, family=poisson(), offset = logtime))$coefficients

qgcomp.noboot(deaths~x2 + offset(logtime), expnms="x2", data = dat, family=poisson(), q=NULL)
qgcomp.noboot(deaths~x2+x1 + offset(logtime), expnms=c('x1', "x2"), data = dat, family=poisson(), q=4)


# error by scoping for using the offset parameter
#qgcomp.noboot(deaths~x2, expnms="x2", data = dat, family=poisson(), q=4, offset = logtime)

# this doesn't fail, doesn't give correct intercept (offset is forgotten) and 
# variance looks wrong
#qgcomp.boot(deaths~x2 + offset(logtime), expnms="x2", data = dat, family=poisson(), parallel=TRUE,
#            q=NULL, B=500, MCsize = 1000)



# matching bootstrapped with non-bootstrapped version
dgm <- function(N){
 dat <- data.frame(id=1:N) 
 dat <- within(dat, {
     u = 0
     x1 = runif(N)*4 + u
     x2 = runif(N)*4 + u
     x3 = runif(N)*4 + u
     x4 = runif(N)*4 + u
     x5 = runif(N)*4 + u
     x6 = runif(N)*4 + u
     y = rpois(N, exp(-5+.3*x1+.3*x2))
})
 dat[,c('y', paste0("x", 1:6))]
}

Xnm = c(paste0("x", 1:6))

dat = dgm(200)
m1 = qgcomp.noboot(y~., expnms=Xnm, data = dat, family=gaussian(), q=4)
m1a = qgcomp(y~., expnms=Xnm, data = dat, family=gaussian(), q=4)
m2 = qgcomp.boot(  y~., expnms=Xnm, data = dat, family=gaussian(), q=4, B=5, parallel=TRUE, MCsize = 50)
print(coef(m1), digits=10)
print(coef(m2$msmfit), digits=10)

#' \donttest{
#' m1 = qgcomp.noboot(y~., expnms=Xnm, data = dat, family=poisson(), q=4)
#' m2 = qgcomp.boot(  y~., expnms=Xnm, data = dat, family=poisson(), q=4, B=5, parallel=TRUE, MCsize = 50)
#' print(coef(m1), digits=10)
#' print(coef(m2$msmfit), digits=10)
#' 
#' 
#' repit <- function(i){
#'   dat = dgm(1000)
#'   m1 = qgcomp.noboot(y~., expnms=Xnm, data = dat, family=poisson(), q=4)
#'   m2 = qgcomp.boot(  y~., expnms=Xnm, data = dat, family=poisson(), q=4, B=5, parallel=TRUE, MCsize = 100000)
#'   res = c(m1$coef, m1$var.coef, 1*(m1$pval>0.05), with(m1, ci.coef[1]<2 & ci.coef[2]>2), m2$coef, m2$var.coef, 1*(m2$pval>0.05), with(m2, ci.coef[2,1]<2 & ci.coef[2,2]>2))
#'   names(res) <- c("psiint", "psi", "varint", "var",  "powint", "pow",  "cover", "b.psiint", "b.psi", "b.varint", "b.var", "b.powint", "b.pow", "b.cover")
#'   res
#' }
#' 
#' 
#' #res = mclapply(1:1000, repit)
#' res = lapply(1:2, repit)
#' res = simplify2array(res)
#' 
#' # equality within toleraance
#' stopifnot(all.equal(res["psiint",],res["b.psiint",], tolerance=sqrt(0.01)))
#' stopifnot(all.equal(res["psi",],res["b.psi",], tolerance=sqrt(0.01)))
#' 
#' 
#' # bootstrap and regular variance good
#' repit2 <- function(i){
#'   dat = dgm(500)
#'   m1 = qgcomp.noboot(y~., expnms=c("x1", "x2"), data = dat, family=poisson(), q=4)
#'   m2 = qgcomp.boot(  y~., expnms=c("x1", "x2"), data = dat, family=poisson(), q=4, B=5, parallel=TRUE, MCsize = 1000)
#'   res = c(m1$coef, m1$var.coef, 1*(m1$pval>0.05), with(m1, ci.coef[1]<2 & ci.coef[2]>2), m2$coef, m2$var.coef, 1*(m2$pval>0.05), with(m2, ci.coef[2,1]<2 & ci.coef[2,2]>2))
#'   names(res) <- c("psiint", "psi", "varint", "var",  "powint", "pow",  "cover", "b.psiint", "b.psi", "b.varint", "b.var", "b.powint", "b.pow", "b.cover")
#'   res
#' }
#' 
#' #res = lapply(1:500, repit2)
#' #res = simplify2array(res)
#' 
#' 
#' #stopifnot(all.equal(res["var",],res["b.var",], tolerance=sqrt(0.01)))
#' 
#' }
  cat("done")
 
