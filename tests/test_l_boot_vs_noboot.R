cat("# boot vs no boot test\n")
library("qgcomp")

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
     y = rnorm(N, x1+x2, 2)
})
 dat[,c('y', paste0("x", 1:6))]
}

Xnm = c(paste0("x", 1:6))

repit <- function(i){
  dat = dgm(50)
  m1 = qgcomp.noboot(y~., expnms=Xnm, data = dat, family=gaussian(), q=4)
  m2 = qgcomp.boot(  y~., expnms=Xnm, data = dat, family=gaussian(), q=4, B=2, parallel=FALSE)
  res = c(m1$coef, m1$var.coef, 1*(m1$pval>0.05), with(m1, ci.coef[1]<2 & ci.coef[2]>2), m2$coef, m2$var.coef, 1*(m2$pval>0.05), with(m2, ci.coef[2,1]<2 & ci.coef[2,2]>2))
  names(res) <- c("psiint", "psi", "varint", "var",  "powint", "pow",  "cover", "b.psiint", "b.psi", "b.varint", "b.var", "b.powint", "b.pow", "b.cover")
  res
}


#res = mclapply(1:1000, repit)
res = lapply(1:2, repit)
res = simplify2array(res)

# equality within toleraance
stopifnot(all.equal(res["psiint",],res["b.psiint",]))
stopifnot(all.equal(res["psi",],res["b.psi",]))

#' \dontest{
#' # bootstrap and regular variance good
#' repit2 <- function(i){
#'   dat = dgm(500)
#'   m1 = qgcomp.noboot(y~., expnms=c("x1", "x2"), data = dat, family=gaussian(), q=4)
#'   m2 = qgcomp.boot(  y~., expnms=c("x1", "x2"), data = dat, family=gaussian(), q=4, B=5, parallel=TRUE)
#'   res = c(m1$coef, m1$var.coef, 1*(m1$pval>0.05), with(m1, ci.coef[1]<2 & ci.coef[2]>2), m2$coef, m2$var.coef, 1*(m2$pval>0.05), with(m2, ci.coef[2,1]<2 & ci.coef[2,2]>2))
#'   names(res) <- c("psiint", "psi", "varint", "var",  "powint", "pow",  "cover", "b.psiint", "b.psi", "b.varint", "b.var", "b.powint", "b.pow", "b.cover")
#'   res
#' }
#' 
#' res = lapply(1:2, repit2)
#' res = simplify2array(res)
#' 
#' 
#' stopifnot(all.equal(res["var",],res["b.var",], tolerance=sqrt(0.01)))
#' }


cat("done")
