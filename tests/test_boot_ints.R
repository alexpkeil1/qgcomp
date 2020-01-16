cat("# boot ints test\n")
# edge case with user defined breaks in qgcomp boot where msm.fit doesn't reach all
# the quantiles
library(qgcomp)


z <- cbind(rnorm(100),rnorm(100),rnorm(100),rnorm(100))

y <- -1 + z %*% c(2,1,0.5,0) + z^2 %*% c(.2,.2,0.5,0) + rnorm(100)
Zn = paste0("z.", 1:4)
qft <- as.formula(
  paste("y ~",paste(Zn, collapse = "+"), "+", paste0(Zn, '*', Zn, collapse = "+"))
)

dat = data.frame(y=y,z=z)

#first define breaks
qdat <- quantize(dat, expnms = Zn, q = 6, breaks = NULL)

# now use these breaks
qdat2 <- quantize(dat, expnms = Zn, q=NULL, breaks=qdat$breaks)

# should be true
stopifnot(identical(qdat,qdat2))

res1 <- qgcomp.boot(qft, degree=2, expnms = Zn, data = dat, 
                    q = 6, 
                    B=20, seed = 12312)
# should be true
stopifnot(identical(res1$breaks, qdat$breaks))

res2 <- qgcomp.boot(qft, degree=2, expnms = Zn, data = dat, 
                    breaks = qdat$breaks, 
                    B=20, seed = 12312)
# should be true
stopifnot(identical(res2$breaks, qdat$breaks))

# should have identical x axes

# should be true
stopifnot(identical(res1$msmfit$df.null, res2$msmfit$df.null))

cat("done")
