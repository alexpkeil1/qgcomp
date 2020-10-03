cat("# basics test\n")
# se_comb
COV = matrix(c(.1, .2, .2, .1), nrow=2)
colnames(COV) <- c("x1", "x2")
stopifnot(sqrt(sum(COV))==qgcomp:::se_comb(covmat=COV, expnms = c("x1", "x2")))


#vc_comb
colnames(COV)[1] <- c("(Intercept)")
stopifnot(COV==qgcomp:::vc_comb(aname="(Intercept)", c("x2"), covmat=COV, grad=1.0))


# grad.poly # anything better here?
for(deg in 1:3){
  stopifnot(all(!is.na(qgcomp:::grad.poly(intvals=c(1,2,3), degree=deg))))
}


# stats
set.seed(50)
# linear model
dat = qgcomp:::.dgm_quantized()
ft = qgcomp::qgcomp.noboot(f=y ~ x1 + x2 + x3 + x4, expnms = c('x1', 'x2'), data=dat, q=2, family=gaussian(), bayes=TRUE)

df.residual(ft)
vcov(ft)
AIC(ft)
BIC(ft)
logLik(ft)
anova(ft)
confint(ft) # not working currently

predict(ft)
predict(ft, newdata = dat[1,])
