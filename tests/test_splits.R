cat("sample splits testing")
library(qgcomp)
set.seed(123223)
dat = qgcomp:::.dgm_quantized(N=1000, coef=c(0.25,-0.25,0,0), ncor=1)
cor(dat)
# overall fit (more or less null due to counteracting exposures)
(overall <- qgcomp.noboot(f=y~., q=NULL, expnms=c("x1", "x2", "x3", "x4"), data=dat))

# partial effects using 40% training/60% validation split
trainidx <- sample(1:nrow(dat), round(nrow(dat)*0.4))
valididx <- setdiff(1:nrow(dat),trainidx)
traindata = dat[trainidx,]
validdata = dat[valididx,]
splitres <- qgcomp:::qgcomp.partials(fun="qgcomp.noboot", f=y~., q=NULL, 
    traindata=traindata,validdata=validdata, expnms=c("x1", "x2", "x3", "x4"))
splitres
# under the null, both should give null results
set.seed(123223)
dat <- qgcomp:::.dgm_quantized(N=1000, coef=c(0,0,0,0), ncor=1)
# 40% training/60% validation
trainidx2 <- sample(1:nrow(dat), round(nrow(dat)*0.4))
valididx2 <- setdiff(1:nrow(dat),trainidx2)
traindata2 <- dat[trainidx2,]
validdata2 <- dat[valididx2,]
splitres2 <- qgcomp:::qgcomp.partials(fun="qgcomp.noboot", f=y~., 
   q=NULL, traindata=traindata2,validdata=validdata2, expnms=c("x1", "x2", "x3", "x4"))
splitres2
