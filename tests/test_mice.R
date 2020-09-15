cat("# multiple imputation through chained equations test\n")
library("qgcomp")

N = 100
set.seed(123)
dat <- data.frame(y=runif(N), x1=runif(N), x2=runif(N), z=runif(N))
true = qgcomp.noboot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), 
                     data=dat, q=2, family=gaussian())
mdat <- dat
mdat$x1 = ifelse(mdat$x1>0.5, mdat$x1, NA)
mdat$x2 = ifelse(mdat$x2>0.75, mdat$x2, NA)
true <- qgcomp.noboot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), 
                    data=dat, q=2, family=gaussian())
cc <- qgcomp.noboot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), 
                    data=mdat[complete.cases(mdat),], q=2, family=gaussian())
                    

cdat = mdat
cdat$x1 = ifelse(is.na(cdat$x1), 0.5/sqrt(2), cdat$x1)
#data = cdat
ff <- function(data){
  nms = names(data)
  j = which(nms == "x2")
  f <- function(){
    nms = names(data)
    res = mice.impute.leftcenslognorm(y=cdat$x2, 
                              ry=!is.na(cdat$x2), 
                              x=cdat[,c("x1", "y", "z")], 
                              wy=is.na(cdat$x2), 
                              lod=NULL, 
                              debug=TRUE)
    res
  }
  f()
}

# works when LOD is not specified and based on minimum non-missing value
ff(cdat)


f0 <- function(data, j){
  print(sys.parent())
  mice.impute.leftcenslognorm(y=cdat$x2, 
                                    ry=!is.na(cdat$x2), 
                                    x=cdat[,c("x1", "y", "z")], 
                                    wy=is.na(cdat$x2), 
                                    lod=c(NA, 0.5, 0.75, NA), 
                                    debug=TRUE)
}

f1 <- function(data, j){
  print(sys.parent())
  #print(eval(as.name("data"), envir = parent.frame(n=1)))
  f0(data, j=j)
}

f2 <- function(data, j){
  print(sys.parent())
  f1(data, j)
}
f3 <- function(data, j){
  print(sys.parent())
  f2(data, j)
}
f4 <- function(data, j){
  print(sys.parent())
  f3(data, j)
}
# none of these appears to work because "j" is never found
f1(cdat, 3)
f2(cdat, 3)
f3(cdat, 3)
f4(cdat, 3)


# # note the following example imputes from the wrong parametric model and is expected to be biased
# # as a result
# library("mice")
# library("survival")
# set.seed(1231)
# impdat = mice(data = mdat, 
#               method = c("", "leftcenslognorm", "leftcenslognorm", ""),
#               lod=c(NA, 0.5, 0.75, NA), debug=FALSE, maxit = 10, m = 50)
# qc.fit.imp <- list(
#   call = call("qgcomp.noboot(y~., expnms = c('x1', 'x2'), family=gaussian())"),
#   call1 = impdat$call,
#   nmis = impdat$nmis,
#   analyses = lapply(1:50, function(x) qgcomp.noboot(y~., expnms = c("x1", "x2"),
#                                                    data=complete(impdat, x), family=gaussian(), bayes=FALSE))
# )
# qc.fit.imp2 <- list(
#   call = call("qgcomp.noboot(y~., expnms = c('x1', 'x2'), family=gaussian())"),
#   call1 = impdat$call,
#   nmis = impdat$nmis,
#   analyses = lapply(1:50, function(x) lm(y~., data=complete(impdat, x)))
# )
# 
# summary(complete(impdat, 1))
# obj <- pool(as.mira(qc.fit.imp))# are results at a given seed numerically stable across versions?
# 
# 
# summary(obj)
# #                 estimate  std.error statistic       df      p.value
# #  (Intercept)  0.67956728 0.09189112  7.395353 64.56892 3.582918e-10
# #  psi1        -0.09551483 0.04876708 -1.958592 49.71615 5.578245e-02
# 
# true
# #              Estimate Std. Error Lower CI  Upper CI t value
# # (Intercept)  0.582729   0.068033  0.44939 0.716071 1.663e-13
# # psi1        -0.099867   0.086641 -0.26968 0.069947    0.2519
# 
# cc
# #             Estimate Std. Error Lower CI  Upper CI t value
# # (Intercept)  0.58390    0.11158  0.36520  0.802600  0.0008
# # psi1        -0.29313    0.12952 -0.54697 -0.039282  0.0534
# 
# 
# obj2 <- pool(as.mira(qc.fit.imp2))# are results at a given seed numerically stable across versions?
# true = lm(y~., data=dat)
# cc = lm(y~., data=cdat)
# 
# summary(obj2)
# #                estimate  std.error  statistic       df      p.value
# # (Intercept)  0.88780523 0.21350133  4.1583123 52.88431 0.0001180606
# # x1          -0.19722820 0.14363657 -1.3731057 77.01935 0.1737039541
# # x2          -0.35849252 0.24710409 -1.4507754 47.94448 0.1533539928
# # z           -0.08929342 0.09794892 -0.9116325 91.46288 0.3643587309
# 
# 
# summary(true)
# #             Estimate Std. Error t value Pr(>|t|)    
# # (Intercept)  0.62516    0.09538   6.554 2.79e-09 ***
# # x1          -0.10569    0.10982  -0.962    0.338    
# # x2          -0.07144    0.09854  -0.725    0.470    
# # z           -0.07624    0.09822  -0.776    0.440    
# 
# summary(cc)
# #             Estimate Std. Error t value Pr(>|t|)
# # (Intercept)  -0.4987     0.9012  -0.553    0.586
# # x1           -0.1079     0.3138  -0.344    0.734
# # x2            1.1194     0.8712   1.285    0.214
# # z            -0.0159     0.1771  -0.090    0.929
# 
