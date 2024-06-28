cat("# ee test\n")
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

#checknames = qgcomp:::checknames
#.qgcomp_object = qgcomp:::.qgcomp_object

res = qgcomp.glm.ee(f=deaths~x1 + x2 + I(x2^2)+ I(x1^2), offset="logtime", expnms=c("x1", "x2"), data = dat, family=poisson(), q=7) # offset not working
pointwisebound.noboot(res)
plot(res, modelband = FALSE, pointwisebars = FALSE, flexfit = TRUE)
print(res)
summary(res)
print(res$fit)
print(res$msmfit)
predict(res)
predict(res$msmfit)
summary(glm(deaths ~ x2, data = dat, family=poisson(), offset = logtime))$coefficients

res = qgcomp.noboot(deaths~x2 + offset(logtime), expnms="x2", data = dat, family=poisson(), q=NULL)
rr = qgcomp.noboot(deaths~x2, expnms="x2", data = dat, family=poisson(), q=NULL)
res = qgcomp.glm.ee(deaths~x2, offset="logtime", expnms=c("x2"), data = dat, family=poisson(), q=NULL) # not working
rr = qgcomp.noboot(deaths~x2+x1, expnms=c('x1', "x2"), data = dat, family=poisson(), q=NULL)
#res2 = qgcomp.glm.ee(deaths~x2, offset="logtime", expnms=c("x2"), data = dat, family=poisson(), q=4) # not working
res2 = qgcomp.glm.ee(deaths~x2+x1, expnms=c('x1', "x2"), data = dat, family=poisson(), q=NULL)
summary(res)

res = qgcomp.noboot(deaths~x2+x1 + offset(logtime), expnms=c('x1', "x2"), data = dat, family=quasipoisson(), q=4)
res = qgcomp.glm.ee(deaths~x2+x1 + offset(logtime), expnms=c('x1', "x2"), data = dat, family=poisson(), q=4)
#res = qgcomp.glm.ee(deaths~x2+x1 , expnms=c('x1', "x2"), data = dat, family=quasipoisson(), q=4)


#TODO: offset