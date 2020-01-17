cat("# bayesqgcomp test\n")
#devtools::install_github("alexpkeil1/qgcomp@dev")
library("qgcomp")
data("metals", package="qgcomp")

Xnm <- c(
  'arsenic','barium','cadmium','calcium','chromium','copper',
  'iron','lead','magnesium','manganese','mercury','selenium','silver',
  'sodium','zinc'
)

# continuous outcome, example with perfect collinearity
metals$cadmium2 = metals$cadmium
Xnm2 = c(Xnm, "cadmium2")
res = try(qc.fit <- qgcomp.noboot(y~.,dat=metals[,c(Xnm2, 'y')], family=gaussian()))
stopifnot(class(res)=="try-error")
# error

res = try(qc.fit <- qgcomp.noboot(y~.,dat=metals[,c(Xnm2, 'y')], family=gaussian(), bayes=TRUE))
stopifnot(class(res)=="qgcompfit")

# compare to results with colinear exposure removed
res = try(qc.fit2 <- qgcomp.noboot(y~.,dat=metals[,c(Xnm, 'y')], family=gaussian()))
stopifnot(class(res)=="qgcompfit")



## slight shrinkage with heavily non-linear models that don't converge
## non-Bayesian
#results5 = qgcomp(y~. + .^2 + .^3 + arsenic*cadmium,
#                  expnms=Xnm,
#                  metals[,c(Xnm, 'y')], family=gaussian(), q=10, B=10, 
#                  seed=125, degree=3)
#
#print(results5)
#plot(results5)
#results5$bootsamps # samples are wack
#
## Bayesian (takes a few mins)
#results5bayes = qgcomp.boot(y~. + .^2 + .^3 + arsenic*cadmium,
#                  expnms=Xnm,
#                  metals[,c(Xnm, 'y')], family=gaussian(), q=10, B=10, 
#                  seed=125, degree=3, 
#                  bayes=TRUE)
#
## note p-values are not computed due to negative
## degrees of freedom (can fix by using a z-statistic instead of a t-statistic)
#print(results5bayes)
#plot(results5bayes)
#results5bayes$bootsamps # samples are reasonable


cat("done")
