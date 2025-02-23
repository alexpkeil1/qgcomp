.esteq_qgc <- function(theta, family, Y, X, Xint, Xmsm, weights,rr=FALSE, offset=0, delta=-Inf, ...){
  # family specific estimating equation based on input matrixes
  # used for "A" part of sandwich variance: V = solve(A) %*% B %&% t(solve(A))
  fam = family$family
  lnk = family$link
  binlink = ifelse(rr, "logitlog", "logit")
  FUN <- switch(fam, 
                gaussian = .esteq_qgclin,
                tobit = .esteq_qgctobit,
                poisson = .esteq_qgcpoisson,
                #binomial = .esteq_qgclogit(theta, Y, X, Xint, Xmsm, weights),
                binomial = switch(binlink,
                                  logit = .esteq_qgclogit,
                                  logitlog = .esteq_qgclogitlog
                ),
                .esteq_error
  )
  FUN(theta=theta, Y=Y, X=X, Xint=Xint, Xmsm=Xmsm, weights=weights, offset, delta=-Inf, ...)
}

.esteq_qgcdf <- function(f, data, theta, family, intvals, expnms, hasintercept, weights, degree=1,rr=FALSE,offset=0, delta=-Inf, ...){
  # family specific estimating equation based on input dataframes
  # used for "B" part of sandwich variance: V = solve(A) %*% B %&% t(solve(A)) because it facilitates
  # individual level calculations. 
  fam = family$family
  lnk = family$link
  X = model.matrix(f, data) # model.frame
  Y = model.response(data)
  Xint = as.matrix(do.call(rbind,lapply(intvals, function(x) {data[,expnms] = x; model.matrix(f,data)})))
  Xmsm = poly(Xint[, expnms[1]], degree=degree, raw=TRUE) # intercept and constant exposure
  if(hasintercept){
    Xmsm = cbind(Xint[,colnames(Xint)[1]], Xmsm)
  }
  binlink = ifelse(rr, "logitlog", "logit")
  FUN <- switch(fam, 
                gaussian = .esteq_qgclin,
                tobit = .esteq_qgctobit,
                poisson = .esteq_qgcpoisson,
                #binomial = .esteq_qgclogit(theta, Y, X, Xint, Xmsm, weights),
                binomial = switch(binlink,
                                  logit = .esteq_qgclogit,
                                  logitlog = .esteq_qgclogitlog
                ),
                .esteq_error
  )
  FUN(theta=theta, Y=Y, X=X, Xint=Xint, Xmsm=Xmsm, weights=weights, offset=offset, delta=-Inf, ...)
}


.esteq_error <- function(...){
  stop(paste0("Family (distribution) not (yet) supported"))
}

.esteq_qgctobit <- function(theta, Y, X, Xint, Xmsm, weights, delta=-Inf,offset=0){
  stop("Tobit not yet functioning")
  # linear model gradient-based estimating equations
  npmsm = ncol(Xmsm)
  np = ncol(X)
  intweights = rep(weights, times=nrow(Xmsm)/nrow(X))
  stopifnot(length(theta)==np+npmsm+1)
  # estimating functions
  fun1 = do.call(c,lapply(1:np, function(x){
    mu = X %*% theta[1:np]
    eps = delta - mu
    sigma = theta[np+npmsm+1]
    t(weights * (Y - mu)* (Y>delta)) %*% X[,x]  + 
      t(weights * (Y<=delta) * dnorm(eps, 0, sigma)/pnorm(eps, 0, sigma)) %*% X[,x]
  } ))
  fun2 = do.call(c,lapply(1:npmsm, function(x) {
    mu = Xmsm %*% theta[(np+1):(np+npmsm)]
    t(intweights * (Xint %*% theta[1:np] - mu)) %*% Xmsm[,x]
  }))
  fun3 = do.call(c,lapply(1:np, function(x){
    mu = X %*% theta[1:np]
    eps = delta - mu
    sigma = theta[np+npmsm+1]
    t(weights * (Y - mu)* (Y>delta)) %*% X[,x]
  } ))
  c(
    fun1,fun2
  )
}

.esteq_qgclin <- function(theta, Y, X, Xint, Xmsm, weights, delta=-Inf,offset=0){
  # linear model gradient-based estimating equations
  npmsm = ncol(Xmsm)
  np = ncol(X)
  intweights = rep(weights, times=nrow(Xmsm)/nrow(X))
  stopifnot(length(theta)==np+npmsm)
  # estimating functions
  fun1 = do.call(c,lapply(1:np, function(x){
    mu = X %*% theta[1:np]
    t(weights * (Y - mu)) %*% X[,x]
  } ))
  fun2 = do.call(c,lapply(1:npmsm, function(x) {
    mu = Xmsm %*% theta[(np+1):(np+npmsm)]
    t(intweights * (Xint %*% theta[1:np] - mu)) %*% Xmsm[,x]
    }))
  c(
    fun1,fun2
  )
}

.esteq_qgcpoisson <- function(theta, Y, X, Xint, Xmsm, weights, delta=-Inf, offset=0){
  # Poisson model gradient-based estimating equations
  npmsm = ncol(Xmsm)
  np = ncol(X)
  intweights = rep(weights, times=nrow(Xmsm)/nrow(X))
  stopifnot(length(theta)==np+npmsm)
  # estimating functions
  fun1 = do.call(c,lapply(1:np, function(x){
    mu = X %*% theta[1:np]
    t(weights * (Y - exp(mu))) %*% X[,x]
  }))
  fun2 = do.call(c,lapply(1:npmsm, function(x){
    ypred = exp(Xint %*% theta[1:np])
    mu = Xmsm %*% theta[(np+1):(np+npmsm)] 
    t(intweights * (ypred  - exp(mu))) %*% Xmsm[,x]
  } ))
  c(
    fun1,fun2
  )
}

.esteq_qgclogit <- function(theta, Y, X, Xint, Xmsm, weights, delta=-Inf, offset=0){
  # logistic model gradient-based estimating equations
  if(length(offset==0)) offset = 0
  npmsm = ncol(Xmsm)
  np = ncol(X)
  intweights = rep(weights, times=nrow(Xmsm)/nrow(X))
  stopifnot(length(theta)==np+npmsm)
  # estimating functions
  fun1 = do.call(c,lapply(1:np, function(x){
    mu = X %*% theta[1:np] 
    t(weights * (Y - .expit(mu))) %*% X[,x]
  }))
  fun2 = do.call(c,lapply(1:npmsm, function(x){
    ypred = .expit(Xint %*% theta[1:np])
    mu = Xmsm %*% theta[(np+1):(np+npmsm)]
    t(intweights * (ypred  - .expit(mu))) %*% Xmsm[,x]
  } ))
  c(
    fun1,fun2
  )
}

.esteq_qgclogitlog <- function(theta, Y, X, Xint, Xmsm, weights, delta=-Inf, offset=0){
  # logistic model gradient-based estimating equations with log-linear MSM
  npmsm = ncol(Xmsm)
  np = ncol(X)
  intweights = rep(weights, times=nrow(Xmsm)/nrow(X))
  stopifnot(length(theta)==np+npmsm)
  # estimating functions
  fun1 = do.call(c,lapply(1:np, function(x){
    mu = X %*% theta[1:np]
    t(weights * (exp(-mu) * Y - 1 + Y)/(1+exp(-mu))) %*% X[,x]
    #sum((exp(-mu) * Y - 1 + Y)/(1+exp(-mu)) * X[,x])
  }))
  fun2 = do.call(c,lapply(1:npmsm, function(x){
    ypred = .expit(Xint %*% theta[1:np])
    mu = Xmsm %*% theta[(np+1):(np+npmsm)]
    t(intweights * (ypred - exp(mu))/(1-exp(mu))) %*% Xmsm[,x]
  } ))
  c(
    fun1,fun2
  )
}

.esteq_qgclog <- function(theta, Y, X, Xint, Xmsm, weights, delta=-Inf, offset=0){
  # logistic model gradient-based estimating equations with log-linear MSM
  npmsm = ncol(Xmsm)
  np = ncol(X)
  intweights = rep(weights, times=nrow(Xmsm)/nrow(X))
  stopifnot(length(theta)==np+npmsm)
  # estimating functions
  fun1 = do.call(c,lapply(1:np, function(x){
    mu = X %*% theta[1:np]
    t(weights * (Y - exp(mu))/(1-exp(mu))) %*% X[,x]
    #sum((exp(-mu) * Y - 1 + Y)/(1+exp(-mu)) * X[,x])
  }))
  fun2 = do.call(c,lapply(1:npmsm, function(x){
    ypred = exp(Xint %*% theta[1:np])
    mu = Xmsm %*% theta[(np+1):(np+npmsm)]
    t(intweights * (ypred - exp(mu))/(1-exp(mu))) %*% Xmsm[,x]
  } ))
  c(
    fun1,fun2
  )
}

# 
# .esteq_qgclingeex = function(data, .intvals, .expnms, .hasintercept, .weights, .degree=1){
#   # following conventions of geex package, not used in package but was used for testing
#     X = model.matrix(f, data) # model.frame
#     Y = model.response(data)
#     Xint = as.matrix(do.call(rbind,lapply(.intvals, function(x) {data[,.expnms] = x; model.matrix(f,data)})))
#     Xmsm = Xint[,c(ifelse(.hasintercept, colnames(Xint)[1], NULL), .expnms[1])] # intercept and constant exposure
#     Xmsm = poly(Xint[, expnms[1]], degree=.degree, raw=TRUE) # intercept and constant exposure
#     if(.hasintercept){
#       Xmsm = cbind(Xint[,colnames(Xint)[1]], Xmsm)
#     }
#     function(theta, .Y=Y, .X=X, .Xint=Xint, .Xmsm=Xmsm){
#     .esteq_qgclin(theta, .Y, .X, .Xint, .Xmsm, .weights)
#   }
# }


#' @importFrom stats gaussian
tobit <- function(){
  dist = stats::gaussian()
  dist$family = "tobit"
  dist
  }



qgcomp.glm.ee <- function(
    f,
    data,
    expnms=NULL,
    q=4,
    breaks=NULL,
    id=NULL,
    weights,
    offset=NULL,
    alpha=0.05,
    rr=TRUE,
    degree=1,
    seed=NULL,
    includeX=TRUE,
    verbose=TRUE,
    ...
){
  
  #' @title Quantile g-computation for continuous and binary outcomes
  #' @aliases qgcomp.ee
  #'
  #' @description This function estimates a dose-response parameter representing a one quantile
  #' increase in a set of exposures of interest. This model estimates the parameters of a marginal
  #' structural model (MSM) based on g-computation with quantized exposures. Note: this function
  #' allows non-linear and non-additive effects of individual components of the exposure, as well as
  #' non-linear joint effects of the mixture via polynomial basis functions, which increase the
  #' computational computational burden due to the need for non-parametric bootstrapping.
  #' 
  #' Estimating equation methodology is used as the underlying estimation scheme. This allows that observations can be correlated, and is fundementally identical to some implementations of "generalized estimating equations" or GEE. Thus, it allows for a more general set of longitudinal data structures than does the qgcomp.glm.noboot function, because it allows that the outcome can vary over time within an individual. Interpretation of parameters is similar to that of a GEE: this function yields population average estimates of the effect of exposure over time. 
  #' 
  #' Note: GEE (and this function, by extension) does not automatically address all problems of longitudinal data, such as lag structures. Those are up to the investigator to specify correctly.
  #' 
  #' Note: qgcomp.ee is an equivalent function (slated for deprecation)
  #'
  #' @details Estimates correspond to the average expected change in the
  #'  (log) outcome per quantile increase in the joint exposure to all exposures
  #'  in `expnms'. Test statistics and confidence intervals are based on
  #'  a non-parametric bootstrap, using the standard deviation of the bootstrap
  #'  estimates to estimate the standard error. The bootstrap standard error is
  #'  then used to estimate Wald-type confidence intervals. Note that no bootstrapping
  #'  is done on estimated quantiles of exposure, so these are treated as fixed
  #'  quantities
  #'  
  #'  Note that there is an argument to "family" that is optional, but it defaults to gaussian
  #'  per the glm function. Current options allow "gaussian" (or gaussian()) or "binomial." This
  #'  function defaults to canonical links (binomial -> logistic link, gaussian ->
  #'  identity link). However, setting family="binomial" and rr=TRUE uses a 
  #'  marginal structural model with a log-link (and an underlying conditional model with
  #'  a logistic link). This mimics the behavior of \code{\link[qgcomp]{qgcomp.glm.boot}}
  #'  which maximizes backwards compatibility but also is a useful default to avoid
  #'  convergence issues when using the log-link for conditional models.
  #'
  #' @param f R style formula
  #' @param data data frame
  #' @param expnms character vector of exposures of interest
  #' @param q NULL or number of quantiles used to create quantile indicator variables
  #' representing the exposure variables. If NULL, then gcomp proceeds with un-transformed
  #' version of exposures in the input datasets (useful if data are already transformed,
  #' or for performing standard g-computation)
  #' @param breaks (optional) NULL, or a list of (equal length) numeric vectors that
  #' characterize the minimum value of each category for which to
  #' break up the variables named in expnms. This is an alternative to using 'q'
  #' to define cutpoints.
  #' @param id (optional) NULL, or variable name indexing individual units of
  #' observation (only needed if analyzing data with multiple observations per
  #' id/cluster). Note that qgcomp.glm.noboot will not produce cluster-appropriate
  #' standard errors. qgcomp.glm.ee can be used for this, which will use bootstrap
  #' sampling of clusters/individuals to estimate cluster-appropriate standard
  #' errors via bootstrapping.
  #' @param weights NULL or "case weights" - sampling weights representing a proportional representation in the analysis data
  #' @param offset Not yet implemented
  #' \code{\link[stats]{glm}} or \code{\link[arm]{bayesglm}}
  #' @param alpha alpha level for confidence limit calculation
  #' @param rr logical: if using binary outcome and rr=TRUE, qgcomp.glm.ee will
  #'   estimate risk ratio rather than odds ratio
  #' @param degree polynomial bases for marginal model (e.g. degree = 2
  #'  allows that the relationship between the whole exposure mixture and the outcome
  #'  is quadratic (default = 1).
  #' @param seed integer or NULL: random number seed for replicable bootstrap results
  #' @param includeX (logical) should the design/predictor matrix be included in the output via the fit and msmfit parameters?
  #' @param verbose (logical) give extra messages about fits
  #' @param ... arguments to glm (e.g. family)
  #' @family qgcomp_methods
  #' @return a qgcompfit object, which contains information about the effect
  #'  measure of interest (psi) and associated variance (var.psi), as well
  #'  as information on the model fit (fit) and information on the
  #'  marginal structural model (msmfit) used to estimate the final effect
  #'  estimates.
  #' @concept variance mixtures
  #' @import stats 
  #' @importFrom rootSolve multiroot
  #' @importFrom numDeriv jacobian
  #' @export
  #' @examples
  #' set.seed(30)
  #' # continuous outcome
  #' dat <- data.frame(y=rnorm(100), x1=runif(100), x2=runif(100), z=runif(100))
  #' # Conditional linear slope
  #' qgcomp.glm.noboot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=4, family=gaussian())
  #' # Marginal linear slope (population average slope, for a purely linear,
  #' #  additive model this will equal the conditional)
  #' qgcomp.glm.ee(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=4,
  #'   family=gaussian()) 
  #' qgcomp.glm.ee(f=y ~ x1 + x2 + I(x1*x2) + z, expnms = c('x1', 'x2'), data=dat, q=4,
  #'   family=gaussian()) 
  #' # no intercept model
  #' qgcomp.glm.ee(f=y ~ -1+z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=4,
  #'   family=gaussian()) 
  #'
  #'  # Note that these give different answers! In the first, the estimate is conditional on Z,
  #'  # but in the second, Z is marginalized over via standardization. The estimates
  #'  # can be made approximately the same by centering Z (for linear models), but
  #'  # the conditional estimate will typically have lower standard errors.
  #'  dat$z = dat$z - mean(dat$z)
  #'
  #' # Conditional linear slope
  #' qgcomp.glm.noboot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=4, family=gaussian())
  #' # Marginal linear slope (population average slope, for a purely linear,
  #' #  additive model this will equal the conditional)
  #'
  #' qgcomp.glm.ee(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=4,
  #'   family=gaussian()) 
  #'  \dontrun{
  #'
  #' # Population average mixture slope which accounts for non-linearity and interactions
  #' # Note this is one use case for the estimating equations approach: replacing bootstrapping
  #' qgcomp.glm.ee(y ~ z + x1 + x2 + I(x1^2) + I(x2*x1), family="gaussian",
  #'  expnms = c('x1', 'x2'), data=dat, q=4)
  #' qgcomp.glm.boot(y ~ z + x1 + x2 + I(x1^2) + I(x2*x1), family="gaussian",
  #'  expnms = c('x1', 'x2'), data=dat, q=4, B=1000)
  #'
  #' # generally non-linear/non-addiive underlying models lead to non-linear mixture slopes
  #' dat$y = dat$y + dat$x1*dat$x2
  #' qgcomp.glm.ee(y ~ z + x1 + x2 + I(x1^2) + I(x2*x1), family="gaussian",
  #'  expnms = c('x1', 'x2'), data=dat, q=4, deg=2)
  #' qgcomp.glm.boot(y ~ z + x1 + x2 + I(x1^2) + I(x2*x1), family="gaussian",
  #'  expnms = c('x1', 'x2'), data=dat, q=4, deg=2, B=1000)
  #'
  #' # binary outcome
  #' dat <- data.frame(y=rbinom(50,1,0.5), x1=runif(50), x2=runif(50), z=runif(50))
  #'
  #' # Conditional mixture OR
  #' qgcomp.glm.noboot(y ~ z + x1 + x2, family="binomial", expnms = c('x1', 'x2'),
  #'   data=dat, q=2)
  #'
  #' #Marginal mixture OR (population average OR - in general, this will not equal the
  #' # conditional mixture OR due to non-collapsibility of the OR)
  #' qgcomp.glm.boot(y ~ z + x1 + x2, family="binomial", expnms = c('x1', 'x2'),
  #'   data=dat, q=2, B=300, MCsize=5000, rr=FALSE)
  #' qgcomp.glm.ee(y ~ z + x1 + x2, family="binomial", expnms = c('x1', 'x2'),
  #'   data=dat, q=2, rr=FALSE)
  #'
  #' # Population average mixture RR
  #' qgcomp.glm.boot(y ~ z + x1 + x2, family="binomial", expnms = c('x1', 'x2'),
  #'   data=dat, q=2, B=300, MCsize=5000, rr=TRUE)
  #' qgcomp.glm.ee(y ~ z + x1 + x2, family="binomial", expnms = c('x1', 'x2'),
  #'   data=dat, q=2, rr=TRUE)
  #' # getting the RR with the poisson trick
  #' qgcomp.glm.ee(y ~ z + x1 + x2, family="poisson", expnms = c('x1', 'x2'),
  #'   data=dat, q=2)
  #' qgcomp.glm.boot(y ~ z + x1 + x2, family="poisson", expnms = c('x1', 'x2'),
  #'   data=dat, q=2, B=300, MCsize=5000)
  #'
  #' # Population average mixture RR, indicator variable representation of x2
  #' # note that I(x==...) operates on the quantile-based category of x,
  #' # rather than the raw value
  #' res = qgcomp.glm.ee(y ~ z + x1 + I(x2==1) + I(x2==2) + I(x2==3),
  #'   family="binomial", expnms = c('x1', 'x2'), data=dat, q=4, rr=TRUE)
  #' res$fit
  #' plot(res)
  #'
  #' # now add in a non-linear MSM
  #' res2a = qgcomp.glm.boot(y ~ z + x1 + I(x2==1) + I(x2==2) + I(x2==3),
  #'   family="binomial", expnms = c('x1', 'x2'), data=dat, q=4, rr=TRUE,
  #'   degree=2)
  #' res2 = qgcomp.glm.ee(y ~ z + x1 + I(x2==1) + I(x2==2) + I(x2==3),
  #'   family="binomial", expnms = c('x1', 'x2'), data=dat, q=4, rr=TRUE,
  #'   degree=2)
  #'   
  #' # conditional model estimates (identical point estimates and similar std. errors)
  #' summary(res2a$fit)$coefficients
  #' res2$fit
  #' # msm estimates (identical point estimates and different std. errors)
  #' summary(res2a$msmfit)$coefficients  # correct point estimates, incorrect standard errors
  #' res2  # correct point estimates, correct standard errors
  #' plot(res2)
  #' # Log risk ratio per one IQR change in all exposures (not on quantile basis)
  #' dat$x1iqr <- dat$x1/with(dat, diff(quantile(x1, c(.25, .75))))
  #' dat$x2iqr <- dat$x2/with(dat, diff(quantile(x2, c(.25, .75))))
  #' # note that I(x>...) now operates on the untransformed value of x,
  #' # rather than the quantized value
  #' res2 = qgcomp.glm.ee(y ~ z + x1iqr + I(x2iqr>0.1) + I(x2>0.4) + I(x2>0.9),
  #'   family="binomial", expnms = c('x1iqr', 'x2iqr'), data=dat, q=NULL, rr=TRUE,
  #'   degree=2)
  #' res2
  #'
  #'
  #' # weighted model
  #' N=5000
  #' dat4 <- data.frame(id=seq_len(N), x1=runif(N), x2=runif(N), z=runif(N))
  #' dat4$y <- with(dat4, rnorm(N, x1*z + z, 1))
  #' dat4$w=runif(N) + dat4$z*5
  #' qdata = quantize(dat4, expnms = c("x1", "x2"), q=4)$data
  #' # first equivalent models with no covariates
  #' qgcomp.glm.noboot(f=y ~ x1 + x2, expnms = c('x1', 'x2'), data=dat4, q=4, family=gaussian())
  #' qgcomp.glm.noboot(f=y ~ x1 + x2, expnms = c('x1', 'x2'), data=dat4, q=4, family=gaussian(),
  #'               weights=w)
  #'
  #' set.seed(13)
  #' qgcomp.glm.ee(f=y ~ x1 + x2, expnms = c('x1', 'x2'), data=dat4, q=4, family=gaussian())
  #' qgcomp.glm.ee(f=y ~ x1 + x2, expnms = c('x1', 'x2'), data=dat4, q=4, family=gaussian(),
  #'             weights=w)
  #' # using the correct model
  #' set.seed(13)
  #' qgcomp.glm.ee(f=y ~ x1*z + x2, expnms = c('x1', 'x2'), data=dat4, q=4, family=gaussian(),
  #'             weights=w, id="id")
  #' (qgcfit <- qgcomp.glm.ee(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat4, q=4,
  #'                        family=gaussian(), weights=w))
  #' qgcfit$fit
  #' summary(glm(y ~ z + x1 + x2, data = qdata, weights=w))
  #' }
  oldq = NULL
  if(is.null(seed)) seed = round(runif(1, min=0, max=1e8))
  
  newform <- terms(f, data = data)
  hasintercept = as.logical(attr(newform, "intercept"))
  class(newform) <- "formula"
  
  cc = match.call(expand.dots = TRUE)
  
  if(is.null(cc$delta)){
    delta=-Inf
  } else{
    delta = eval(cc$delta)
  }
  if(!is.null(cc$family) && cc$family != "tobit"){
    testfit <- glm(y~., data=data.frame(y=c(0,1)), eval(cc$family))
    family = testfit$family
  } 
  if(!is.null(cc$family) && cc$family == "tobit"){
    family  = tobit()
  }
  if(is.null(cc$family)){
    family=gaussian()
  }
  #family = testfit$family
  
#  famlist = c("binomial", "gaussian", "poisson")
#  if(!(family$family %in% famlist))
#    stop(paste0("Distribution (family) `", family$family, "` not (yet?) supported."))
  
  
  nobs = nrow(data)
  origcall <- thecall <- match.call(expand.dots = FALSE)
  names(thecall) <- gsub("f", "formula", names(thecall))
  m <- match(c("formula", "data", "weights", "offset"), names(thecall), 0L)
  hasweights = ifelse("weights" %in% names(thecall), TRUE, FALSE)
  thecall <- thecall[c(1L, m)]
  thecall$drop.unused.levels <- TRUE
  
  thecall[[1L]] <- quote(stats::model.frame)
  thecalle <- eval(thecall, parent.frame())
  if(hasweights){
    data$weights <- as.vector(model.weights(thecalle))
  } else data$weights = rep(1, nobs)
  
  if (is.null(expnms)) {
    
    #expnms <- attr(terms(f, data = data), "term.labels")
    expnms <- attr(newform, "term.labels")
    
    message("Including all model terms as exposures of interest\n")
  }
  lin = checknames(expnms)
  if(!lin) stop("Model appears to be non-linear and I'm having trouble parsing it:
                  please use `expnms` parameter to define the variables making up the exposure")
  if (!is.null(q) & !is.null(breaks)){
    # if user specifies breaks, prioritize those
    oldq = q
    q <- NULL
  }
  if (!is.null(q) | !is.null(breaks)){
    ql <- quantize(data, expnms, q, breaks)
    qdata <- ql$data
    br <- ql$breaks
    if(is.null(q)){
      # rare scenario with user specified breaks and q is left at NULL
      nvals <- length(br[[1]])-1
    } else{
      nvals <- q
    }
    intvals <- (seq_len(nvals))-1
  } else {
    # if( is.null(breaks) & is.null(q)) # also includes NA
    qdata <- data
    # if no transformation is made (no quantiles, no breaks given)
    # then draw distribution values from quantiles of all the exposures
    # pooled together
    # : allow user specification of this
    nvals = length(table(unlist(data[,expnms])))
    if(nvals < 10){
      if(verbose) message("\nNote: using all possible values of exposure as the
              intervention values\n")
      p = length(expnms)
      intvals <- as.numeric(names(table(unlist(data[,expnms]))))
      
      br <- lapply(seq_len(p), function(x) c(-1e16, intvals[2:nvals]-1e-16, 1e16))
    }else{
      if(verbose) message("\nNote: using quantiles of all exposures combined in order to set
          proposed intervention values for overall effect (25th, 50th, 75th %ile)
        You can ensure this is valid by scaling all variables in expnms to have similar ranges.")
      intvals = as.numeric(quantile(unlist(data[,expnms]), c(.25, .5, .75)))
      br <- NULL
    }
  }
  if(is.null(id)) {
    id <- "id__"
    qdata$id__ <- seq_len(dim(qdata)[1])
  }
  if(is.null(offset)) {
    offset <- "offset__"
    qdata$offset__ <- rep(0,dim(qdata)[1])
  }
  
  
  
  # simultaneous estimation of conditional and 
  basevars = get_all_vars(newform, data=qdata) # variables before being processed by model.frame
  #modframe = model.frame(newform, data=basevars)
  modframe = model.frame(newform, data=qdata)
  X = model.matrix(newform, modframe)
  Y = model.response(modframe)
  # nesting model.frame within the model.matrix function seems to be necessary to get interaction terms to propogate after setting exposures
#  Xint = as.matrix(do.call(rbind,lapply(intvals, function(x) {modframe[,expnms] = x; model.matrix(newform,model.frame(newform, modframe))})))
  #Xint = as.matrix(do.call(rbind,lapply(intvals, function(x) {modframe[,expnms] = x; model.matrix(newform,modframe)})))
  #Xint = as.matrix(do.call(rbind,lapply(intvals, function(x) {mf2 = modframe; mf2[,expnms] = x; model.matrix(newform,model.frame(newform, data=mf2))}))) # works in non-linear setting
  Xint = as.matrix(model.matrix(newform, do.call(rbind,lapply(intvals, function(x) {mf2 = basevars; mf2[,expnms] = x; model.frame(newform, data=mf2)})))) # works in non-linear setting
  Xmsm = poly(Xint[, expnms[1]], degree=degree, raw=TRUE) # intercept and constant exposure
  if(hasintercept){
    Xmsm = cbind(Xint[,colnames(Xint)[1]], Xmsm)
  }

  # point estimates
  np = ncol(X)
  npmsm = ncol(Xmsm)
  
  startvals <- function(family,Y,X,np,npmsm,offset=0){
    fam = family$family
    res = switch(fam,
           binomial = c(log(mean(Y))-(mean(offset)), rep(0, np-1), log(mean(Y))-(mean(offset)), rep(0, npmsm-1)),
           poisson = c(log(mean(Y))-(mean(offset)), rep(0, np-1), log(mean(Y))-(mean(offset)), rep(0, npmsm-1)),
           tobit = c(rep(0, np+npmsm),1),
           rep(0, np+npmsm)
    )
    res
  }
  parminits = startvals(family,Y,X,np,npmsm)
  #.esteq_qgc(parminits, family=family, Y=Y, X=X, Xint=Xint, Xmsm=Xmsm, weights=qdata$weights, rr=FALSE)
  eqfit <- rootSolve::multiroot(.esteq_qgc, start=parminits, family=family, Y=Y, X=X, Xint=Xint, Xmsm=Xmsm, weights=qdata$weights, rr=rr, offset=qdata$offset__, delta=delta)
  # "bread" of the sandwich covariance
  A = numDeriv::jacobian(func=.esteq_qgc, x=eqfit$root, family=family, Y=Y, X=X, Xint=Xint, Xmsm=Xmsm, weights=qdata$weights,method="Richardson", rr=rr,offset=qdata$offset__, delta=delta)
  #
  
  # "meat" of the sandwich covariance
  uid =   unique(qdata[,id,drop=TRUE])
  psii = lapply(uid, function(x){
    selidx = which(qdata[,id,drop=TRUE] == x);
    .esteq_qgcdf(newform, data=modframe[selidx,,drop=FALSE], theta=eqfit$root, family, intvals, expnms, hasintercept, weights=qdata$weights[selidx], degree=degree, rr=rr,offset=qdata$offset__[selidx], delta=delta) 
    } )
  Bi = lapply(psii, function(x) x%*%t(x))
  n = length(Y)
  B = Bi[[1]]
  for(i in 2:length(Bi)){
    B = B + Bi[[i]]
  }
  
  # sandwich covariance
  ibread = solve(A)
  (fullcovmat = ibread %*% B %*% t(ibread))

  condidx = 1:np
  msmidx = (np+1):(np+npmsm)
  estb <- as.numeric(eqfit$root[msmidx])
  nobs <- dim(qdata)[1]
  covmat <- fullcovmat[msmidx,msmidx,drop=FALSE]
  seb <- sqrt(diag(covmat))
  cnms = c(paste0("psi", seq_len(length(msmidx)-1)))
  if(hasintercept)
    cnms = c("(Intercept)", cnms)
  colnames(covmat) <- rownames(covmat) <- names(estb) <- cnms
  
  tstat <- estb / seb
  df <- nobs - length(attr(terms(f, data = data), "term.labels")) - 1 - degree # df based on obs - gcomp terms - msm terms
  pval <- 2 - 2 * pt(abs(tstat), df = df)
  pvalz <- 2 - 2 * pnorm(abs(tstat))
  ci <- cbind(estb + seb * qnorm(alpha / 2), estb + seb * qnorm(1 - alpha / 2))
  # qgcomp 'weights' not applicable in this setting, generally (i.e. if using this function for non-linearity,
  #   then weights will vary with level of exposure)
  if (!is.null(oldq)){
    q = oldq
  }
  allest = eqfit$root
  names(allest) <- c(colnames(X), cnms)
  rownames(fullcovmat) <- colnames(fullcovmat) <- names(allest)
  
  msmfamily = family
  if(msmfamily$family == "binomial" & rr == TRUE){
    msmfamily = binomial(link="log")
  }
    
  fit = list(formula = newform, est=allest[condidx], vcov=fullcovmat[condidx,condidx], family=family, type="conditional")
  msmfit = list(est=allest[msmidx], vcov=fullcovmat[msmidx,msmidx], family=msmfamily, type="msm")
  if(includeX){
    fit$X = X
    msmfit$X = Xmsm
  }
  
  attr(fit, "class") <- c("eefit", attr(fit, "class"))
  attr(msmfit, "class") <- c("eefit", attr(msmfit, "class"))
  
  psiidx = seq_len(degree)+ifelse(hasintercept,1,0)
  qx <- qdata[, expnms]
  names(qx) <- paste0(names(qx), "_q")
  res <- .qgcomp_object(
    qx = qx, fit=fit, msmfit = msmfit,
    psi = estb[psiidx], var.psi = seb[psiidx] ^ 2, 
    covmat.psi=covmat[psiidx,psiidx, drop=FALSE], 
    ci = ci[psiidx,],
    coef = estb, var.coef = seb ^ 2, covmat.coef=covmat, ci.coef = ci,
    expnms=expnms, q=q, breaks=br, degree=degree,
    bootstrap=FALSE,
    y.expected = fit$family$linkinv(Xint %*% coef(fit)), index=Xint[,expnms[1]], # predictions from conditional fit at intervention data
    y.expectedmsm=predict(msmfit),
    bread = A, meat = B,
    covmat.all_robust = fullcovmat,
    alpha=alpha,
    call=origcall,
    hasintercept=hasintercept
  )
  if(msmfit$family$family=='gaussian'){
    res$tstat <- tstat
    res$df <- df
    res$pval <- pval
  }
  else{
    res$zstat <- tstat
    res$pval <- pvalz
  }
  attr(res, "class") <- c("eeqgcompfit", attr(res, "class"))
  res
}


#' @rdname qgcomp.glm.ee
#' @export
qgcomp.ee <- qgcomp.glm.ee
