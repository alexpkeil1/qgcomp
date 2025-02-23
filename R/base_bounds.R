

.pointwise.multinomial <- function(q, py, se.diff, alpha, pwr){
  stop("Not yet implemented")
  # mean, mean differences
  data.frame(quantile= (seq_len(q)) - 1,
             quantile.midpoint=((seq_len(q)) - 1 + 0.5)/(q),
             hx = py,
             mean.diff = py - py[pwr],
             se.diff = se.diff, # standard error on link scale
             ll.diff =  py - py[pwr] + qnorm(alpha/2) * se.diff,
             ul.diff =  py - py[pwr] + qnorm(1-alpha/2) * se.diff
  )
}



.logit <- function(p) log(p) - log(1-p)
.expit <- function(mu) 1/(1+exp(-mu))
.safelog <- function(x,eps=1e-10) ifelse(x==0, log(x+eps), log(x))

# confidence intervals for pointwise comparisons #
.pointwise.lin <- function(q, py, se.diff, alpha, pwr){
  # mean, mean differences
  pw = data.frame(quantile= (seq_len(q)) - 1,
             quantile.midpoint=((seq_len(q)) - 1 + 0.5)/(q),
             hx = py,
             mean.diff = py - py[pwr],
             se.diff = se.diff, # standard error on link scale
             ll.diff =  py - py[pwr] + qnorm(alpha/2) * se.diff,
             ul.diff =  py - py[pwr] + qnorm(1-alpha/2) * se.diff
  )
  pw$ll.linpred = pw$hx - pw$mean.diff + pw$ll.diff
  pw$ul.linpred = pw$hx - pw$mean.diff + pw$ul.diff
  pw
}



.pointwise.log <- function(q, py, se.diff, alpha, pwr){
  # risk, risk ratios / prevalence ratios
  pw = data.frame(quantile= (seq_len(q)) - 1,
             quantile.midpoint=((seq_len(q)) - 1 + 0.5)/(q),
             hx = py,
             rr = exp(py - py[pwr]),
             se.lnrr = se.diff, # standard error on link scale
             ll.rr = exp(py - py[pwr] + qnorm(alpha/2) * se.diff),
             ul.rr = exp(py - py[pwr] + qnorm(1-alpha/2) * se.diff)
  )
  pw$ll.linpred = pw$hx - log(pw$rr) + log(pw$ll.rr)
  pw$ul.linpred = pw$hx - log(pw$rr) + log(pw$ul.rr)
  pw
}



.pointwise.logit <- function(q, py, se.diff, alpha, pwr){
  # odds, odds ratios
  pw = data.frame(quantile= (seq_len(q)) - 1,
             quantile.midpoint=((seq_len(q)) - 1 + 0.5)/(q),
             hx = py, # log odds
             or = exp(py - py[pwr]),
             se.lnor = se.diff, # standard error on link scale
             ll.or = exp(py - py[pwr] + qnorm(alpha/2) * se.diff),
             ul.or = exp(py - py[pwr] + qnorm(1-alpha/2) * se.diff)
  )
  pw$ll.linpred = pw$hx - log(pw$or) + log(pw$ll.or)
  pw$ul.linpred = pw$hx - log(pw$or) + log(pw$ul.or)
  pw
}


.pointwise.zi <- function(q, py, se.diff, alpha, pwr, bootY){
  # marginal risk, risk ratios / prevalence ratios
  lby = .safelog(bootY)
  ff = bootY[pwr,]
  rr = py / py[pwr]
  #rrdist = sweep(bootY, 2, ff, "/")
  lrrdist = sweep(lby,2, lby[pwr,], "-")
  rrdist = exp(lrrdist)
  se.lnrr = apply(lrrdist,1, sd)
  data.frame(quantile= (seq_len(q)) - 1,
             quantile.midpoint=((seq_len(q)) - 1 + 0.5)/(q),
             ey = py, # expected (marginal) rate
             rr = rr,
             se.lnrr = se.lnrr, # standard error on link scale
             ul.rr = exp(.safelog(rr) + qnorm(1-alpha/2)*se.lnrr),
             ll.rr = exp(.safelog(rr) +  qnorm(alpha/2)*se.lnrr),
             pctul.rr = apply(rrdist,1,quantile, .975),
             pctll.rr = apply(rrdist,1,quantile, .025)
  )
}

.pointwise.lin.boot <- function(q, py, se.diff, alpha, pwr){
  pw = .pointwise.lin(q, py, se.diff, alpha, pwr)
  pw$ll.linpred = pw$hx - pw$mean.diff + pw$ll.diff
  pw$ul.linpred = pw$hx - pw$mean.diff + pw$ul.diff
  pw
}

.pointwise.log.boot <- function(q, py, se.diff, alpha, pwr){
  pw = .pointwise.log(q, py, se.diff, alpha, pwr)
  pw$ll.linpred = pw$hx - log(pw$rr) + log(pw$ll.rr)
  pw$ul.linpred = pw$hx - log(pw$rr) + log(pw$ul.rr)
  pw
}

.pointwise.logit.boot <- function(q, py, se.diff, alpha, pwr){
  pw = .pointwise.logit(q, py, se.diff, alpha, pwr)
  pw$ll.linpred = pw$hx - log(pw$or) + log(pw$ll.or)
  pw$ul.linpred = pw$hx - log(pw$or) + log(pw$ul.or)
  pw
}

.pointwise.zi.boot <- .pointwise.zi




# confidence intervals for model slopes  #
.modelwise.lin <- function(q, py, se.diff, alpha, ll, ul){
  # mean, mean differences
  data.frame(quantile= (seq_len(q)) - 1,
             quantile.midpoint=((seq_len(q)) - 1 + 0.5)/(q),
             hx = py,
             m = py,
             se.pw = se.diff, # standard error on link scale
             ll.pw = py + qnorm(alpha/2) * se.diff,
             ul.pw = py + qnorm(1-alpha/2) * se.diff,
             ll.simul= ifelse(is.na(ll), NA, ll),
             ul.simul=ifelse(is.na(ul), NA, ul)
  )
}

.modelwise.log <- function(q, py, se.diff, alpha, ll, ul){
  # risk
  data.frame(quantile= (seq_len(q)) - 1,
             quantile.midpoint=((seq_len(q)) - 1 + 0.5)/(q),
             hx = py, # log risk
             r = exp(py), # risk
             se.pw = se.diff, # standard error on link scale
             ll.pw = exp(py + qnorm(alpha/2) * se.diff),# risk
             ul.pw = exp(py + qnorm(1-alpha/2) * se.diff), # risk
             ll.simul= ifelse(is.na(ll), NA, exp(ll)),
             ul.simul=ifelse(is.na(ul), NA, exp(ul))
  )
}

.modelwise.logit <- function(q, py, se.diff, alpha, ll, ul){
  # odds
  data.frame(quantile= (seq_len(q)) - 1,
             quantile.midpoint=((seq_len(q)) - 1 + 0.5)/(q),
             hx = py, # log odds
             o = exp(py), # odds
             se.pw = se.diff, # standard error on link scale
             ll.pw = exp(py + qnorm(alpha/2) * se.diff),
             ul.pw = exp(py + qnorm(1-alpha/2) * se.diff), # odds
             ll.simul= ifelse(is.na(ll), NA, exp(ll)),
             ul.simul=ifelse(is.na(ul), NA, exp(ul))
  )
}

.modelwise.zi <- function(q, py, se.diff, alpha, ll, ul, bootY){
  # marginal expected rates
  lby = .safelog(bootY)
  se.lpy = apply(lby,1, sd)
  data.frame(quantile= (seq_len(q)) - 1,
             quantile.midpoint=((seq_len(q)) - 1 + 0.5)/(q),
             hx = py, # expected (marginal) rate
             r = py, # expected (marginal) rate
             se.pw = se.lpy, # standard error on link scale
             ll.pw = exp(.safelog(py) + qnorm(alpha/2) * se.lpy),
             ul.pw = exp(.safelog(py) + qnorm(1-alpha/2) * se.lpy),
             ll.simul = NA,
             ul.simul = NA
  )
}



############ pointwise boot ############
#' @title Estimating pointwise comparisons for qgcomp.glm.boot objects
#'
#' @description Calculates: expected outcome (on the link scale), mean difference (link scale)
#' and the standard error of the mean difference (link scale) for pointwise comparisons
#'
#' @details The comparison of interest following a qgcomp fit is often comparisons of model
#' predictions at various values of the joint-exposures (e.g. expected outcome at all exposures
#' at the 1st quartile vs. the 3rd quartile). The expected outcome at a given joint exposure,
#' and marginalized over non-exposure covariates (W), is
#' given as E(Y^s|S) = sum_w E_w(Y|S,W)Pr(W) = sum_i E(Y_i|S) where Pr(W) is the emprical distribution of
#' W and S takes on integer values 0 to q-1.
#' Thus, comparisons are of the type
#' E_w(Y|S=s) - E_w(Y|S=s2) where s and s2 are two different values of the joint exposures (e.g. 0 and 2).
#' This function yields E_w(Y|S) as well as E_w(Y|S=s) - E_w(Y|S=p) where s is any value of S and p is
#' the value chosen via "pointwise ref" - e.g. for binomial variables this will equal the risk/
#' prevalence difference at all values of S, with the referent category S=p-1. The standard error
#' of E(Y|S=s) - E(Y|S=p) is calculated from the bootstrap covariance matrix of E_w(Y|S), such that
#' the standard error for E_w(Y|S=s) - E_w(Y|S=p) is given by
#'
#' Var(E_w(Y|S=s)) + - Var(E_w(Y|S=p)) - 2*Cov(E_w(Y|S=p), - E_w(Y|S=s))
#'
#' This is used to create pointwise confidence intervals. Note that this differs slightly from the
#' \code{\link[qgcomp]{pointwisebound.noboot}} function, which estimates the variance of the conditional
#' regression line given by E(Y|S,W=w), where w is a vector of medians of W (i.e. predictions
#' are made at the median value of all covariates).
#'
#' @param x "qgcompfit" object from `qgcomp.glm.boot`,
#' @param alpha alpha level for confidence intervals
#' @param pointwiseref referent quantile (e.g. 1 uses the lowest joint-exposure category as
#' the referent category for calculating all mean differences/standard deviations)
#' @return A data frame containing
#'  \describe{
#'  \item{linpred: }{The linear predictor from the marginal structural model}
#'  \item{rr/or/mean.diff: }{The canonical effect measure (risk ratio/odds ratio/mean difference) for the marginal structural model link}
#'  \item{se....: }{the stndard error of the effect measure}
#'  \item{ul..../ll....: }{Confidence bounds for the effect measure, and bounds centered at the linear predictor (for plotting purposes)}
#' }
#' @seealso \code{\link[qgcomp]{qgcomp.glm.boot}}, \code{\link[qgcomp]{pointwisebound.noboot}}
#' @export
#' @examples
#' set.seed(12)
#' \dontrun{
#' n=100
#' # non-linear model for continuous outcome
#' dat <- data.frame(x1=(x1 <- runif(100)), x2=runif(100), x3=runif(100), z=runif(100),
#'                   y=runif(100)+x1+x1^2)
#' ft <- qgcomp.glm.boot(y ~ z + x1 + x2 + x3, expnms=c('x1','x2','x3'), data=dat, q=10)
#' pointwisebound.boot(ft, alpha=0.05, pointwiseref=3)
#' }
pointwisebound.boot <- function(x, pointwiseref=1, alpha=0.05){
  if(!x$bootstrap || inherits(x, "survqgcompfit")){
    stop("This function does not work with this type of qgcomp fit")
  }
  if(inherits(x, "qgcompmultfit")){
    stop("Not yet implemented")
  }
  #link = x$fit$family$link
  link = x$msmfit$family$link
  if(inherits(x, "ziqgcompfit")){
    link = "zi"
    npsi = 2*length(coef(x)[[1]])
    maxcidx=1
    for(modtype in names(x$psi)){
      cidx = grep(paste0("^",modtype), names(unlist(x$msmfit$coefficients)))
      maxcidx = max(cidx, maxcidx)
    }
    bootY = x$bootsamps[-seq_len(maxcidx),] # rowwise
  }
  pwr = pointwiseref+0 # may break this in the future
  py = tapply(x$y.expectedmsm, x$index, mean)
  ycovmat = x$cov.yhat # bootstrap covariance matrix of E(y|x) from MSM
  pw.diff = c(0,diff(py))
  pw.vars = numeric(length(pw.diff))
  pw.vars[pwr] = 0
  pw.idx = seq_len(length(py))[-pwr]
  for(j in pw.idx){
    yc = ycovmat[c(pwr,j),c(pwr,j)]
    #pw.vars[j] = 2*sum(diag(yc)) - sum(yc) # shortcut to subtracting covariances
    pw.vars[j] = c(1,-1) %*% yc %*% c(1,-1)
  }
  res = switch(link,
               identity = .pointwise.lin.boot(x$q, py, sqrt(pw.vars), alpha, pointwiseref),
               log = .pointwise.log.boot(x$q, log(py), sqrt(pw.vars), alpha, pointwiseref),
               logit = .pointwise.logit.boot(x$q, .logit(py), sqrt(pw.vars), alpha, pointwiseref),
               zi = .pointwise.zi.boot(x$q, py, NULL, alpha, pwr, bootY)
  )
  fix = which(names(res)=="hx")
  names(res)[fix] = "linpred"
  attr(res, "link") = link
  res
}


############ pointwise noboot ############
#' @title Estimating pointwise comparisons for qgcomp.glm.noboot objects
#'
#' @description Calculates: expected outcome (on the link scale), mean difference (link scale)
#' and the standard error of the mean difference (link scale) for pointwise comparisons
#'
#' @details The comparison of interest following a qgcomp fit is often comparisons of model
#' predictions at various values of the joint-exposures (e.g. expected outcome at all exposures
#' at the 1st quartile vs. the 3rd quartile). The expected outcome at a given joint exposure
#' and at a given level of non-exposure covariates (W=w) is
#' given as E(Y|S,W=w), where S takes on integer values 0 to q-1. Thus, comparisons are of the type
#' E(Y|S=s,W=w) - E(Y|S=s2,W=w) where s and s2 are two different values of the joint exposures (e.g. 0 and 2).
#' This function yields E(Y|S,W=w) as well as E(Y|S=s,W=w) - E(Y|S=p,W=w) where s is any value of S and p is
#' the value chosen via "pointwise ref" - e.g. for binomial variables this will equal the risk/
#' prevalence difference at all values of S, with the referent category S=p-1. For the non-boostrapped
#' version of quantile g-computation (under a linear model). Note that w is taken to be the referent level of covariates
#' so that if meaningful values of E(Y|S,W=w) and E(Y|S=s,W=w) - E(Y|S=p,W=w) are desired,
#' then it is advisable to set the referent levels of W to meaningful values. This
#' can be done by, e.g. centering continuous age so that the predictions are made
#' at the population mean age, rather than age 0.
#'
#' Note that function only works with standard "qgcompfit" objects from `qgcomp.glm.noboot`
#' or `qgcomp.glm.ee` (so it doesn't work
#' with zero inflated, hurdle, or Cox models)
#'
#' Variance for the overall effect estimate is given by:
#' \eqn{transpose(G) Cov(\beta)  G}
#'
#' Where the "gradient vector" G is given by
#' \deqn{G = [\partial(f(\beta))/\partial\beta_1 = 1,
#'   ...,
#'   \partial(f(\beta))/\partial\beta_3k= 1]
#'   }
#' \eqn{f(\beta) = \sum_i^p \beta_i},  and \eqn{\partial y/ \partial x} denotes
#' the partial derivative/gradient. The vector G takes on values that equal the difference in quantiles of
#' S for each pointwise comparison (e.g. for a comparison of the 3rd vs the 5th category,
#' G is a vector of 2s)
#'
#' This variance is used to create pointwise confidence intervals via a normal approximation:
#' (e.g. upper 95% CI = psi + variance*1.96)
#'
#' @param x "qgcompfit" object from `qgcomp.glm.noboot`,
#' @param alpha alpha level for confidence intervals
#' @param pointwiseref referent quantile (e.g. 1 uses the lowest joint-exposure category as
#' the referent category for calculating all mean differences/standard deviations)
#' @return A data frame containing
#'  \describe{
#'  \item{hx: }{The "partial" linear predictor \eqn{\beta_0 + \psi\sum_j X_j^q w_j}, or the effect of the mixture + intercept after
#'  conditioning out any confounders. This is similar to the h(x) function in bkmr. This is not a full prediction of the outcome, but
#'  only the partial outcome due to the intercept and the confounders}
#'  \item{rr/or/mean.diff: }{The canonical effect measure (risk ratio/odds ratio/mean difference) for the marginal structural model link}
#'  \item{se....: }{the stndard error of the effect measure}
#'  \item{ul..../ll....: }{Confidence bounds for the effect measure}
#' }
#' @seealso \code{\link[qgcomp]{qgcomp.glm.noboot}}, \code{\link[qgcomp]{pointwisebound.boot}}
#' @export
#' @examples
#' set.seed(12)
#' \dontrun{
#' n = 100
#' dat <- data.frame(x1=(x1 <- runif(n)), x2=(x2 <- runif(n)), 
#'                   x3=(x3 <- runif(n)), z=(z <- runif(n)),
#'                   y=rnorm(n)+x1 + x2 - x3 +z)
#' # linear model for continuous outcome
#' ft <- qgcomp.glm.noboot(y ~ z + x1 + x2 + x3, 
#'        expnms=c('x1','x2','x3'), data=dat, q=10)
#' ft2 <- qgcomp.glm.boot(y ~ z + x1 + x2 + x3, 
#'         expnms=c('x1','x2','x3'), data=dat, q=10)
#' pointwisebound.noboot(ft, alpha=0.05, pointwiseref=3)
#' pointwisebound.boot(ft2, alpha=0.05, pointwiseref=3)
#' dat <- data.frame(x1=(x1 <- runif(n)), x2=(x2 <- runif(n)), 
#'                    x3=(x3 <- runif(n)), z=(z <- runif(n)),
#'                   y=rbinom(n, 1, 1/(1+exp(-(x1 + x2 - x3 +z)))))
#' # glms for binary outcome, centering covariate to a potentially more meaningful value
#' dat$zcen = dat$z - mean(dat$z)
#' ft <- qgcomp.glm.noboot(y ~ zcen + x1 + x2 + x3, 
#'         expnms=c('x1','x2','x3'), data=dat, q=10, family=binomial())
#' ft2 <- qgcomp.glm.boot(y ~ zcen + x1 + x2 + x3, 
#'         expnms=c('x1','x2','x3'), data=dat, q=10, family=binomial())
#' pointwisebound.noboot(ft, alpha=0.05, pointwiseref=3)
#' pointwisebound.boot(ft2, alpha=0.05, pointwiseref=3)
#' dat$z = as.factor(sample(1:3, n, replace=TRUE))
#' ftf <- qgcomp.glm.noboot(y ~ zcen + x1 + x2 + x3, 
#'        expnms=c('x1','x2','x3'), data=dat, q=10, family=binomial())
#' pointwisebound.noboot(ftf, alpha=0.05, pointwiseref=3)
#' }
pointwisebound.noboot <- function(x, alpha=0.05, pointwiseref=1){
  if(x$bootstrap || inherits(x, "ziqgcompfit") || inherits(x, "survqgcompfit")){
    stop("This function is only for qgcomp.glm.noboot objects")
  }
    #sediff <- function(qdiffj){
  #  # gives standard deviation of pointwise difference on log scale (conditional)
  #  se_comb(x$expnms, vcov(x$fit), grad = qdiffj)
  #}
  sediff <- function(qdiffj, labcolon=""){
    # gives standard deviation of pointwise difference on log scale (conditional)
    se_comb(paste0(labcolon,x$expnms), vcov(x$fit), grad = qdiffj)
  }
  q = x$q
  if(is.null(q))    
    stop("This function only works if q is not NULL")

  link = x$fit$family$link
  #link = x$msmfit$family$link
  #coefs = x$fit$coefficients
  coefs = coef(x$fit)
  px = length(x$expnms)
  diffs = c(rev(seq_len(q-1)), 0:(q-1))
  qdiff = diffs[(q:(q*2-1)) - pointwiseref+1]
  
  if(inherits(x, "qgcompmultfit")){
    stop("Not yet implemented")
    # sketch of code for multfits
    zvars = coefs[-which(colnames(coefs) %in% c("(Intercept)", x$expnms)),]
    msmdesign = cbind(rep(2:(x$nlevels+1), each=q), 0:(q-1))
    labs = x$labs
    cf = coef(x)
    lbp = log(.get_baseline_prob_multinom(x))
    py = numeric(nrow(msmdesign))
    for(r in seq_along(length(py))){
      designrow = msmdesign[r,]
      lab = labs[designrow[1]]
      whichco = paste0(lab, c(".(Intercept)", ".psi"))
      logrr = cf[whichco] %*% c(1, designrow[2])
      py[r] = lbp + logrr
    }
    # py
    se.diffl = lapply(labs[-1], function(l) vapply(qdiff, sediff, 0,labcolon=paste0(l,":")))
    se.diff = do.call(c,se.diffl)
    #cbind(msmdesign, py, rep(qdiff, length(labs[-1])), se.diff)
    
  } else{
    zvars = coefs[-which(names(coefs) %in% c("(Intercept)", x$expnms))]
    xtest = 0:(q-1)
    if(x$hasintercept){
      xtest = cbind(rep(1, q), xtest)
    }
    if(!is.null(x$degree) && x$degree>1){
      extra = do.call(cbind,lapply(2:x$degree, function(x) (0:(q-1))^x))
      xtest = cbind(xtest, extra)
    }
    py = xtest %*% c(coef(x))
    se.diff = vapply(qdiff, sediff, 0)
  }
  res = switch(link,
               identity = .pointwise.lin(q, py, se.diff, alpha, pointwiseref),
               log = .pointwise.log(q, py, se.diff, alpha, pointwiseref),
               logit = .pointwise.logit(q, py, se.diff, alpha, pointwiseref)
  )
  fix = which(names(res)=="hx")
  names(res)[fix] = "linpred"
  attr(res, "link") = link
  res
}


############ model boot ############
#' @title Estimating qgcomp regression line confidence bounds
#'
#' @description Calculates: expected outcome (on the link scale), and upper and lower
#'  confidence intervals (both pointwise and simultaneous)
#'
#' @details This method leverages the bootstrap distribution of qgcomp model coefficients
#' to estimate pointwise regression line confidence bounds. These are defined as the bounds
#' that, for each value of the independent variable X (here, X is the joint exposure quantiles)
#' the 95% bounds (for example) for the model estimate of the regression line E(Y|X) are expected to include the
#' true value of E(Y|X) in 95% of studies. The "simultaneous" bounds are also calculated, and the 95%
#' simultaneous bounds contain the true value of E(Y|X) for all values of X in 95% of studies. The
#' latter are more conservative and account for the multiple testing implied by the former. Pointwise
#' bounds are calculated via the standard error for the estimates of E(Y|X), while the simultaneous
#' bounds are estimated using the bootstrap method of Cheng (reference below). All bounds are large
#' sample bounds that assume normality and thus will be underconservative in small samples. These
#' bounds may also inclue illogical values (e.g. values less than 0 for a dichotomous outcome) and
#' should be interpreted cautiously in small samples.
#'
#'
#' Reference:
#'
#' Cheng, Russell CH. "Bootstrapping simultaneous confidence bands."
#' Proceedings of the Winter Simulation Conference, 2005.. IEEE, 2005.
#'
#' @param x "qgcompfit" object from `qgcomp.glm.boot`,
#' @param alpha alpha level for confidence intervals
#' @param pwonly logical: return only pointwise estimates (suppress simultaneous estimates)
#' @return A data frame containing
#'  \describe{
#'  \item{linpred: }{The linear predictor from the marginal structural model}
#'  \item{r/o/m: }{The canonical measure (risk/odds/mean) for the marginal structural model link}
#'  \item{se....: }{the stndard error of linpred}
#'  \item{ul..../ll....: }{Confidence bounds for the effect measure, and bounds centered at the canonical measure (for plotting purposes)}
#' }
#' The confidence bounds are either  "pointwise" (pw) and "simultaneous" (simul) confidence
#' intervals at each each quantized value of all exposures.
#' @seealso \code{\link[qgcomp]{qgcomp.glm.boot}}
#' @export
#' @examples
#' set.seed(12)
#' \dontrun{
#' dat <- data.frame(x1=(x1 <- runif(50)), x2=runif(50), x3=runif(50), z=runif(50),
#'                   y=runif(50)+x1+x1^2)
#' ft <- qgcomp.glm.boot(y ~ z + x1 + x2 + x3, expnms=c('x1','x2','x3'), data=dat, q=5)
#' modelbound.boot(ft, 0.05)
#' }
modelbound.boot <- function(x, alpha=0.05, pwonly=FALSE){
  if(!x$bootstrap || inherits(x, "survqgcompfit")){
    stop("This function does not work with this type of qgcomp fit")
  }
  if(inherits(x, "qgcompmultfit")){
    stop("Not yet implemented")
  }
  link = x$msmfit$family$link
  if( inherits(x, "ziqgcompfit")){
    pwonly=TRUE
    link = "zi"
    maxcidx=1
    for(modtype in names(x$psi)){
      cidx = grep(paste0("^",modtype), names(unlist(x$msmfit$coefficients)))
      maxcidx = max(cidx, maxcidx)
    }
    bootY = x$bootsamps[-seq_len(maxcidx),] # rowwise
  }
  py = tapply(x$y.expectedmsm, x$index, mean)
  ycovmat = x$cov.yhat # bootstrap covariance matrix of E(y|x) from MSM
  pw.vars = diag(ycovmat)

  if(!pwonly){
    coef.boot = x$bootsamps
    boot.err = t(coef.boot - x$coef)
    iV = solve(qr(x$covmat.coef, tol=1e-20))
    #chi.boots = sapply(seq_len(nrow(boot.err)), function(i) boot.err[i,] %*% iV %*% boot.err[i,])
    chi.boots = vapply(seq_len(nrow(boot.err)), function(i) boot.err[i,] %*% iV %*% boot.err[i,], 0.0)
    chicrit = qchisq(1-alpha, length(x$coef))
    C.set = t(coef.boot[,which(chi.boots<chicrit)])
    fx = function(coef, q=x$q, degree=x$degree){
      reps = numeric(q)
      for(q_ in 0:(q-1)){
        cd = 1
        iv = integer(0)
        while(length(iv)<degree){
          iv = c(iv, q_^cd)
          cd = cd + 1
        }
        reps[q_+1] = c(1, iv) %*% coef
      }
      reps
    }
    fullset = apply(C.set, 1, fx)
    ll = apply(fullset, 1, min)
    ul = apply(fullset, 1, max)
  } else{
    ll=NA
    ul=NA
  }
  res = switch(link,
               identity = .modelwise.lin(x$q, py, sqrt(pw.vars), alpha, ll, ul),
               log = .modelwise.log(x$q, log(py), sqrt(pw.vars), alpha, ll, ul),
               logit = .modelwise.logit(x$q, .logit(py), sqrt(pw.vars), alpha, ll, ul),
               zi = .modelwise.zi(x$q, py, NULL, alpha, ll, ul, bootY)
  )
  fix = which(names(res)=="hx")
  names(res)[fix] = "linpred"
  attr(res, "link") = link
  res
}




modelbound.boot_old <- function(x, alpha=0.05, pwonly=FALSE){
  if(!x$bootstrap) stop("This function is only for qgcomp.glm.boot objects")
  #x = qgcomp.glm.boot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=4,family=gaussian(), B=10000, parallel = TRUE)
  py = tapply(x$y.expectedmsm, x$index, mean)
  ycovmat = x$cov.yhat # bootstrap covariance matrix of E(y|x) from MSM
  pw.vars = diag(ycovmat)

  if(!pwonly){
    coef.boot = x$bootsamps
    boot.err = t(coef.boot - x$coef)
    iV = solve(qr(x$covmat.coef, tol=1e-20))
    #chi.boots = sapply(seq_len(nrow(boot.err)), function(i) boot.err[i,] %*% iV %*% boot.err[i,])
    chi.boots = vapply(seq_len(nrow(boot.err)), function(i) boot.err[i,] %*% iV %*% boot.err[i,], 0.0)
    chicrit = qchisq(1-alpha, length(x$coef))
    C.set = t(coef.boot[,which(chi.boots<chicrit)])
    fx = function(coef, q=x$q, degree=x$degree){
      reps = numeric(q)
      for(q_ in 0:(q-1)){
        cd = 1
        iv = integer(0)
        while(length(iv)<degree){
          iv = c(iv, q_^cd)
          cd = cd + 1
        }
        reps[q_+1] = c(1, iv) %*% coef
      }
      reps
    }
    fullset = apply(C.set, 1, fx)
    ll = apply(fullset, 1, min)
    ul = apply(fullset, 1, max)
  } else{
    ll=NA
    ul=NA
  }
  data.frame(quantile = seq_len(x$q)-1,
             quantile.midpoint=(seq_len(x$q)-1+0.5)/(x$q),
             y.expected = py,
             se.pw=sqrt(pw.vars),
             ll.pw=py + qnorm(alpha/2)*sqrt(pw.vars),
             ul.pw=py + qnorm(1-alpha/2)*sqrt(pw.vars),
             ll.simul=ll,
             ul.simul=ul
  )
}


############ old functions ############
pointwisebound.boot_old <- function(x, alpha=0.05, pointwiseref=1){
  if(!x$bootstrap) stop("This function is only for qgcomp.glm.boot objects")
  # : error catching for other functions

  pwr = pointwiseref+0 # may break this in the future
  py = tapply(x$y.expectedmsm, x$index, mean)
  ycovmat = x$cov.yhat # bootstrap covariance matrix of E(y|x) from MSM
  pw.diff = c(0,diff(py))
  pw.vars = numeric(length(pw.diff))
  pw.vars[pwr] = 0
  pw.idx = seq_len(length(py))[-pwr]
  for(j in pw.idx){
    yc = ycovmat[c(pwr,j),c(pwr,j)]
    #pw.vars[j] = 2*sum(diag(yc)) - sum(yc) # shortcut to subtracting covariances
    pw.vars[j] = c(1,-1) %*% yc %*% c(1,-1)
  }
  data.frame(quantile = seq_len(x$q)-1,
             quantile.midpoint=(seq_len(x$q)-1+0.5)/(x$q),
             y.expected = py,
             mean.diff=py-py[pwr],
             se.diff=sqrt(pw.vars),
             ul.pw = py + qnorm(1-alpha/2)*sqrt(pw.vars),
             ll.pw = py + qnorm(alpha/2)*sqrt(pw.vars)
  )
}



pointwisebound.noboot_old <- function(x, alpha=0.05, pointwiseref = 1){
  if(is.null(x$fit$family$family) || x$fit$family$family == "cox"){
    stop("pointwisebound.noboot not implemented for this method")
  }
  q = x$q
  link = x$fit$family$link
  coefs = x$fit$coefficients
  modmat = model.matrix(x$fit)
  zvars = coefs[-which(names(coefs) %in% c("(Intercept)", x$expnms))]
  av = apply(modmat[,names(zvars), drop=FALSE],2, median)
  # predict at median of non exposure variables
  py = cbind(rep(1, q), 0:(q-1), matrix(rep(av, q), byrow=TRUE, nrow = q)) %*% c(coef(x), zvars)
  px = length(x$expnms)
  sediff <- function(qdiff){
    se_comb(x$expnms, vcov(x$fit), grad = rep(qdiff, px))
  }
  diffs = c(rev(seq_len(q-1)), seq_len(q)-1)
  qdiff = diffs[(q:(q*2-1)) - pointwiseref+1]
  #se.diff = sapply(qdiff, sediff)
  se.diff = vapply(qdiff, sediff, 0.0)
  data.frame(quantile= (seq_len(x$q)) - 1,
             quantile.midpoint=(seq_len(x$q) - 1 + 0.5)/(x$q),
             #y.expected = sapply(seq_len(length(py)), function(i) switch(link,
             y.expected = vapply(seq_len(length(py)), function(i) switch(link,
                                                                  identity = (py[i]),
                                                                  log=exp(py[i]),
                                                                  logit=(1/(1+exp(-py[i])))), 0.0),
             #mean.diff = sapply(seq_len(length(py)), function(i) switch(link,
             mean.diff = vapply(seq_len(length(py)), function(i) switch(link,
                                                                 identity = (py[i] - py[pointwiseref]),
                                                                 log=exp(py[i]) - exp(py[pointwiseref]),
                                                                 logit=(1/(1+exp(-(py[i])))) - (1/(1+exp(-(py[pointwiseref]))))), 0.0),
             se.diff = se.diff,
             #ul.pw = sapply(seq_len(length(py)), function(i) switch(link,
             ul.pw = vapply(seq_len(length(py)), function(i) switch(link,
                                                             identity = (py[i] + qnorm(1-alpha/2) * se.diff[i]),
                                                             log=exp((py[i] + qnorm(1-alpha/2) * se.diff[i])),
                                                             logit=(1/(1+exp(-(py[i] + qnorm(1-alpha/2) * se.diff[i]))))), 0.0) ,
             #ll.pw = sapply(seq_len(length(py)), function(i) switch(link,
             ll.pw = vapply(seq_len(length(py)), function(i) switch(link,
                                                             identity = (py[i] + qnorm(alpha/2) * se.diff[i]),
                                                             log=exp((py[i] + qnorm(alpha/2) * se.diff[i])),
                                                             logit=(1/(1+exp(-(py[i] + qnorm(alpha/2) * se.diff[i]))))), 0.0)
  )
}

