#qgcomp_surv.R: quantile g-computation methods for survival analysis

####### Cox survival curves
#' Survival curve data from a qgcomp survival fit
#' 
#' It is often of interest to examine survival curves from qgcomp.cox.boot 
#' models. They can be useful for checking assumptions about how well
#' the marginal structural model conforms to the underlying conditional
#' model, such that the overall fit approximates the non-linearity 
#' in the underlying model. This function will yield survival curves,
#' but no measures of uncertainty.
#'
#' @param x a `qgcompfit` object from \code{\link[qgcomp]{qgcomp.cox.boot}} 
#' @param ... not used
#'
#' @return a list of data.frames:
#' \describe{
#'  \item{**mdfpop**}{Average Survival curve (survival, time) based on marginal structural model, averaged over the population at every quantile of exposure}
#'  \item{**cdfpop**}{Population average survival curve (survival, time) based on the underlying conditional model}
#'  \item{**mdfq**}{Survival curves (survival, time) for each quantile based on marginal structural model}
#'  \item{**cdfq**}{Survival curves (survival, time) for each quantile based on underlying conditional model}
#' }
#' @md
#' @import survival
#' @export
#' @examples
#' set.seed(50)
#' N=200
#' dat <- data.frame(time=(tmg <- pmin(.1,rweibull(N, 10, 0.1))), 
#'                   d=1.0*(tmg<0.1), x1=runif(N), x2=runif(N), z=runif(N))
#' expnms=paste0("x", 1:2)
#' f = survival::Surv(time, d)~x1 + x2
#' (fit1 <- survival::coxph(f, data = dat))
#' (obj <- qgcomp.cox.noboot(f, expnms = expnms, data = dat))
#' \dontrun{
#' (obj2 <- qgcomp.cox.boot(f, expnms = expnms, data = dat, B=10, MCsize=20000))
#' curves = cox.survcurve.boot(obj2)
#' rbind(head(curves$mdfq),tail(curves$mdfq))
#' }
#' 
#'
qgcomp.survcurve.boot <- function(x, ...){
  if(is.null(x$q)){
    q = length(x$breaks[[1]])-1
    warning("q is NULL in the qgcomp fit: be sure that the given breaks between categories are meaningful")
  } else{
    q = x$q
  }
  namespaceImport("survival")
  rootdat <- as.data.frame(x$fit$x)
  psidat <- data.frame(psi=0)
  rootfun <- function(idx, df){
    df[,x$expnms] <- idx
    df
  }
  rootfun2 <- function(idx, df){
    df[,"psi"] <- idx
    df[,"psi1"] <- idx
    df[,"psi2"] <- idx^2
    df[,"psi3"] <- idx^3
    df[,"psi4"] <- idx^4
    df
  }
  newmarg = lapply(0:(q-1), rootfun2, df=psidat)
  margdf = data.frame(do.call("rbind", newmarg))
  newcond = lapply(0:(q-1), rootfun, df=rootdat)
  conddf = data.frame(do.call("rbind", newcond))
  msmobj = survival::survfit(x$msmfit, newdata=margdf)
  gcompobj = survival::survfit(x$fit, newdata=conddf)
  #
  mdfl = lapply(seq_len(q), function(zz) with(survival::survfit(x$msmfit, newdata=newmarg[[zz]]), data.frame(time=time, surv=surv, q=zz)))
  cdfl = lapply(seq_len(q), function(zz) with(survival::survfit(x$fit, newdata=newcond[[zz]][1,]), data.frame(time=time, surv=surv, q=zz)))
  mdfq = do.call(rbind, mdfl)
  cdfq = do.call(rbind, cdfl)
  mdf = with(msmobj, data.frame(time=time, surv=apply(surv, 1, mean)))
  cdf = with(gcompobj, data.frame(time=time, surv=apply(surv, 1, mean)))
  list(
    mdfpop = mdf, # 
    cdfpop = cdf,
    mdfq = mdfq,
    cdfq = cdfq
  )
}


coxmsm_fit <- function(
  f, qdata, intvals, expnms, main=TRUE, degree=1, id=NULL, weights, cluster=NULL, MCsize=10000, ...){
  #' @title Marginal structural Cox model (MSM) fitting within quantile g-computation
  #' @description this is an internal function called by \code{\link[qgcomp]{qgcomp.cox.noboot}},
  #'  \code{\link[qgcomp]{qgcomp.cox.boot}}, and \code{\link[qgcomp]{qgcomp.cox.noboot}},
  #'  but is documented here for clarity. Generally, users will not need to call
  #'  this function directly.
  #' @details This function first computes expected outcomes under hypothetical
  #' interventions to simultaneously set all exposures to a specific quantile. These
  #' predictions are based on g-computation, where the exposures are `quantized',
  #' meaning that they take on ordered integer values according to their ranks,
  #' and the integer values are determined by the number of quantile cutpoints used.
  #' The function then takes these expected outcomes and fits an additional model
  #' (a marginal structural model) with the expected outcomes as the outcome and
  #' the intervention value of the exposures (the quantile integer) as the exposure.
  #' Under causal identification assumptions and correct model specification,
  #' the MSM yields a causal exposure-response representing the incremental
  #' change in the expected outcome given a joint intervention on all exposures.
  #' @param f an R formula representing the conditional model for the outcome, given all
  #' exposures and covariates. Interaction terms that include exposure variables
  #' should be represented via the \code{\link[base]{AsIs}} function. Offset
  #' terms can be included via \code{Surv(time,event) ~ exposure + offset(z)}
  #' @param qdata a data frame with quantized exposures (as well as outcome and other covariates)
  #' @param intvals sequence, the sequence of integer values that the joint exposure 
  #' is 'set' to for estimating the msm. For quantile g-computation, this is just 
  #' 0:(q-1), where q is the number of quantiles of exposure.
  #' @param expnms a character vector with the names of the columns in qdata that represent
  #' the exposures of interest (main terms only!)
  #' @param main logical, internal use: produce estimates of exposure effect (psi)
  #'  and expected outcomes under g-computation and the MSM
  #' @param degree polynomial bases for marginal model (e.g. degree = 2
  #'  allows that the relationship between the whole exposure mixture and the outcome
  #'  is quadratic. Default=1)
  #' @param id (optional) NULL, or variable name indexing individual units of 
  #' observation (only needed if analyzing data with multiple observations per 
  #' id/cluster)
  #' @param weights "case weights" - passed to the "weight" argument of 
  #' \code{\link[survival]{coxph}}
  #' @param cluster not yet implemented
  #' @param MCsize integer: sample size for simulation to approximate marginal 
  #'  hazards ratios
  #' @param ... arguments to coxph (e.g. ties)
  #' @seealso \code{\link[qgcomp]{qgcomp.cox.boot}}, and \code{\link[qgcomp]{qgcomp.cox.noboot}}
  #' @concept variance mixtures
  #' @import survival
  #' @export
  #' @examples
  #' set.seed(50)
  #' dat <- data.frame(time=(tmg <- pmin(.1,rweibull(50, 10, 0.1))), d=1.0*(tmg<0.1), 
  #'                   x1=runif(50), x2=runif(50), z=runif(50))
  #' expnms=paste0("x", 1:2)
  #' qdata  = quantize(dat, expnms)$data
  #' f = survival::Surv(time, d)~x1 + x2
  #' fit <- survival::coxph(f, data = qdata, y=TRUE, x=TRUE)
  #' r1 = qdata[1,,drop=FALSE]
  #' times = survival::survfit(fit, newdata=r1, se.fit=FALSE)$time
  #' (obj <- coxmsm_fit(f, qdata, intvals=c(0,1,2,3), expnms, main=TRUE, degree=1, 
  #'    id=NULL, MCsize=100))
  #' #dat2 <- data.frame(psi=seq(1,4, by=0.1))
  #' #summary(predict(obj))
  #' #summary(predict(obj, newdata=dat2))
  msmweights <- NULL

  newform <- terms(f, data = qdata)
  class(newform) <- "formula"

  nobs = nrow(qdata)
  thecall <- match.call(expand.dots = FALSE)
  names(thecall) <- gsub("f", "formula", names(thecall))
  names(thecall) <- gsub("qdata", "data", names(thecall))
  m <- match(c("formula", "data", "weights", "offset"), names(thecall), 0L)
  #m <- match(c("f", "data", "weights", "offset"), names(thecall), 0L)
  hasweights = ifelse("weights" %in% names(thecall), TRUE, FALSE)
  thecall <- thecall[c(1L, m)]
  thecall$drop.unused.levels <- TRUE
  
  thecall[[1L]] <- quote(stats::model.frame)
  thecalle <- eval(thecall, parent.frame())
  if(hasweights){
    qdata$weights <- as.vector(model.weights(thecalle))
  } else qdata$weights = rep(1, nobs)

  if(is.null(id)) {
    id = "id__"
    qdata$id__ = seq_len(dim(qdata)[1])
  }
  # conditional outcome regression fit
  environment(newform) <- list2env(list(qdata=qdata))
  fit <- coxph(newform, data = qdata, x=TRUE, y=TRUE, 
               weights=weights,# cluster=cluster, 
               ...)
  if(nrow(qdata) != length(fit$linear.predictors)){
    cat("nobs:", nobs)
    cat("; n preds:", length(predmat[[1]]))
    stop("You have missing data. qgcomp.cox.boot only works 
             on complete case data. Observations with missing 
             values of outcomes, exposures or covariates should be removed 
             from the dataset prior to fitting.")
  }
  
  ## get predictions (set exposure to 0,1,...,q-1)
  if(is.null(intvals)){
    intvals = seq_len(length(table(qdata[expnms[1]]))) - 1
  }
  ymat = fit$y
  tval = grep("stop|time",colnames(ymat) , value=TRUE)
  stop = as.numeric(ymat[,tval])
  times = sort(-sort(-unique(stop))[-1])
  newids <- data.frame(temp=sort(sample(unique(qdata[,id, drop=TRUE]), MCsize, 
                                         #probs=weights, #weighted bootstrap sampling works with fixed weights, but not time-varying weights
                                         replace = TRUE
  )))
  names(newids) <- id
  newdata <- merge(qdata,newids, by=id, all.x=FALSE, all.y=TRUE)[seq_len(MCsize),] # this might trim one id, keep MCsize large
  
  predit <- function(idx){
    newdata[,expnms] <- idx
    # predictions under hypothetically removing competing risks
    # assuming censoring at random and no late entry
    pfit = survfit(fit, newdata=newdata, se.fit=FALSE)
    if(any(diff(pfit$n.risk)>0)){
      warning("Late entry/counting process style data is detected 
              in call to qgcomp.cox.boot. Note:
              this function is not valid with time varying
              exposures/covariates in the model if using counting process
              style data (>1 record per individual)."
              )
    }
    ch = pfit$cumhaz
    h1 = ch[1,]
    haz = rbind(h1, apply(ch, 2, diff))
    dh = dim(haz)
    ui = matrix(runif(prod(dh)), nrow=dh[1], ncol=dh[2])
    dall = 1.0*(haz > ui)
    dt1 = apply(dall, 2, which.max)
    dt2 = apply(dall, 2, which.min)
    tidx = ifelse(dt1==1 & dt2==1,length(times),dt1)
    d = ifelse(tidx < length(times),1, 0)
    time = times[tidx]
    mean(d)
    Surv(time,d)
  }
  predmat = lapply(intvals, predit)
  msmdat <- data.frame(
    Ya = do.call("c", predmat),
    psi = rep(intvals, each=MCsize),
    #msmweights = newdata$weights
    msmweights = rep(newdata$weights, times=length(intvals))
  )
  msmforms = paste0("Ya ~ ", 
                    "poly(psi, degree=",degree,", raw=TRUE)"
  )
  msmform = as.formula(msmforms)
  
  msmfit <- coxph(msmform, data=msmdat, x=TRUE, y=TRUE, 
                  weights=msmweights
                  )
  coxfam = list(family='cox', link='log', linkfun=log)
  class(coxfam) = "family"
  res = list(fit=fit, msmfit=msmfit)
  if(main) {
    res$Ya = msmdat$Ya   # expected outcome under joint exposure, by gcomp
    #res$Yamsm = exp(-predict(msmfit, newdata = msmdat[,], type="expected")) # not yet implemented
    res$A =  msmdat$psi # joint exposure (0 = all exposures set category with 
    # upper cut-point as first quantile)
  }
  #res$fit[['family']] = coxfam # kludge for print function
  #res$msmfit[['family']] = coxfam # kludge for print function
  res$fit[['family']] = cox() # kludge for print function
  res$msmfit[['family']] = cox() # kludge for print function
  class(res$msmfit) = c("coxmsmfit",class(res$msmfit))
  res
}

coxmsm.fit = coxmsm_fit

#predict.coxmsmfit <- function(msmfit, newdata=NULL, ...){
#  if(is.null(newdata)){
#    pfit = survfit(msmfit, se.fit=FALSE)
#  } else{
#    pfit = survfit(msmfit, newdata=newdata, se.fit=FALSE)
#  }
#  bh = pfit$cumhaz
#  tms = pfit$time
#  deltat = c(tms[1], diff(tms))
#  haz = rbind(bh[1,], apply(bh, 2, diff))
#  Ya = apply(deltat*haz,2, sum)
#  Ya
#}

qgcomp.cox.noboot <- function (f, data, expnms = NULL, q = 4, breaks = NULL,
                               id=NULL, weights, cluster=NULL, alpha=0.05,...) {
  #' @title Quantile g-computation for survival outcomes under linearity/additivity
  #'  
  #'
  #' @description This function performs quantile g-computation in a survival
  #' setting. The approach estimates the covariate-conditional hazard ratio for 
  #' a joint change of 1 quantile in each exposure variable specified in expnms
  #' parameter
  #' 
  #' @details For survival outcomes (as specified using methods from the 
  #' survival package), this yields a conditional log hazard ratio representing  
  #' a change in the expected conditional hazard (conditional on covariates)
  #' from increasing every exposure by 1 quantile. In general, this quantity 
  #' quantity is not equivalent to g-computation estimates. Hypothesis test
  #' statistics and 95% confidence intervals are based on using the delta
  #' estimate variance of a linear combination of random variables.
  #' 
  #' @param f R style survival formula, which includes \code{\link[survival]{Surv}}
  #'   in the outcome definition. E.g. \code{Surv(time,event) ~ exposure}. Offset
  #'   terms can be included via \code{Surv(time,event) ~ exposure + offset(z)}
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
  #' id/cluster)
  #' @param weights "case weights" - passed to the "weight" argument of 
  #' \code{\link[survival]{coxph}}
  #' @param cluster not yet implemented
  #' @param alpha alpha level for confidence limit calculation
  #' @param ... arguments to coxph
  #' @seealso \code{\link[qgcomp]{qgcomp.cox.boot}}, \code{\link[qgcomp]{qgcomp.glm.boot}}, 
  #'   and \code{\link[qgcomp]{qgcomp}}
  #' @family qgcomp_methods
  #' @return a qgcompfit object, which contains information about the effect
  #'  measure of interest (psi) and associated variance (var.psi), as well
  #'  as information on the model fit (fit) and information on the 
  #'  weights/standardized coefficients in the positive (pos.weights) and 
  #'  negative (neg.weights) directions.
  #' @concept variance mixtures
  #' @import survival
  #' @export
  #' @examples
  #' set.seed(50)
  #' N=200
  #' dat <- data.frame(time=(tmg <- pmin(.1,rweibull(N, 10, 0.1))), 
  #'                 d=1.0*(tmg<0.1), x1=runif(N), x2=runif(N), z=runif(N))
  #' expnms=paste0("x", 1:2)
  #' f = survival::Surv(time, d)~x1 + x2
  #' (fit1 <- survival::coxph(f, data = dat))
  #' (obj <- qgcomp.cox.noboot(f, expnms = expnms, data = dat))
  #' \dontrun{
  #' 
  #' # weighted analysis
  #' dat$w = runif(N)
  #' qdata = quantize(dat, expnms=expnms)
  #' (obj2 <- qgcomp.cox.noboot(f, expnms = expnms, data = dat, weight=w))
  #' obj2$fit
  #' survival::coxph(f, data = qdata$data, weight=w)
  #' 
  #' # not run: bootstrapped version is much slower
  #' (obj2 <- qgcomp.cox.boot(f, expnms = expnms, data = dat, B=200, MCsize=20000))
  #' }

  of <- f
  newform <- terms(f, data = data)
  class(newform) <- "formula"

  nobs = nrow(data)
  origcall <- thecall <- match.call(expand.dots = FALSE)
  names(thecall) <- gsub("f", "formula", names(thecall))
  m <- match(c("f", "data", "weights", "offset"), names(thecall), 0L)
  #m <- match(c("f", "data", "weights", "offset"), names(thecall), 0L)
  hasweights = ifelse("weights" %in% names(thecall), TRUE, FALSE)
  thecall <- thecall[c(1L, m)]
  thecall$drop.unused.levels <- TRUE
  
  thecall[[1L]] <- quote(stats::model.frame)
  thecalle <- eval(thecall, parent.frame())
  if(hasweights){
    data$weights <- as.vector(model.weights(thecalle))
  } else data$weights = rep(1, nobs)

  if (is.null(expnms)) {
    message("Including all model terms as exposures of interest")
  
    #expnms <- attr(terms(f, data = data), "term.labels")
    expnms <- attr(newform, "term.labels")
  
  }
  if (!is.null(q) | !is.null(breaks)) {
    ql <- quantize(data, expnms, q, breaks)
    qdata <- ql$data
    br <- ql$breaks
  }
  else {
    qdata <- data
    br <- breaks
  }
  # original fit
  environment(newform) <- list2env(list(qdata=qdata))
  fit <- coxph(newform, data = qdata,
               weights=weights, 
               #cluster=cluster,
               ...)
  coxfam = list(family='cox', link='log', linkfun=log)
  class(coxfam) = "family"
  fit[['family']] = coxfam # kludge for print function
  mod <- summary(fit)
  estb <- sum(mod$coefficients[expnms, 1])
  covMat = fit$var
  colnames(covMat) <- names(mod$coefficients[, 1])
  seb <- se_comb(expnms, covmat = covMat)
  tstat <- estb/seb
  #df <- mod$df.null - length(expnms)
  #pval <- 2 - 2 * pt(abs(tstat), df = df)
  pvalz <- 2 - 2 * pnorm(abs(tstat))
  ci <- c(estb + seb * qnorm(alpha/2), estb + seb * qnorm(1 - 
                                                            alpha/2))
  wcoef <- fit$coefficients[expnms]
  names(wcoef) <- gsub("_q", "", names(wcoef))
  poscoef <- which(wcoef > 0)
  negcoef <- which(wcoef <= 0)
  pos.weights <- abs(wcoef[poscoef])/sum(abs(wcoef[poscoef]))
  neg.weights <- abs(wcoef[negcoef])/sum(abs(wcoef[negcoef]))
  pos.psi <- sum(wcoef[poscoef])
  neg.psi <- sum(wcoef[negcoef])
  #nmpos = names(pos.weights)
  #nmneg = names(neg.weights)
  #se.pos.psi <- se_comb(nmpos, covmat = covMat)
  #se.neg.psi <- se_comb(nmneg, covmat = covMat)
  qx <- qdata[, expnms]
  names(qx) <- paste0(names(qx), "_q")
  names(estb) = "psi1"
  res <- .qgcomp_object(
    qx = qx, fit = fit, 
    psi = estb, var.psi = seb^2, covmat.psi = seb^2, ci = ci, 
    coef = estb, var.coef = seb^2, covmat.coef = seb^2, ci.coef = ci, 
    expnms = expnms, q = q, breaks = br, degree = 1, 
    pos.psi = pos.psi, neg.psi = neg.psi, 
    pos.weights = sort(pos.weights, decreasing = TRUE), 
    neg.weights = sort(neg.weights, decreasing = TRUE), 
    pos.size = sum(abs(wcoef[poscoef])), neg.size = sum(abs(wcoef[negcoef])), 
    zstat = tstat, pval = pvalz, alpha=alpha, call=origcall, hasintercept=FALSE
  )
  attr(res, "class") <- c("survqgcompfit", attr(res, "class"))
  res
}

qgcomp.cox.boot <- function(f, data, expnms=NULL, q=4, breaks=NULL, 
                                       id=NULL, weights, cluster=NULL, alpha=0.05, B=200, MCsize=10000, degree=1, 
                                       seed=NULL, parallel=FALSE,  parplan = FALSE,
 ...
){# bayes=FALSE,rr=TRUE, 
  #' @title Quantile g-computation for survival outcomes
  #'  
  #' @description This function yields population average effect estimates for 
  #'   (possibly right censored) time-to event outcomes
  #'  
  #' @details `qgcomp.cox.boot' estimates the
  #'  log(hazard ratio) per quantile increase in the joint exposure to all exposures 
  #'  in `expnms'. This function uses g-computation to estimate the parameters of a
  #'  marginal structural model for the population average effect of increasing all
  #'  exposures in `expnms' by a single quantile. This approach involves specifying 
  #'  an underlying conditional outcome model, given all exposures of interest (possibly
  #'  with non-linear basis function representations such as splines or product terms)
  #'  and confounders or covariates of interest. This model is fit first, which is used
  #'  to generate expected outcomes at each quantile of all exposures, which is then
  #'  used in a second model to estimate a population average dose-response curve that
  #'  is linear or follows a simple polynomial function. See section on MCSize below
  #'  
  #'  Test statistics and confidence intervals are based on 
  #'  a non-parametric bootstrap, using the standard deviation of the bootstrap
  #'  estimates to estimate the standard error. The bootstrap standard error is 
  #'  then used to estimate Wald-type confidence intervals. Note that no bootstrapping
  #'  is done on estimated quantiles of exposure, so these are treated as fixed
  #'  quantities
  #'  
  #'  MCSize is crucial to get accurate point estimates. In order to get marginal
  #'  estimates of the population hazard under different values of the joint exposure
  #'  at a given quantile for all exposures in `expnms`, `qgcomp.cox.boot` uses
  #'  Monte Carlo simulation to generate outcomes implied by the underlying conditional model
  #'  and then fit a separate (marginal structural) model to those outcomes. In order to get
  #'  accurate results that don't vary much from run-to-run of this approach, MCsize
  #'  must be set large enough so that results are stable across runs according to a pre-determined
  #'  precision (e.g. 2 significant digits). 
  #'
  #' @param f R style survival formula, which includes \code{\link[survival]{Surv}}
  #'   in the outcome definition. E.g. \code{Surv(time,event) ~ exposure}. Offset
  #'   terms can be included via \code{Surv(time,event) ~ exposure + offset(z)}
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
  #' id/cluster). Note that qgcomp.cox.noboot will not produce cluster-appropriate
  #' standard errors. qgcomp.cox.boot can be used for this, which will use bootstrap
  #' sampling of clusters/individuals to estimate cluster-appropriate standard
  #' errors via bootstrapping.
  #' @param weights "case weights" - passed to the "weight" argument of 
  #' \code{\link[survival]{coxph}}
  #' @param cluster not yet implemented
  #' @param alpha alpha level for confidence limit calculation
  #' @param B integer: number of bootstrap iterations (this should typically be >=200,
  #'  though it is set lower in examples to improve run-time).
  #' @param degree polynomial bases for marginal model (e.g. degree = 2
  #'  allows that the relationship between the whole exposure mixture and the outcome
  #'  is quadratic.
  #' @param MCsize integer: sample size for simulation to approximate marginal 
  #'  hazards ratios (if < sample size, then set to sample size). Note that large
  #'  values will slow down the fitting, but will result in higher accuracy - if you 
  #'  run the function multiple times you will see that results vary due to simulation
  #'  error. Ideally, MCsize would be set such that simulation error is negligible
  #'  in the precision reported (e.g. if you report results to 2 decimal places, then
  #'  MCsize should be set high enough that you consistenty get answers that are the same
  #'  to 2 decimal places).
  #' @param seed integer or NULL: random number seed for replicable bootstrap results
  #' @param parallel logical (default FALSE): use future package to speed up bootstrapping
  #' @param parplan (logical, default=FALSE) automatically set future::plan to plan(multisession) (and set to existing plan, if any, after bootstrapping)
  #' @param ... arguments to coxph
  #' @family qgcomp_methods
  #' @return a qgcompfit object, which contains information about the effect
  #'  measure of interest (psi) and associated variance (var.psi), as well
  #'  as information on the model fit (fit) and information on the 
  #'  marginal structural model (msmfit) used to estimate the final effect
  #'  estimates.
  #' @concept variance mixtures
  #' @import stats survival
  #' @export
  #' @examples
  #' set.seed(50)
  #' N=200
  #' dat <- data.frame(time=(tmg <- pmin(.1,rweibull(N, 10, 0.1))), 
  #'                 d=1.0*(tmg<0.1), x1=runif(N), x2=runif(N), z=runif(N))
  #' expnms=paste0("x", 1:2)
  #' f = survival::Surv(time, d)~x1 + x2
  #' (fit1 <- survival::coxph(f, data = dat))
  #' (obj <- qgcomp.cox.noboot(f, expnms = expnms, data = dat))
  #' \dontrun{
  #' # not run (slow when using boot version to proper precision)
  #' (obj2 <- qgcomp.cox.boot(f, expnms = expnms, data = dat, B=10, MCsize=20000))
  #' 
  #' # weighted analysis
  #' 
  #' # using future package, marginalizing over confounder z
  #' (obj3 <- qgcomp.cox.boot(survival::Surv(time, d)~x1 + x2 + z, expnms = expnms, data = dat, 
  #'                          B=1000, MCsize=20000, parallel=TRUE, parplan=TRUE))
  #' # non-constant hazard ratio, non-linear terms
  #' (obj4 <- qgcomp.cox.boot(survival::Surv(time, d)~factor(x1) + splines::bs(x2) + z, 
  #'                          expnms = expnms, data = dat, 
  #'                          B=1000, MCsize=20000, parallel=FALSE, degree=1))
  #'                          
  #' # weighted analysis
  #' dat$w = runif(N)
  #' (objw1 <- qgcomp.cox.noboot(f, expnms = expnms, data = dat, weights=w))
  #' (objw2 <- qgcomp.cox.boot(f, expnms = expnms, data = dat, weights=w, B=5, MCsize=20000))
  #' }
  requireNamespace("survival")
  
  if(is.null(seed)) seed = round(runif(1, min=0, max=1e8))

  newform <- terms(f, data = data)
  class(newform) <- "formula"

  nobs = nrow(data)
  origcall <- thecall <- match.call(expand.dots = FALSE)
  names(thecall) <- gsub("f", "formula", names(thecall))
  m <- match(c("f", "data", "weights", "offset"), names(thecall), 0L)
  #m <- match(c("f", "data", "weights", "offset"), names(thecall), 0L)
  hasweights = ifelse("weights" %in% names(thecall), TRUE, FALSE)
  thecall <- thecall[c(1L, m)]
  thecall$drop.unused.levels <- TRUE
  
  thecall[[1L]] <- quote(stats::model.frame)
  thecalle <- eval(thecall, parent.frame())
  if(hasweights){
    data$weights <- as.vector(model.weights(thecalle))
  } else data$weights = rep(1, nobs)

  if (is.null(expnms)) {
    message("Including all model terms as exposures of interest")
  
    #expnms <- attr(terms(f, data = data), "term.labels")
    expnms <- attr(newform, "term.labels")
  
  }
  lin = checknames(expnms)
  if(!lin) stop("Model appears to be non-linear and I'm having trouble parsing it: 
                please use `expnms` parameter to define the variables making up the exposure")
  if (!is.null(q) & !is.null(breaks)){
    # if user specifies breaks, prioritize those
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
    intvals <- seq_len(nvals)-1
  } else {
    # if( is.null(breaks) & is.null(q)) # also includes NA
    qdata <- data
    # if no transformation is made (no quantiles, no breaks given)
    # then draw distribution values from quantiles of all the exposures
    # pooled together
    # : allow user specification of this
    message("\nNote: using quantiles of all exposures combined in order to set 
          proposed intervention values for overall effect (25th, 50th, 75th %ile)
        You can ensure this is valid by scaling all variables in expnms to have similar ranges.")
    intvals = as.numeric(quantile(unlist(data[,expnms]), c(.25, .5, .75)))
    br <- NULL
  }
  if(is.null(id)) {
    id <- "id__"
    qdata$id__ <- seq_len(dim(qdata)[1])
  }
  if(dim(qdata)[1]>MCsize) MCsize = dim(qdata)[1]
  ###
  environment(newform) <- list2env(list(qdata=qdata))
  msmfit <- coxmsm_fit(newform, qdata, intvals, expnms, main=TRUE,degree=degree, 
                       weights=weights,
                       # cluster = cluster,
                       id=id, MCsize=MCsize, ...)
  # main estimate  
  estb <- as.numeric(msmfit$msmfit$coefficients)
  #bootstrap to get std. error
  #nobs <- dim(qdata)[1]
  nids <- length(unique(qdata[,id, drop=TRUE]))
  starttime = Sys.time()
  psi.only <- function(i=1, f=f, qdata=qdata, intvals=intvals, expnms=expnms, degree=degree,
                       weights=weights, #cluster = cluster,
                       nids=nids, id=id, MCsize=MCsize, ...){
    if(i==2 & !parallel){
      timeiter = as.numeric(Sys.time() - starttime)
      if((timeiter*B/60)>0.5) message(paste0("Expected time to finish: ", round(B*timeiter/60, 2), " minutes \n"))
    }
    bootids <- data.frame(temp=sort(sample(unique(qdata[,id, drop=TRUE]), nids, replace = TRUE)))
    names(bootids) <- id
    qdata_ <- merge(qdata,bootids, by=id, all.x=FALSE, all.y=TRUE)
    newft = coxmsm_fit(f, qdata_, intvals, expnms, main=FALSE, degree, id, 
                       weights=weights, #cluster = cluster,
                       MCsize=MCsize, ...)
    as.numeric(newft$msmfit$coefficients)
  }
  set.seed(seed)
  if(parallel){
    #Sys.setenv(R_FUTURE_SUPPORTSMULTICORE_UNSTABLE="quiet")
    if (parplan) {
        oplan <- future::plan(strategy = future::multisession)
        on.exit(future::plan(oplan), add = TRUE)
      }
    bootsamps <- future.apply::future_lapply(X=seq_len(B), FUN=psi.only,
                                    f=newform, qdata=qdata, intvals=intvals, 
                                    expnms=expnms, degree=degree, nids=nids, id=id,
                                    future.seed=TRUE,
                                    weights=qdata$weights, MCsize=MCsize, ...)
    
  }else {
    bootsamps <- lapply(X=seq_len(B), FUN=psi.only,
                        f=newform, qdata=qdata, intvals=intvals, 
                        expnms=expnms, degree=degree, nids=nids, id=id, 
                        weights=weights, MCsize=MCsize, ...)
  }
  bootsamps = do.call("cbind", bootsamps)
  if(is.null(dim(bootsamps))) {
    seb <- sd(bootsamps)
    covmat <- var(bootsamps)
    names(covmat) <- 'psi1'
  }else{
    seb <- apply(bootsamps, 1, sd)
    covmat <- cov(t(bootsamps))
    colnames(covmat) <- rownames(covmat) <- names(estb) <- paste0("psi", seq_len(nrow(bootsamps)))
  }
  tstat <- estb / seb
  pvalz <- 2 - 2 * pnorm(abs(tstat))
  ci <- cbind(estb + seb * qnorm(alpha / 2), estb + seb * qnorm(1 - alpha / 2))
  # 'weights' not applicable in this setting, generally (i.e. if using this function 
  #   for non-linearity, then weights will vary with level of exposure)
  qx <- qdata[, expnms]
  res <- .qgcomp_object(
    qx = qx, fit = msmfit$fit, msmfit = msmfit$msmfit, 
    psi = estb, var.psi = seb ^ 2, covmat.psi=covmat, ci = ci,
    coef = estb, var.coef = seb ^ 2, covmat.coef=covmat, ci.coef = ci, 
    expnms=expnms, q=q, breaks=br, degree=degree,
    bootstrap=TRUE,
    y.expected=msmfit$Ya, 
    y.expectedmsm=msmfit$Yamsm, 
    index=msmfit$A,
    bootsamps = bootsamps,
    alpha=alpha, call=origcall, hasintercept=FALSE
  )
  if(msmfit$fit$family$family=='cox'){
    res$zstat <- tstat
    res$pval <- pvalz
  } else{
    stop("MSM fit is not a cox model, which is an unexpected bug. 
         Send any relevant info to akeil@unc.edu")
  }
  attr(res, "class") <- c("survqgcompfit", attr(res, "class"))
  res
}
