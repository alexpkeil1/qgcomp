#qgcomp_surv.R: quantile g-computation methods for survival analysis

qgcomp.cox.noboot <- function (f, data, expnms = NULL, q = 4, breaks = NULL,
                               id=NULL, alpha=0.05,
          ...) {
  #' @title estimation of quantile g-computation fit for a survival outcome
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
  #' id/cluster)
  #' @param alpha alpha level for confidence limit calculation
  #' @param ... arguments to glm (e.g. family)
  #' @seealso \code{\link[qgcomp]{qgcomp.boot}}, and \code{\link[qgcomp]{qgcomp}}
  #' @return a qgcompfit object, which contains information about the effect
  #'  measure of interest (psi) and associated variance (var.psi), as well
  #'  as information on the model fit (fit) and information on the 
  #'  weights/standardized coefficients in the positive (pweights) and 
  #'  negative (nweight) directions.
  #' @keywords variance, mixtures
  #' @import survival
  #' @export
  #' @examples
  #' runif(1)
  if (is.null(expnms)) {
    cat("Including all model terms as exposures of interest")
    expnms <- attr(terms(f, data = data), "term.labels")
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
  fit <- coxph(f, data = qdata, ...)
  coxfam = list(family='cox', link='log', linkfun=log)
  class(coxfam) = "family"
  fit[['family']] = coxfam # kludge for print function
  mod <- summary(fit)
  estb <- sum(mod$coefficients[expnms, 1])
  covMat = fit$var
  colnames(covMat) <- names(mod$coefficients[, 1])
  seb <- se_comb(expnms, covmat = covMat)
  tstat <- estb/seb
  df <- mod$df.null - length(expnms)
  pval <- 2 - 2 * pt(abs(tstat), df = df)
  pvalz <- 2 - 2 * pnorm(abs(tstat))
  ci <- c(estb + seb * qnorm(alpha/2), estb + seb * qnorm(1 - 
                                                            alpha/2))
  wcoef <- fit$coefficients[expnms]
  names(wcoef) <- gsub("_q", "", names(wcoef))
  poscoef <- which(wcoef > 0)
  negcoef <- which(wcoef <= 0)
  pweights <- abs(wcoef[poscoef])/sum(abs(wcoef[poscoef]))
  nweights <- abs(wcoef[negcoef])/sum(abs(wcoef[negcoef]))
  pos.psi <- sum(wcoef[poscoef])
  neg.psi <- sum(wcoef[negcoef])
  nmpos = names(pweights)
  nmneg = names(nweights)
  se.pos.psi <- se_comb(nmpos, covmat = covMat)
  se.neg.psi <- se_comb(nmneg, covmat = covMat)
  qx <- qdata[, expnms]
  names(qx) <- paste0(names(qx), "_q")
  res <- list(qx = qx, fit = fit, psi = estb, var.psi = seb^2, 
              ci = ci, expnms = expnms, q = q, breaks = br, degree = 1, 
              pos.psi = pos.psi, neg.psi = neg.psi, pweights = sort(pweights, 
                   decreasing = TRUE), nweights = sort(nweights, decreasing = TRUE), 
              psize = sum(abs(wcoef[poscoef])), nsize = sum(abs(wcoef[negcoef])), 
              bootstrap = FALSE, zstat = tstat, pval = pvalz)
  attr(res, "class") <- "qgcompfit"
  res
}

coxmsm.fit <- function(
  f, qdata, intvals, expnms, rr=TRUE, main=TRUE, degree=1, id=NULL, 
  ...){
  #' @title fitting marginal structural model (MSM) based on g-computation with
  #' quantized exposures
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
  #' @param f an r formula representing the conditional model for the outcome, given all
  #' exposures and covariates. Interaction terms that include exposure variables
  #' should be represented via the \code{\link[base]{I}} function
  #' @param qdata a data frame with quantized exposures
  #' @param intvals sequence, the sequence of integer values that the joint exposure 
  #' is 'set' to for estimating the msm. For quantile g-computation, this is just 
  #' 0:(q-1), where q is the number of quantiles of exposure.
  #' @param expnms a character vector with the names of the columns in qdata that represent
  #' the exposures of interest (main terms only!)
  #' @param rr logical, estimate log(risk ratio) (family='binomial' only)
  #' @param main logical, internal use: produce estimates of exposure effect (psi)
  #'  and expected outcomes under g-computation and the MSM
  #' @param degree polynomial basis function for marginal model (e.g. degree = 2
  #'  allows that the relationship between the whole exposure mixture and the outcome
  #'  is quadratic. Default=1)
  #' @param id (optional) NULL, or variable name indexing individual units of 
  #' observation (only needed if analyzing data with multiple observations per 
  #' id/cluster)
  #' @param ... arguments to coxph (e.g. family)
  #' @seealso \code{\link[qgcomp]{qgcomp.cox.boot}}, and \code{\link[qgcomp]{qgcomp.cox.noboot}}
  #' @keywords variance, mixtures
  #' @import survival
  #' @examples
  #' runif(1)
    # not yet implemented
  {
    id = "id__"
    qdata$id__ = 1:dim(qdata)[1]
  }
  # conditional outcome regression fit
  #fit <- glm(f, data = qdata[,!(names(qdata) %in% id)], ...)
  fit <- coxph(f, data = qdata[,!(names(qdata) %in% id)], ...)
  coxfam = list(family='cox', link='log', linkfun=log)
  class(coxfam) = "family"
  fit[['family']] = coxfam # kludge for print function
  bh = basehaz(fit, centered=TRUE)
  bh$bh = c(0, diff(bh$hazard))
  
  if(!is.null(f[2][[1]][4][[1]])){
    envar = f[2][[1]][2][[1]]
    exvar = f[2][[1]][3][[1]]
    entry = as.numeric(with(qdata, eval(envar)))
  } else{
    exvar = f[2][[1]][2][[1]]
    entry = rep(0, nrow(qdata))
  }
  # 
  bd = merge(bh[, c('time', 'bh')], qdata)
  bd$lasttime = (bd$time==with(bd, eval(exvar)))
  #bd = dplyr::filter(bd, time>=with(bd, eval(envar)))
  bd = bd[bd$time>=with(bd, eval(envar)),]
  bd$lhr = predict(fit, newdata = bd, type = 'lp')
  bd$hazard = with(bd, exp(log(bh) + lhr))
  
  # large sample with replacement (mainly for sampling entry times)
  M = 1000
  dd = data.frame(bigid = 1:M, id = sample(bd[,id], size = M, replace=TRUE))
  # create multiple copies with exposures
  
  ########
  ## get predictions (set exposure to 0,1,...,q-1)
  #if(is.null(intvals)){
  #  intvals = (1:length(table(qdata[expnms[1]]))) - 1
  #}
  #predit <- function(idx){
  #  newdata <- qdata
  #  newdata[,expnms] <- idx
  #  suppressWarnings(predict(fit, newdata=newdata, type='response'))
  #}
  #predmat = lapply(intvals, predit)
  ########
  bigd = merge(dd, bd, by.x="id", by.y=id)

  # simulate from conditional fit
  bigd$event = rbinom(nrow(bigd), 1, bigd$hazard)
  bigd = bigd[bigd$event==1 | bigd$lasttime,]
  dd$endtime = with(bigd, tapply(eval(exvar), bigid, min))
  dd$begtime = with(bigd, tapply(eval(envar), bigid, min))
  dd$event   = with(bigd, tapply(event, bigid, max))
  
  # fit using simulated
  
  
  #if(fit$family$family=="gaussian") rr=FALSE
  ### 
  # get predictions (set exposure to 0,1,...,q-1)
  if(is.null(intvals)){
    intvals = (1:length(table(qdata[expnms[1]]))) - 1
  }
  predit <- function(idx){
    newdata <- qdata
    newdata[,expnms] <- idx
    suppressWarnings(predict(fit, newdata=newdata, type='response'))
  }
  predmat = lapply(intvals, predit)
  # fit MSM using g-computation estimates of expected outcomes under joint 
  #  intervention
  nobs <- dim(qdata)[1]
  msmdat <- data.frame(
    Ya = unlist(predmat),
    psi = rep(intvals, each=nobs))
  # to do: allow functional form variations for the MSM via specifying the model formula
  if(!rr) suppressWarnings(msmfit <- glm(Ya ~ poly(psi, degree=degree, raw=TRUE), data=msmdat,...))
  if(rr)  suppressWarnings(msmfit <- glm(Ya ~ poly(psi, degree=degree, raw=TRUE), data=msmdat, family=binomial(link='log'), start=rep(-0.0001, degree+1)))
  res = list(fit=fit, msmfit=msmfit)
  if(main) {
    res$Ya = msmdat$Ya   # expected outcome under joint exposure, by gcomp
    res$Yamsm = predict(msmfit, type='response')
    res$A =  msmdat$psi # joint exposure (0 = all exposures set category with 
    # upper cut-point as first quantile)
  }
  res
}


qgcomp.cox.boot <- function (f, data, expnms = NULL, q = 4, breaks = NULL,
                               B=10, id=NULL, alpha=0.05,
          ...) {
  #' @title estimation of quantile g-computation fit for a survival outcome
  #'  
  #'
  #' @description This function performs quantile g-computation in a survival
  #' setting. 
  #' 
  #' @details For survival outcomes ...
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
  #' @param B Number of bootstrap iterations (default is 10)
  #' @param id (optional) NULL, or variable name indexing individual units of 
  #' observation (only needed if analyzing data with multiple observations per 
  #' id/cluster)
  #' @param alpha alpha level for confidence limit calculation
  #' @param ... arguments to glm (e.g. family)
  #' @seealso \code{\link[qgcomp]{qgcomp.boot}}, and \code{\link[qgcomp]{qgcomp}} 
  #'  for time-fixed outcomes (cross-sectional or cohort design with outcomes measured
  #'  at the end of follow-up) 
  #' @return a qgcompfit object, which contains information about the effect
  #'  measure of interest (psi) and associated variance (var.psi), as well
  #'  as information on the model fit (fit) and information on the 
  #'  weights/standardized coefficients in the positive (pweights) and 
  #'  negative (nweight) directions.
  #' @keywords variance, mixtures
  #' @import survival
#  #' @export
  #' @examples
  #' runif(1)
  # dummy function
}