#qgcomp_long.R: longitudinal quantile g-computation methods
#library(qgcomp)

qgcomp.gee.noboot <- function(f, data, expnms=NULL, q=4, breaks=NULL, id=NULL, 
                            corstr="independence", alpha=0.05, ...){
  #' @title estimation of longitudinal quantile g-computation fit 
  #' (continuous outcome) or conditional quantile odds ratio (binary outcome)
  #'
  #' @description This function estimates a joint effect of multiple exposures
  #'  
  #' 
  #' @details For continuous outcomes, under a linear model with no 
  #' interaction terms, this is equivalent to g-computation of the effect of
  #' increasing every exposure by 1 quantile. For binary outcomes
  #' outcomes, this yields a conditional log odds ratio representing the 
  #' change in the expected conditional odds (conditional on covariates)
  #' from increasing every exposure by 1 quantile. In general, the latter 
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
  #' @param id (required) NULL, or variable name indexing individual units of 
  #' observation (only needed if analyzing data with multiple observations per 
  #' id/cluster)
  #' @param corstr character string with the type of correlation structure implied
  #' by the working covariance matrix: "independence", "exchangeable", "ar1", 
  #' "unstructured" and "userdefined" (See \code{\link[geepack]{geeglm}})
  #' @param alpha alpha level for confidence limit calculation
  #' @param ... arguments to glm (e.g. family)
  #' @seealso \code{\link[qgcomp]{qgcomp.boot}}, and \code{\link[qgcomp]{qgcomp}}
  #' @return a qgcompfit object, which contains information about the effect
  #'  measure of interest (psi) and associated variance (var.psi), as well
  #'  as information on the model fit (fit) and information on the 
  #'  weights/standardized coefficients in the positive (pweights) and 
  #'  negative (nweight) directions.
  #' @keywords variance, mixtures
  #' @import stats geepack
  #' @export
  #' @examples
  #' set.seed(50)
  #' geedat <- data.frame(i = round((1:100)/4), x1=runif(50), x2=runif(50), z=runif(50))
  #' geedat <- within(geedat, {y=runif(50)+i/10})
  #' qgcomp.noboot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), 
  #'    data=geedat, id='i', q=2)
  #' qgcomp.gee.noboot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), 
  #'    data=geedat, id='i', q=2, corstr="independence")
  #' qgcomp.gee.noboot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), 
  #'    data=geedat, id='i', q=2, corstr="exchangeable")
    id___ = NULL
    if (is.null(expnms)) {
      cat("Including all model terms as exposures of interest")
      expnms <- attr(terms(f, data = data), "term.labels")
    }
    if (!is.null(q) | !is.null(breaks)){
      ql <- quantize(data, expnms, q, breaks)
      qdata <- ql$data
      br <- ql$breaks
    } else{
      qdata <- data
      br <- breaks
      }
    if(is.null(id)) {
      # not yet implemented
      id = "id__"
      qdata$id__ = 1:dim(qdata)[1]
    }
    ### 
    # TODO warning about missing data
    ###
    # conditional outcome regression fit
    # workaround for scoping issue in geeglm, but ugly
    qdata$id___ = qdata[,id, drop=TRUE]
    fit <- geeglm(f, data = qdata, 
                  corstr = corstr, id = id___,
                  ...)
    mod <- summary(fit)
    estb <- sum(mod$coefficients[expnms,1, drop=TRUE])
    nm = rownames(mod$coefficients)
    covmat = mod$cov.scaled
    rownames(covmat) <- colnames(covmat) <- nm
    seb <- se_comb(expnms, covmat = covmat)
    tstat <- estb / seb
    df <- fit$df.residual
    pval <- 2 - 2 * pt(abs(tstat), df = df)
    pvalz <- 2 - 2 * pnorm(abs(tstat))
    ci <- c(estb + seb * qnorm(alpha / 2), estb + seb * qnorm(1 - alpha / 2))
    # 'weights'
    wcoef <- fit$coefficients[expnms]
    names(wcoef) <- gsub("_q", "", names(wcoef))
    poscoef <- which(wcoef > 0)
    negcoef <- which(wcoef <= 0)
    pweights <- abs(wcoef[poscoef]) / sum(abs(wcoef[poscoef]))
    nweights <- abs(wcoef[negcoef]) / sum(abs(wcoef[negcoef]))
    # 'post-hoc' positive and negative estimators 
    # similar to constrained gWQS
    pos.psi <- sum(wcoef[poscoef])
    neg.psi <- sum(wcoef[negcoef])
    nmpos <- names(pweights)
    nmneg <- names(nweights)
    se.pos.psi <- se_comb(nmpos, covmat = mod$cov.scaled)
    se.neg.psi <- se_comb(nmneg, covmat = mod$cov.scaled)
    qx <- qdata[, expnms]
    names(qx) <- paste0(names(qx), "_q")
    res <- list(
      qx = qx, fit = fit, psi = estb, var.psi = seb ^ 2, covmat.psi = c('psi1' = seb^2), 
      ci = ci, expnms=expnms, q=q, breaks=br, degree=1,
      pos.psi = pos.psi, neg.psi = neg.psi,
      pweights = sort(pweights, decreasing = TRUE),
      nweights = sort(nweights, decreasing = TRUE), 
      psize = sum(abs(wcoef[poscoef])),
      nsize = sum(abs(wcoef[negcoef])),
      bootstrap=FALSE
    )
      if(fit$family$family=='gaussian'){
        res$tstat <- tstat
        res$df <- df
        res$pval <- pval
      }
      if(fit$family$family=='binomial'){
        res$zstat <- tstat
        res$pval <- pvalz
      }
    attr(res, "class") <- "qgcompfit"
    res
  }


msm.gee.fit <- function(f, qdata, intvals, expnms, rr=TRUE, main=TRUE, degree=1, 
                        corstr="independence", id=NULL, ...){
  #' @title fitting marginal structural model (MSM) based on g-computation with
  #' quantized exposures
  #' @description this is an internal function called by \code{\link[qgcomp]{qgcomp}},
  #'  \code{\link[qgcomp]{qgcomp.boot}}, and \code{\link[qgcomp]{qgcomp.noboot}},
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
  #' @param expnms a character vector with the names of the columns in qdata that represent
  #' the exposures of interest (main terms only!)
  #' @param intvals sequence, the sequence of integer values that the joint exposure 
  #' is 'set' to for estimating the msm. For quantile g-computation, this is just 
  #' 0:(q-1), where q is the number of quantiles of exposure.
  #' @param rr logical, estimate log(risk ratio) (family='binomial' only)
  #' @param main logical, internal use: produce estimates of exposure effect (psi)
  #'  and expected outcomes under g-computation and the MSM
  #' @param degree polynomial basis function for marginal model (e.g. degree = 2
  #'  allows that the relationship between the whole exposure mixture and the outcome
  #'  is quadratic. Default=1)
  #' @param corstr character string with the type of correlation structure implied
  #' by the working covariance matrix: "independence", "exchangeable", "ar1", 
  #' "unstructured" and "userdefined" (See \code{\link[geepack]{geeglm}})
  #' @param id (optional) NULL, or variable name indexing individual units of 
  #' observation (only needed if analyzing data with multiple observations per 
  #' id/cluster)
  #' @param ... arguments to glm (e.g. family)
  #' @seealso \code{\link[qgcomp]{qgcomp.boot}}, and \code{\link[qgcomp]{qgcomp}}
  #' @keywords variance, mixtures
  #' @import stats geepack
  #' @examples
  #' set.seed(50)
  #' geedat <- data.frame(i = round((1:100)/4), x1=runif(50), x2=runif(50), z=runif(50))
  #' geedat <- within(geedat, {y=runif(50)+i/10})
  #' X <- c('x1', 'x2')
  #' qdat <- quantize(geedat, X, q=4)$data
  #' # Note that gee under indepence cov. matrix yields a glm
  #' mod <- qgcomp:::msm.gee.fit(f=y ~ z + x1 + x2 + I(x1*x2), 
  #'         expnms = c('x1', 'x2'), qdata=qdat, intvals=1:4,
  #'         corstr="independence", degree=2, id='i')
  #' mod2 <- qgcomp:::msm.fit(f=y ~ z + x1 + x2 + I(x1*x2), 
  #'         expnms = c('x1', 'x2'), qdata=qdat, intvals=1:4,
  #'         degree=2)
  #' mod3 <- qgcomp:::msm.gee.fit(f=y ~ z + x1 + x2 + I(x1*x2), 
  #'         expnms = c('x1', 'x2'), qdata=qdat, intvals=1:4,
  #'         corstr="exchangeable", degree=2)
  #' summary(mod$fit) # outcome regression model
  #' summary(mod2$fit) # outcome regression model
  #' summary(mod$msmfit) # msm fit (variance not valid - must be obtained via bootstrap)
  #' summary(mod2$msmfit) # msm fit (variance not valid - must be obtained via bootstrap)
  #' summary(mod3$msmfit) # msm fit (variance not valid - must be obtained via bootstrap)
    id___ = NULL
    if(is.null(id)) {
      id <- "id__"
      qdata$id__ <- 1:dim(qdata)[1]
    }
    # conditional outcome regression fit
    # workaround for scoping issue in geeglm, but ugly
    qdata$id___ = qdata[,id, drop=TRUE]
    fit <- geeglm(f, data = qdata, 
                  corstr = corstr, id = id___,
                  ...)
    #fit <- gee(f, data = qdata, 
    #              corstr = corstr, id = id,
    #              ...)
    if(fit$family$family=="gaussian") rr=FALSE
    ### 
    # get predictions (set exposure to 0,1,...,q-1)
    if(is.null(intvals)){
      intvals <- (1:length(table(qdata[expnms[1]]))) - 1
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
      psi = rep(intvals, each=nobs),
      id = rep(qdata[,id, drop=TRUE], times=length(intvals))
      )
    # to do: allow functional form variations for the MSM via specifying the model formula
    if(!rr){
      suppressWarnings(
       msmfit <-  geeglm(Ya ~ poly(psi, degree=degree, raw=TRUE), data=msmdat, 
                  corstr = corstr, id = id, ...)
        )
    }
    if(rr){
      suppressWarnings(
        glmfit <- glm(Ya ~ poly(psi, degree=degree, raw=TRUE), data=msmdat, 
                  family=binomial(link='log'), start=rep(-0.0001, degree+1))
        )
        start <- as.numeric(glmfit$coefficients)
      suppressWarnings(        
        msmfit <-  geese(Ya ~ poly(psi, degree=degree, raw=TRUE), data=msmdat, 
                  corstr = corstr, id = id, 
                  family=binomial(link='log'), b=start)
        )
    }
    res <- list(fit=fit, msmfit=msmfit)
    if(main) {
      res$Ya <- msmdat$Ya   # expected outcome under joint exposure, by gcomp
      res$Yamsm <- predict(msmfit, type='response')
      res$A <- msmdat$psi # joint exposure (0 = all exposures set category with 
       # upper cut-point as first quantile)
    }
    res
}


qgcomp.gee.boot <- function(f, data, expnms=NULL, q=4, breaks=NULL, 
                            corstr = "independence",
                            id=NULL, alpha=0.05, B=200, 
                        rr=TRUE, degree=1, seed=NULL, ...){
  #' @title estimation of quantile g-computation fit, using bootstrap confidence intervals
  #'  
  #' @description This function yields population average effect estimates for 
  #'   both continuous and binary outcomes
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
  #' @param corstr character string with the type of correlation structure implied
  #' by the working covariance matrix: "independence", "exchangeable", "ar1", 
  #' "unstructured" and "userdefined" (See \code{\link[geepack]{geeglm}})
  #' @param id (optional) NULL, or variable name indexing individual units of 
  #' observation (only needed if analyzing data with multiple observations per 
  #' id/cluster)
  #' @param alpha alpha level for confidence limit calculation
  #' @param B integer: number of bootstrap iterations (this should typically be
  #' >=200, though it is set lower in examples to improve run-time).
  #' @param rr logical: if using binary outcome and rr=TRUE, qgcomp.boot will 
  #'   estimate risk ratio rather than odds ratio
  #' @param degree polynomial basis function for marginal model (e.g. degree = 2
  #'  allows that the relationship between the whole exposure mixture and the outcome
  #'  is quadratic.
  #' @param seed integer or NULL: random number seed for replicable bootstrap results
  #' @param ... arguments to glm (e.g. family)
  #' @seealso \code{\link[qgcomp]{qgcomp.noboot}}, and \code{\link[qgcomp]{qgcomp}}
  #' @return a qgcompfit object, which contains information about the effect
  #'  measure of interest (psi) and associated variance (var.psi), as well
  #'  as information on the model fit (fit) and information on the 
  #'  marginal structural model (msmfit) used to estimate the final effect
  #'  estimates.
  #' @keywords variance, mixtures
  #' @import stats geepack
  #' @export
  #' @examples
  #' set.seed(30)
  #' geedat <- data.frame(i = round((1:50)/4), x1=runif(50), x2=runif(50), z=runif(50))
  #' geedat <- within(geedat, {y=rbinom(50, 1, 1/(1+exp(-log((i+1)/10))))})
  #' qgcomp.boot(y ~ z + x1 + x2, family="binomial", expnms = c('x1', 'x2'), 
  #'   data=geedat, q=2, B=10, id='i')
  #' mod = qgcomp.gee.boot(f = y ~ z + x1 + x2, family="binomial", expnms = c('x1', 'x2'), 
  #'   data=geedat, q=2, B=100, id='i', corstr="independence")
  #' mod2 = qgcomp.gee.boot(y ~ z + x1 + x2 + I(x1*x2), family="binomial", expnms = c('x1', 'x2'), 
  #'   data=geedat, q=2, B=100, id='i', corstr="independence", rr=FALSE)
  #' # compared with standard qgcomp
  #' set.seed(30)
  #' geedat <- data.frame(i = round((1:100)/4), x1=runif(100), x2=runif(100), z=runif(100))
  #' geedat <- within(geedat, {y=rbinom(100, 1, 1/(1+exp(-log((i+1)/10))))})
  #' mod3 <- qgcomp.gee.boot(y ~ z + x1 + x2+ I(x1*x2), family=binomial(), expnms = c('x1', 'x2'), 
  #'   data=geedat, q=4, B=10, rr=FALSE, degree=2, id='i', corstr="exchangeable")
  #' plot(mod3)
  #' mod4 <- qgcomp.boot(y ~ z + x1 + x2+ I(x1*x2), family=binomial(), expnms = c('x1', 'x2'), 
  #'   data=geedat, q=4, B=10, rr=FALSE, degree=2, id='i')
  #' plot(mod4)

      # character names of exposure mixture components
    if(is.null(seed)) seed = round(runif(1, min=0, max=1e8))
    if (is.null(expnms)) {
      cat("Including all model terms as exposures of interest")
      expnms <- attr(terms(f, data = data), "term.labels")
    }
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
      intvals <- (1:nvals)-1
    } else {
      # if( is.null(breaks) & is.null(q)) # also includes NA
      qdata <- data
      # if no transformation is made (no quantiles, no breaks given)
      # then draw distribution values from quantiles of all the exposures
      # pooled together
      # TODO: allow user specification of this
      cat("\nNote: using quantiles of all exposures combined in order to set 
          proposed intervention values for overall effect (25th, 50th, 75th %ile)")
      intvals = as.numeric(quantile(unlist(data[,expnms]), c(.25, .5, .75)))
      br <- NULL
    }
    if(is.null(id)) {
      id <- "id__"
      qdata$id__ <- 1:dim(qdata)[1]
    }
    ###
    msmfit <- msm.gee.fit(f, qdata, intvals, expnms, rr, main=FALSE, degree=degree, corstr=corstr, id=id,
                          ...)
    # main estimate  
    estb <- as.numeric(msmfit$msmfit$coefficients[-1])
    #bootstrap to get std. error
    nobs <- dim(qdata)[1]
    nids <- length(unique(qdata[,id, drop=TRUE]))
    psi.only <- function(i=1, f=f, qdata=qdata, intvals=intvals, expnms=expnms, rr=rr, degree=degree, nids=nids, id=id, ...){
      bootids <- data.frame(temp=sort(sample(unique(qdata[,id, drop=TRUE]), nids, replace = TRUE)))
      names(bootids) <- id
      qdata_ <- merge(qdata,bootids, by=id, all.x=FALSE, all.y=TRUE)
      as.numeric(
        msm.gee.fit(f, qdata_, intvals, expnms, rr, main=FALSE, degree=degree, corstr=corstr, id=id,
                ...)$msmfit$coefficients[-1]
      )
    }
    set.seed(seed)
    bootsamps <- sapply(X=1:B, FUN=psi.only,f=f, qdata=qdata, intvals=intvals, 
                       expnms=expnms, rr=rr, degree=degree, nids=nids, id=id, ...)
    if(is.null(dim(bootsamps))) {
      seb <- sd(bootsamps)
      covmat <- var(bootsamps)
      names(covmat) <- 'psi1'
    }else{
      seb <- apply(bootsamps, 1, sd)
      covmat <- cov(t(bootsamps))
      colnames(covmat) <- rownames(covmat) <- paste0("psi", 1:nrow(bootsamps))
    }
    tstat <- estb / seb
    df <- nobs - length(attr(terms(f, data = data), "term.labels")) - 1 - degree # df based on obs - gcomp terms - msm terms
    pval <- 2 - 2 * pt(abs(tstat), df = df)
    pvalz <- 2 - 2 * pnorm(abs(tstat))
    ci <- cbind(estb + seb * qnorm(alpha / 2), estb + seb * qnorm(1 - alpha / 2))
    # 'weights' not applicable in this setting, generally (i.e. if using this function for non-linearity, 
    #   then weights will vary with level of exposure)
    qx <- qdata[, expnms]
    res <- list(
      qx = qx, fit = msmfit$fit, msmfit = msmfit$msmfit, psi = estb, 
      var.psi = seb ^ 2, covmat.psi=covmat, ci = ci,
      expnms=expnms, q=q, breaks=br, degree=degree,
      pos.psi = NULL, neg.psi = NULL, 
      pweights = NULL,nweights = NULL, psize = NULL,nsize = NULL, bootstrap=TRUE,
      y.expected=msmfit$Ya, y.expectedmsm=msmfit$Yamsm, index=msmfit$A,
      bootsamps = bootsamps
    )
      if(msmfit$fit$family$family=='gaussian'){
        res$tstat <- tstat
        res$df <- df
        res$pval <- pval
      }
      if(msmfit$fit$family$family=='binomial'){
        res$zstat <- tstat
        res$pval <- pvalz
      }
    attr(res, "class") <- "qgcompfit"
    res
}
