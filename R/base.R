
se_comb <- function(expnms, covmat, grad=NULL){
  #calculate standard error of weighted linear combination of random variables
  # given a vector of weights and a covariance matrix (not exported)
  # usage: qgcomp:::se_comb(expnms='x', covmat=summary(lmfit)$cov.scaled)
  #calculate standard error of weighted linear combination of random variables
  #This is simple version of the delta method for a linear combination:
  #  f(x) = x1 + x2 + x3
  # given gradient vector G = 
  #   [d(f(x))/dx1 = 1,
  #   d(f(x))/dx2 = 1,
  #   d(f(x))/dx3 = 1]
  # t(G) %*% cov(x) %*% G = delta variance
  #' vcov = rbind(c(1.2, .9),c(.9, 2.0))
  #' colnames(vcov) <- rownames(vcov) <- expnms <- c("x1", "x2")
  #' qgcomp:::se_comb(expnms, vcov, c(1, 0))^2 # returns the given variance
  #' qgcomp:::se_comb(expnms, vcov, c(1, 1)) # default linear MSM fit: all exposures
  #' # have equal weight
  #' qgcomp:::se_comb(expnms, vcov, c(.3, .1)) # used when one exposure contributes
  #'   # to the overall fit more than others  = d(msmeffect)/dx

  if(!is.matrix(covmat)) {
    nm <- names(covmat)
    covmat = matrix(covmat)
    colnames(covmat) <- nm
  }
  weightvec <- rep(0, dim(covmat)[1])
  # eventual extension: allow non-zero 'weights' such that the intervention
  # could correspond to 1 unit increases in some variables, and < 1 unit increases
  # in others
  if(is.null(grad)) weightvec[which(colnames(as.matrix(covmat)) %in% expnms)] <- 1
  if(!is.null(grad)) weightvec[which(colnames(as.matrix(covmat)) %in% expnms)] <- grad
  var <- weightvec %*% covmat %*% weightvec # delta method
  #var <- sum(wcovmat)
  sqrt(var)[1,,drop=TRUE] # should be a scalar
}

grad.poly <- function(intvals, degree){
  # returns matrix with each column referring
  if(degree==1){
    mat <- matrix(1, nrow=length(intvals), ncol=1)
  }else{
    mat <- matrix(1, nrow=length(intvals), ncol=degree)
    for(d in 2:degree){
      mat[,d] <- d*poly(intvals, degree = d-1, simple = TRUE, raw = TRUE)[,d-1]
    }
  }
  mat
}


quantize <- function (data, expnms, q=4, breaks=NULL) {
  #' @title create variables representing indicator functions with cutpoints defined
  #' by quantiles
  #' @description This function creates categorical variables in place of the
  #' exposure variables named in 'expnms.' For example, a continuous exposure
  #' 'x1' will be replaced in the output data by another 'x1' that takes on values
  #' 0:(q-1), where, for example, the value 1 indicates that the original x1 value
  #' falls between the first and the second quantile.
  #' @details This function is a vectorized version of `quantile_f` from the `gWQS` 
  #' package that also allows the use of externally defined breaks
  #' @param data a data frame
  #' @param expnms a character vector with the names of  the columns to be
  #' quantized
  #' @param q integer, number of quantiles used in creating quantized variables
  #' @param breaks (optional) list of (equal length) numeric vectors that 
  #' characterize the minimum value of each category for which to 
  #' break up the variables named in expnms. This is an alternative to using 'q'
  #' to define cutpoints.
  #' @keywords variance, mixtures
  #' @import stats
  #' @export
  #' @examples
  #' set.seed(1232)
  #' dat = data.frame(y=runif(100), x1=runif(100), x2=runif(100), z=runif(100))
  #' qdata = quantize(data=dat, expnms=c("x1", "x2"), q=4)
  #' table(qdata$data$x1)
  #' table(qdata$data$x2)
  #' summary(dat[c("y", "z")]);summary(qdata$data[c("y", "z")]) # not touched
  #' dat = data.frame(y=runif(100), x1=runif(100), x2=runif(100), z=runif(100))
  #' # using 'breaks' requires specifying min and max (the qth quantile)
  #' # example with theoretical quartiles (could be other relevant values)
  #' qdata2 = quantize(data=dat, expnms=c("x1", "x2"), 
  #'    breaks=list(c(-1e64, .25, .5, .75, 1e64), 
  #'                c(-1e64, .25, .5, .75, 1e64)
  #'                ))
  #' table(qdata2$data$x1)
  #' table(qdata2$data$x2)
    e <- new.env()
    e$retbr <- list()
    qt <- function(i){
      # not exported
        datmat <- as.numeric(unlist(data[, expnms[i]]))
        if(!is.null(breaks)){
          # prioritize breaks if given by user
          br  <- breaks[[i]]
          e$retbr[[i]] <- breaks[[i]]
        }else{
          br <- unique(quantile(datmat, probs = seq(0, 1, by = 1 / q), na.rm = TRUE))
          br[1] <- -1e64
          br[length(br)] <- 1e64
          e$retbr[[i]] <- br 
        }
        cut(datmat, breaks = br, labels = FALSE,
             include.lowest = TRUE) - 1
    }
    data[, expnms] <- sapply(1:length(expnms), qt)
    return(list(data=data, breaks=e$retbr))
}


msm.fit <- function(f, qdata, intvals, expnms, rr=TRUE, main=TRUE, degree=1, id=NULL, ...){
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
  #' @param id (optional) NULL, or variable name indexing individual units of 
  #' observation (only needed if analyzing data with multiple observations per 
  #' id/cluster)
  #' @param ... arguments to glm (e.g. family)
  #' @seealso \code{\link[qgcomp]{qgcomp.boot}}, and \code{\link[qgcomp]{qgcomp}}
  #' @keywords variance, mixtures
  #' @import stats
  #' @examples
  #' set.seed(50)
  #' dat <- data.frame(y=runif(200), x1=runif(200), x2=runif(200), z=runif(200))
  #' X <- c('x1', 'x2')
  #' qdat <- quantize(dat, X, q=4)$data
  #' mod <- qgcomp:::msm.fit(f=y ~ z + x1 + x2 + I(x1*x2), 
  #'         expnms = c('x1', 'x2'), qdata=qdat, intvals=1:4)
  #' summary(mod$fit) # outcome regression model
  #' summary(mod$msmfit) # msm fit (variance not valid - must be obtained via bootstrap)
    if(is.null(id)) {
      id <- "id__"
      qdata$id__ <- 1:dim(qdata)[1]
    }
    # conditional outcome regression fit
    fit <- glm(f, data = qdata[,!(names(qdata) %in% id)], ...)
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
      psi = rep(intvals, each=nobs))
    # to do: allow functional form variations for the MSM via specifying the model formula
    if(!rr) suppressWarnings(msmfit <- glm(Ya ~ poly(psi, degree=degree, raw=TRUE), data=msmdat,...))
    if(rr)  suppressWarnings(msmfit <- glm(Ya ~ poly(psi, degree=degree, raw=TRUE), data=msmdat, family=binomial(link='log'), start=rep(-0.0001, degree+1)))
    res <- list(fit=fit, msmfit=msmfit)
    if(main) {
      res$Ya <- msmdat$Ya   # expected outcome under joint exposure, by gcomp
      res$Yamsm <- predict(msmfit, type='response')
      res$A <- msmdat$psi # joint exposure (0 = all exposures set category with 
       # upper cut-point as first quantile)
    }
    res
}


qgcomp.noboot <- function(f, data, expnms=NULL, q=4, breaks=NULL, id=NULL, alpha=0.05, ...){
  #' @title estimation of quantile g-computation fit (continuous outcome)
  #'  or conditional quantile odds ratio (binary outcome)
  #'
  #' @description This function mimics the output of a weighted quantile sums regression in 
  #' large samples. 
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
  #' @import stats
  #' @export
  #' @examples
  #' set.seed(50)
  #' dat <- data.frame(y=runif(50), x1=runif(50), x2=runif(50), z=runif(50))
  #' qgcomp.noboot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2)
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
    fit <- glm(f, data = qdata[,!(names(qdata) %in% id), drop=FALSE], ...)
    mod <- summary(fit)
    estb <- sum(mod$coefficients[expnms,1, drop=TRUE])
    seb <- se_comb(expnms, covmat = mod$cov.scaled)
    tstat <- estb / seb
    df <- mod$df.null - length(expnms)
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

#TODO: explain (log) better - here and in the noboot
qgcomp.boot <- function(f, data, expnms=NULL, q=4, breaks=NULL, id=NULL, alpha=0.05, B=200, 
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
  #' @import stats
  #' @export
  #' @examples
  #' set.seed(30)
  #' # continuous outcome
  #' dat <- data.frame(y=rnorm(100), x1=runif(100), x2=runif(100), z=runif(100))
  #' # Conditional linear slope
  #' qgcomp.noboot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=4, family=gaussian())
  #' # Marginal linear slope (population average slope, for a purely linear, 
  #' #  additive model this will equal the conditional)
  #' qgcomp.boot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=4, 
  #'   family=gaussian(), B=10) #increase B to at least 200 in actual examples
  #'   
  #' # Population average mixture slope which accounts for non-linearity and interactions
  #' qgcomp.boot(y ~ z + x1 + x2 + I(x1^2) + I(x2*x1), family="gaussian", 
  #'  expnms = c('x1', 'x2'), data=dat, q=4, B=10)
  #'  
  #' # binary outcome
  #' dat <- data.frame(y=rbinom(50,1,0.5), x1=runif(50), x2=runif(50), z=runif(50))
  #' 
  #' # Conditional mixture OR
  #' qgcomp.noboot(y ~ z + x1 + x2, family="binomial", expnms = c('x1', 'x2'), 
  #'   data=dat, q=2)
  #'   
  #' #Marginal mixture OR (population average OR - in general, this will not equal the 
  #' # conditional mixture OR due to non-collapsibility of the OR)
  #' qgcomp.boot(y ~ z + x1 + x2, family="binomial", expnms = c('x1', 'x2'), 
  #'   data=dat, q=2, B=10)
  #'   
  #' # Population average mixture RR
  #' qgcomp.boot(y ~ z + x1 + x2, family="binomial", expnms = c('x1', 'x2'), 
  #'   data=dat, q=2, rr=TRUE, B=10)
  #'   
  #' # Population average mixture RR, indicator variable representation of x2
  #' # note that I(x==...) operates on the quantile-based category of x,
  #' # rather than the raw value
  #' res = qgcomp.boot(y ~ z + x1 + I(x2==1) + I(x2==2) + I(x2==3), 
  #'   family="binomial", expnms = c('x1', 'x2'), data=dat, q=4, rr=TRUE, B=10)
  #' res$fit  
  #' plot(res)
  #' 
  #' # now add in a non-linear MSM
  #' res2 = qgcomp.boot(y ~ z + x1 + I(x2==1) + I(x2==2) + I(x2==3), 
  #'   family="binomial", expnms = c('x1', 'x2'), data=dat, q=4, rr=TRUE, B=10, 
  #'   degree=2)
  #' res2$fit  
  #' res2$msmfit  
  #' plot(res2)
  #' # Log risk ratio per one IQR change in all exposures (not on quantile basis)
  #' dat$x1iqr <- dat$x1/with(dat, diff(quantile(x1, c(.25, .75))))
  #' dat$x2iqr <- dat$x2/with(dat, diff(quantile(x2, c(.25, .75))))
  #' # note that I(x>...) nowoperates on the untransformed value of x,
  #' # rather than the raw value
  #' res2 = qgcomp.boot(y ~ z + x1iqr + I(x2iqr>0.1) + I(x2>0.4) + I(x2>0.9), 
  #'   family="binomial", expnms = c('x1iqr', 'x2iqr'), data=dat, q=NULL, rr=TRUE, B=10, 
  #'   degree=2)
  #' res2
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
    msmfit <- msm.fit(f, qdata, intvals, expnms, rr, main=TRUE,degree=degree, id=id, ...)
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
        msm.fit(f, qdata_, intvals, expnms, rr, main=FALSE, degree, id,
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


qgcomp <- function(f,data=data,family=gaussian(),rr=TRUE,...){
  #' @title estimation of quantile g-computation fit
  #' 
  #' @description This function automatically selects between qgcomp.noboot and qgcomp.boot
  #'  to select the most efficient approach to estimate the average expected 
  #'  change in the (log) outcome per quantile increase in the joint 
  #'  exposure to all exposures in `expnms'
  #'
  #' @param f R style formula
  #' @param data data frame
  #' @param family `gaussian()` or `binomial()`
  #' @param rr logical: if using binary outcome and rr=TRUE, qgcomp.boot will 
  #' estimate risk ratio rather than odds ratio. Note, to get population average 
  #' effect estimates for a binary outcome, set rr=TRUE (default: ORs are generally not
  #' of interest as population average effects, so if rr=FALSE then a conditional
  #' OR will be estimated, which cannot be interpreted as a population average
  #' effect
  #' @param ... arguments to qgcomp.noboot or qgcomp.boot (e.g. q)
  #' @seealso \code{\link[qgcomp]{qgcomp.noboot}} and \code{\link[qgcomp]{qgcomp.boot}}
  #' @return a qgcompfit object, which contains information about the effect
  #'  measure of interest (psi) and associated variance (var.psi), as well
  #'  as information on the model fit (fit) and possibly information on the 
  #'  marginal structural model (msmfit) used to estimate the final effect
  #'  estimates (qgcomp.boot only). If appropriate, weights are also reported.
  #' @keywords variance, mixtures
  #' @import stats
  #' @export
  #' @examples
  #' set.seed(50)
  #' dat <- data.frame(y=runif(50), x1=runif(50), x2=runif(50), z=runif(50))
  #' qgcomp.noboot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2)
  #' qgcomp.boot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, B=10, seed=125)
  #' # automatically selects appropriate method
  #' qgcomp(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2)
  #' # note for binary outcome this will 
  #' dat <- data.frame(y=rbinom(50, 1, 0.5), x1=runif(50), x2=runif(50), z=runif(50))
  #' qgcomp.noboot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, family=binomial())
  #' qgcomp.boot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, B=10, seed=125, 
  #'   family=binomial())
  #' # automatically selects appropriate method
  #' qgcomp(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, family=binomial())
  #' qgcomp(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, family=binomial(), rr=TRUE)
  terms <- attr(terms(f,data=data), 'term.labels')
  doboot <- ifelse(isTRUE(grep("I\\(", terms)>0), TRUE, FALSE)
  if(rr | doboot){
    res <- qgcomp.boot(f=f,data=data,family=family,rr=rr,...)
  }else{
    res <- qgcomp.noboot(f=f,data=data,family=family,...)
  }
  res
}

print.qgcompfit <- function(x, ...){
  #' @title default printing method for a qgcompfit object
  #' 
  #' @description Gives variable output depending on whether `qgcomp.noboot` or `qgcomp.boot`
  #' is called. For `qgcomp.noboot` will output final estimate of joint exposure
  #' effect (similar to the 'index' effect in weighted quantile sums), as well
  #' as estimates of the 'weights' (standardized coefficients). For `qgcomp.boot`,
  #' the marginal effect is given, but no weights are reported since this approach
  #' generally incorporates non-linear models with interaction terms among exposures,
  #' which preclude weights with any useful interpretation.
  #' 
  #' @param x "qgcompfit" object from `qgcomp`, `qgcomp.noboot` or `qgcomp.boot` 
  #' function
  #' @param ... unused
  #' @seealso \code{\link[qgcomp]{qgcomp.noboot}}, \code{\link[qgcomp]{qgcomp.boot}}, and \code{\link[qgcomp]{qgcomp}}
  #' @keywords variance, mixtures
  #' @export
  #' @examples
  #' set.seed(50)
  #' dat <- data.frame(y=runif(50), x1=runif(50), x2=runif(50), z=runif(50))
  #' obj1 <- qgcomp.noboot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2)
  #' obj2 <- qgcomp.boot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, B=10, seed=125)
  #' # does not need to be explicitly called, but included here for clarity
  #' print(obj1)
  #' print(obj2)
  fam <- x$fit$family$family
  if(!is.null(x$psize)) {
    cat(paste0("Scaled effect size (positive direction, sum of positive coefficients = ", signif(x$psize, 3) , ")\n"))
    if (length(x$pweights) > 0) {
      print(x$pweights, digits = 3)
    } else cat("None\n")
    cat("\n")
  }
  if(!is.null(x$nsize)) {
    cat(paste0("Scaled effect size (negative direction, sum of negative coefficients = ", signif(-x$nsize, 3) , ")\n"))
    if (length(x$nweights) > 0) {
      print(x$nweights, digits = 3)
    } else cat("None\n")
    cat("\n")
  }
  if (fam == "binomial"){
    estimand <- 'OR'
    if(x$bootstrap && x$msmfit$family$link=='log') estimand = 'RR'
    cat(paste0("Mixture log(",estimand,")", ifelse(x$bootstrap, " (bootstrap CI)", " (Delta method CI)"), ":\n\n"))
    if(is.null(dim(x$ci))){
      pdat <- cbind(Estimate=x$psi, "Std. Error"=sqrt(x$var.psi), "Lower CI"=x$ci[1], "Upper CI"=x$ci[2], "Z value"=x$zstat, "Pr(>|z|)"=x$pval)
      rownames(pdat) <- paste0('psi',1:length(x$psi))
      printCoefmat(pdat,has.Pvalue=TRUE,tst.ind=5L,signif.stars=FALSE, cs.ind=1L:2)
    } else{
      pdat <- cbind(Estimate=x$psi, "Std. Error"=sqrt(x$var.psi), "Lower CI"=x$ci[,1], "Upper CI"=x$ci[,2], "Z value"=x$zstat, "Pr(>|z|)"=x$pval)
      rownames(pdat) <- paste0('psi',1:length(x$psi))
      printCoefmat(pdat,has.Pvalue=TRUE,tst.ind=5L,signif.stars=FALSE, cs.ind=1L:2)
    }
  }
  if (fam == "gaussian" | fam == "cox"){
    cat(paste0("Mixture slope parameters", ifelse(x$bootstrap, " (bootstrap CI)", " (Delta method CI)"), ":\n\n"))
    if(is.null(dim(x$ci))){
      pdat <- cbind(Estimate=x$psi, "Std. Error"=sqrt(x$var.psi), "Lower CI"=x$ci[1], "Upper CI"=x$ci[2], "t value"=x$tstat, "Pr(>|t|)"=x$pval)
      rownames(pdat) <- paste0('psi',1:length(x$psi))
      printCoefmat(pdat,has.Pvalue=TRUE,tst.ind=5L,signif.stars=FALSE, cs.ind=1L:2)
    } else{
      pdat <- cbind(Estimate=x$psi, "Std. Error"=sqrt(x$var.psi), "Lower CI"=x$ci[,1], "Upper CI"=x$ci[,2], "t value"=x$tstat, "Pr(>|t|)"=x$pval)
      rownames(pdat) <- paste0('psi',1:length(x$psi))
      printCoefmat(pdat,has.Pvalue=TRUE,tst.ind=5L,signif.stars=FALSE, cs.ind=1L:2)
    }
  }
}


plot.qgcompfit <- function(x, suppressprint=FALSE, ...){
  #' @title plot.qgcompfit: default plotting method for a qgcompfit object
  #'
  #' @description Plot a quantile g-computation object. For qgcomp.noboot, this function will
  #' create a butterfly plot of weights. For qgcomp.boot, this function will create
  #' a box plot with smoothed line overlaying that represents a non-parametric
  #' fit of a model to the expected outcomes in the population at each quantile
  #' of the joint exposures (e.g. '1' represents 'at the first quantile for
  #' every exposure')
  #' 
  #' @param x "qgcompfit" object from `qgcomp.noboot` or  `qgcomp.boot` functions
  #' @param suppressprint If TRUE, suppresses the plot, rather than printing it 
  #'   by default (it can be saved as a ggplot2 object and used programatically)
  #'   (default = FALSE)
  #' @param ... unused
  #' @seealso \code{\link[qgcomp]{qgcomp.noboot}}, \code{\link[qgcomp]{qgcomp.boot}}, and \code{\link[qgcomp]{qgcomp}}
  #' @import ggplot2 grid gridExtra
  #' @export
  #' @examples
  #' set.seed(12)
  #' dat <- data.frame(x1=(x1 <- runif(100)), x2=runif(100), x3=runif(100), z=runif(100),
  #'                   y=runif(100)+x1+x1^2)
  #' ft <- qgcomp.noboot(y ~ z + x1 + x2 + x3, expnms=c('x1','x2','x3'), data=dat, q=4)
  #' ft
  #' plot(ft)
  #' # examinging fit
  #' plot(ft$fit, which=1) # residual vs. fitted is not straight line!
  #' 
  #' # using non-linear outcome model
  #' ft2 <- qgcomp.boot(y ~ z + x1 + x2 + x3 + I(x1*x1), expnms=c('x1','x2','x3'), 
  #' data=dat, q=4, B=10)
  #' ft2
  #' plot(ft2$fit, which=1) # much better looking fit diagnostics suggests
  #' # it is better to include interaction term for x
  #' plot(ft2) # the msm predictions don't match up with a smooth estimate
  #' # of the expected outcome, so we should consider a non-linear MSM

  #' # using non-linear marginal structural model
  #' ft3 <- qgcomp.boot(y ~ z + x1 + x2 + x3 + I(x1*x1), expnms=c('x1','x2','x3'), 
  #' data=dat, q=4, B=10, degree=2)
  #' # plot(ft3$fit, which=1) - not run - this is identical to ft2 fit
  #' plot(ft3) # the MSM estimates look much closer to the smoothed estimates
  #' # suggesting the non-linear MSM fits the data better and should be used
  #' # for inference about the effect of the exposure
  ymin <- ymax <- w <- v <- NULL # appease R CMD check
  theme_butterfly_l <- list(theme(
    legend.position = c(0,0), 
    legend.justification = c(0,0),
    legend.background = element_blank(), 
    panel.background = element_blank(), 
    panel.grid.major.x = element_blank(), 
    panel.grid.major.y = element_blank(), 
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(colour = "black"), 
    #axis.text = element_text(colour="black", face="bold", size=14, family="Helvetica"), 
    #axis.title = element_text(size=16, face="bold", family="Helvetica"), 
    axis.text = element_text(colour="black", face="bold", size=14), 
    axis.title = element_text(size=16, face="bold"), 
    legend.key = element_blank(),
    plot.margin = unit(c(t=1, r=0, b=.75, l=0.5), "cm"),
    panel.border = element_blank()))

  theme_butterfly_r <- list(theme(
    panel.background = element_blank(), 
    panel.grid.major.x = element_blank(), 
    panel.grid.major.y = element_blank(), 
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(colour = "black"), 
    #axis.text.x = element_text(colour="black", face="bold", size=14, family="Helvetica"), 
    #axis.title.x = element_text(size=16, face="bold", family="Helvetica"), 
    axis.text.x = element_text(colour="black", face="bold", size=14), 
    axis.title.x = element_text(size=16, face="bold"), 
    axis.ticks.y = element_blank(), 
    axis.text.y = element_blank(), 
    axis.title.y = element_blank(), 
    legend.key = element_blank(),
    plot.margin = unit(c(t=1, r=0.5, b=.75, l=0.0), "cm"),
    panel.border = element_blank()))

  nms = unique(names(sort(c(x$pweights, x$nweights), decreasing = FALSE)))
  
  #vpl <- grid::viewport(width=0.525, height=1, x=0, y=0, just=c("left", "bottom"))
  #vpr <- grid::viewport(width=0.475, height=1, x=0.525, y=0, just=c("left", "bottom"))
  if(!x$bootstrap){
    if(length(x$pweights)==0) x$pweights = x$nweights*0
    if(length(x$nweights)==0) x$nweights = x$pweights*0
    pright <- ggplot() + 
    stat_identity(aes(x=v, y=w), position = "identity", geom="bar", 
                  data=data.frame(w=x$pweights, v=names(x$pweights))) + 
    scale_y_continuous(name="Positive weights", expand=c(0.000,0.000), breaks=c(0.25, 0.5, 0.75)) +
    scale_x_discrete(limits=nms, breaks=nms, labels=nms, drop=FALSE, position="top") +
    geom_hline(aes(yintercept=0)) + 
    coord_flip(ylim=c(0,1)) + 
    theme_butterfly_r
    pleft <- ggplot() + 
    stat_identity(aes(x=v, y=w), position = "identity", geom="bar", 
                  data=data.frame(w=x$nweights, v=names(x$nweights))) + 
    scale_y_reverse(name="Negative weights", expand=c(0.000,0.000), breaks=c(0.25, 0.5, 0.75)) +
    scale_x_discrete(name="Variable", limits=nms, breaks=nms, labels=nms, drop=FALSE) +
    geom_hline(aes(yintercept=0)) + 
    coord_flip(ylim=c(0,1)) + 
    theme_butterfly_l
    if((length(x$nweights)>0 & length(x$pweights)>0)){
      maxstr = max(mapply(nchar, c(names(x$nweights), names(x$pweights))))
      lw = 1+maxstr/20
      p1 <- gridExtra::arrangeGrob(grobs=list(pleft, pright), ncol=2, padding=0.0, widths=c(lw,1))
      grid::grid.newpage()
      grid::grid.draw(p1)
    }
  }
  if(x$bootstrap){
       # variance based on delta method and knowledge that non-linear
       #functions will always be polynomials in qgcomp


       # default plot for bootstrap results (no weights obtained)
   p <- ggplot() 
     if(x$msmfit$family$family=='gaussian'){
       #confidence band
       y = x$y.expectedmsm
       COV = x$covmat.psi
       intvals = as.numeric(names(table(x$index)))
       grad = grad.poly(intvals, x$degree)
       py = tapply(x$y.expectedmsm, x$index, mean)
       varpy = 0*py
       for(i in 1:length(intvals)){
         varpy[i] = se_comb(expnms = paste0('psi', 1:x$degree), 
                        covmat=COV, grad = grad[i,])
       }
       pyup = py + qnorm(.975)*sqrt(varpy)
       pydo = py + qnorm(.025)*sqrt(varpy)
       p <- p + geom_ribbon(aes(x=x,ymin=ymin,ymax=ymax, 
                                fill="Model confidence band"),
                            data=data.frame(ymin=pydo, ymax=pyup, x=intvals/max(intvals))) +
                    geom_line(aes(x=x,y=y, color="Model fit"),
                            data=data.frame(y=y, x=x$index/max(x$index)))
     }
     if(x$msmfit$family$family=='binomial'){
       y = x$y.expectedmsm # probabilities (not on model scale)
       #variance/gradient on model scales
       COV = x$covmat.psi
       intvals = as.numeric(names(table(x$index)))
       grad = grad.poly(intvals, x$degree)
       py = tapply(x$y.expectedmsm, x$index, mean)
       varpy = 0*py   
       for(i in 1:length(intvals)){
         varpy[i] = se_comb(expnms = paste0('psi', 1:x$degree), 
                        covmat=COV, grad = grad[i,])
       }

       if(x$msmfit$family$link=='log'){
         pyup = exp(log(py) + qnorm(.975)*sqrt(varpy))
         pydo = exp(log(py) + qnorm(.025)*sqrt(varpy))       

       }
       if(x$msmfit$family$link=='logit'){
         pyup = 1/(1+exp(-(log(py/(1-py)) + qnorm(.975)*sqrt(varpy))))
         pydo = 1/(1+exp(-(log(py/(1-py)) + qnorm(.025)*sqrt(varpy))))      
       }

       p <- p + geom_ribbon(aes(x=x,ymin=ymin,ymax=ymax, 
                                fill="Model confidence band"),
                            data=data.frame(ymin=pydo, ymax=pyup, x=intvals/max(intvals))) +
                    geom_line(aes(x=x,y=y, color="Model fit"),
                            data=data.frame(y=y, x=x$index/max(x$index)))
     }
     p <- p + geom_smooth(aes(x=x,y=y, color="Smooth fit"),
                          data=data.frame(y=x$y.expected, x=x$index/max(x$index)), 
                          method = 'gam', 
                          formula=y~s(x, k=4,fx=TRUE), se = FALSE) + 
     scale_x_continuous(name=("Joint exposure quantile")) + 
     scale_y_continuous(name="E(outcome)") + 
     scale_fill_grey(name="", start=.9) + 
     scale_colour_grey(name="", start=0.0, end=0.6) + 
     theme_classic()
     if(!suppressprint) print(p)
  }
  #grid.text("Density", x=0.55, y=0.1, gp=gpar(fontsize=14, fontface="bold", fontfamily="Helvetica"))
}

predict.qgcompfit <- function(object, expnms=NULL, newdata=NULL, type="response", ...){
  #' @title predict.qgcompfit: default prediction method for a qgcompfit object
  #'
  #' @description get predicted values from a qgcompfit object, or make predictions
  #' in a new set of data based on the qgcomfit object. Note that when making predictions
  #' from an object from qgcomp.boot, the predictions are made from the g-computation
  #' model rather than the marginal structural model. Predictions from the marginal
  #' structural model can be obtained via \code{\link[qgcomp]{msm.predict}}
  #' 
  #' @param object "qgcompfit" object from `qgcomp.noboot` or  `qgcomp.boot` functions
  #' @param expnms character vector of exposures of interest
  #' @param newdata (optional) new set of data with all predictors from "qgcompfit" object
  #' @param type  (from predict.glm) the type of prediction required. The default 
  #' is on the scale of the linear predictors; the alternative "response" is on 
  #' the scale of the response  variable. Thus for a default binomial model the 
  #' default predictions are of log-odds (probabilities on logit scale) and 
  #' type = "response" gives the predicted probabilities. The "terms" option 
  #' returns a matrix giving the fitted values of each term in the model formula 
  #' on the linear predictor scale.
  #' @param ... arguments to predict.glm
  #' @export
  #' @examples
  #' set.seed(50)
  #' dat <- data.frame(y=runif(50), x1=runif(50), x2=runif(50), z=runif(50))
  #' obj1 <- qgcomp.noboot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2)
  #' obj2 <- qgcomp.boot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, B=10, seed=125)
  #' set.seed(52)
  #' dat2 <- data.frame(y=runif(50), x1=runif(50), x2=runif(50), z=runif(50))
  #' summary(predict(obj1, expnms = c('x1', 'x2'), newdata=dat2))
  #' summary(predict(obj2, expnms = c('x1', 'x2'), newdata=dat2))
 if(is.null(newdata)){
   pred <- predict(object$fit, type=type, ...) 
  }
 if(!is.null(newdata)){
   newqdata <- quantize(newdata, expnms, q=NULL, object$breaks)$data
   pred <- predict(object$fit, newdata=newqdata, type=type, ...) 
 }
  return(pred)
}

msm.predict <- function(object, newdata=NULL){
  #' @title msm.predict: secondary prediction method for the MSM within a qgcompfit 
  #' object. 
  #' 
  #' @description Makes predictions from the MSM (rather than the g-computation 
  #' model) from a "qgcompfit" object. Generally, this should not be used in 
  #' favor of the default \code{\link[qgcomp]{predict.qgcompfit}} function. This
  #' function can only be used following the `qgcomp.boot` function. For the 
  #' `qgcomp.noboot` function, \code{\link[qgcomp]{predict.qgcompfit}} gives 
  #' identical inference to predicting from an MSM.
  #'
  #' @description get predicted values from a qgcompfit object from
  #' \code{\link[qgcomp]{qgcomp.boot}}
  #' 
  #' @param object "qgcompfit" object from `qgcomp.boot` function
  #' @param newdata (optional) new set of data (data frame) with a variable 
  #' called `psi` representing the joint exposure level of all exposures
  #' under consideration
  #' @export
  #' @examples
  #' set.seed(50)
  #' dat <- data.frame(y=runif(50), x1=runif(50), x2=runif(50), z=runif(50))
  #' obj <- qgcomp.boot(y ~ z + x1 + x2 + I(z*x1), expnms = c('x1', 'x2'), data=dat, q=4, B=10, seed=125)
  #' dat2 <- data.frame(psi=seq(1,4, by=0.1))
  #' summary(msm.predict(obj))
  #' summary(msm.predict(obj, newdata=dat2))
 if(!object$bootstrap) stop("only valid for results from qgcomp.boot function")
 if(is.null(newdata)){
   pred <- predict(object$msmfit, type='response') 
  }
 if(!is.null(newdata)){
   pred <- predict(object$msmfit, newdata=newdata, type='response') 
 }
  return(pred)
}

