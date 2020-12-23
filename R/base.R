
se_comb <- function(expnms, covmat, grad=NULL){
  #' @title Calculate standard error of weighted linear combination of random variables
  #' @description This function uses the Delta method to calculate standard errors of linear
  #' functions of variables (similar to `lincom` in Stata). Generally, users will not need to 
  #' call this function directly.
  #' @details This function takes inputs of a set of exposure names (character vector)
  #' and a covariance matrix (with colnames/rownames that contain the full set
  #' of exposure names), as well as a possible `grad` parameter to calculate
  #' the variance of a weighted combination of the exposures in `expnms`, where the
  #' weights are based off of `grad` (which defaults to 1, so that this function
  #' yields the variance of a sum of all variables in `expnms`)
  #' 
  #' Here is simple version of the delta method for a linear combination of 
  #' three model coefficients:
  #' 
  #'  \eqn{f(\beta) = \beta_1 + \beta_2 + \beta_3}
  #' given gradient vector 
  #' \deqn{G = [d(f(\beta))/d\beta_1 = 1,
  #'   d(f(\beta))/d\beta_2 = 1,
  #'   d(f(\beta))/d\beta_3 = 1]}
  #' \eqn{t(G) Cov(\beta)  G} = delta method variance, where t() is the transpose operator
  #' 
  #' @param expnms a character vector with the names of the columns to be
  #' of interest in the covariance matrix for a which a standard error will be
  #' calculated (e.g. same as expnms in qgcomp fit)
  #' @param covmat covariance matrix for parameters, e.g. from a model or 
  #' bootstrap procedure
  #' @param grad the "weight" vector for calculating the contribution of each variable
  #' in expnms to the final standard error. For a linear combination, this is equal 
  #' to a vector of ones (and is set automatically). Or can be calculated via the 
  #' grad.poly procedure, in the case of coming up with proper weights when the combination
  #' of expnms derives from a polynomial function (as in qgcomp.boot with degree>1).
  #
  #' @examples
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

vc_comb <- function(aname="(Intercept)", expnms, covmat, grad=NULL){
  #' @title Calculate covariance matrix between one random variable and a linear combination of 
  #' random variables
  #' @description This function uses the Delta method to calculate a covariance matrix of linear
  #' functions of variables and is used internally in qgcomp. Generally, users will not need to 
  #' call this function directly.
  #' @details This function takes inputs of a name of random variable (character), as
  #' set of exposure names (character vector) and a covariance matrix (with colnames/rownames 
  #' that contain the indepdendent variable and the full set
  #' of exposure names). See \code{\link[qgcomp]{se_comb}} for details on variances of sums
  #' of random variables. Briefly, for variables A, B and C with covariance matrix Cov(A,B,C),
  #' we can calculate the covariance Cov(A,B+C) with the formula Cov(A,B) + Cov(A,C), and
  #' Cov(A,B+C+D) = Cov(A,B) + Cov(A,C) + Cov(A,D), and so on.
  #' 
  #' @param aname character scalar with the name of the first column of interest (e.g. variable
  #' A in the examples given in the details section)
  #' @param expnms a character vector with the names of the columns to be
  #' of interest in the covariance matrix for a which a standard error will be
  #' calculated (e.g. same as expnms in qgcomp fit)
  #' @param covmat covariance matrix for parameters, e.g. from a model or 
  #' bootstrap procedure
  #' @param grad not yet used
  #' 
  #' @return A covariance matrix
  #
  #' @examples
  #' vcov = rbind(c(0.010051348, -0.0039332248, -0.0036965571), 
  #'              c(-0.003933225,  0.0051807876,  0.0007706792),
  #'              c(-0.003696557,  0.0007706792,  0.0050996587))
  #' colnames(vcov) <- rownames(vcov) <- c("(Intercept)", "x1", "x2")
  #' expnms <- rownames(vcov)[2:3]
  #' aname = rownames(vcov)[1]
  #' qgcomp:::vc_comb(aname, expnms, vcov) # returns the given covariance matrix

  if(!is.matrix(covmat)) {
    nm <- names(covmat)
    covmat = matrix(covmat)
    colnames(covmat) <- nm
  }
  weightvec <- rep(0, dim(covmat)[1])
  # eventual extension: allow non-unity 'weights' such that the intervention
  # could correspond to 1 unit increases in some variables, and < 1 unit increases
  # in others
  #if(!is.null(grad)) weightvec[which(colnames(as.matrix(covmat)) %in% expnms)] <- grad
  if(!is.null(grad)) grad = NULL # not yet used
  if(is.null(grad)) weightvec[which(colnames(as.matrix(covmat)) %in% expnms)] <- 1
  outcov = matrix(NA, nrow=2, ncol=2)
  acol = which(colnames(as.matrix(covmat)) %in% aname)
  bcol = which(colnames(as.matrix(covmat)) %in% expnms)
  outcov[1,1] <- covmat[acol,acol]
  outcov[1,2] <- outcov[2,1] <- sum(covmat[acol, bcol])
  outcov[2,2] <- weightvec %*% covmat %*% weightvec # delta method
  outcov
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
  #' @title Quantizing exposure data
  #' @description Create variables representing indicator functions with cutpoints defined
  #' by quantiles. Output a list that includes: 1) a dataset that is a copy of data, 
  #' except that the variables whose names are included in the `expnms` variable are 
  #' transformed to their quantized version and 2) an unnamed list of the quantile cutpoints
  #' that are used for each of the variables that were quantized
  #' 
  #' @details This function creates categorical variables in place of the
  #' exposure variables named in 'expnms.' For example, a continuous exposure
  #' 'x1' will be replaced in the output data by another 'x1' that takes on values
  #' 0:(q-1), where, for example, the value 1 indicates that the original x1 value
  #' falls between the first and the second quantile.
  #' @return A list containing the following fields
  #' \describe{
  #' \item{data}{a quantized version of the original dataframe}
  #' \item{breaks}{a list of the quantile cutpoints used to create the quantized variables which
  #' includes a very small number for the minimum and a very large number for the maximum to avoid
  #' causing issues when using these breaks to quantize new data.}
  #' }
  #' @param data a data frame
  #' @param expnms a character vector with the names of  the columns to be
  #' quantized
  #' @param q integer, number of quantiles used in creating quantized variables
  #' @param breaks (optional) list of (equal length) numeric vectors that 
  #' characterize the minimum value of each category for which to 
  #' break up the variables named in expnms. This is an alternative to using 'q'
  #' to define cutpoints.
  #' @concept variance mixtures
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
    if(length(expnms)==1){
      data[, expnms] <- qt(1)
    }else{
      #data[, expnms] <- sapply(seq_len(length(expnms)), qt)
      data[, expnms] <- vapply(seq_len(length(expnms)), qt, rep(0.0, nrow(data)))
    }
    return(list(data=data, breaks=e$retbr))
}

checknames <- function(terms){
  #' @title Check for valid model terms in a qgcomp fit
  #' @description This is an internal function called by \code{\link[qgcomp]{qgcomp}},
  #'  \code{\link[qgcomp]{qgcomp.boot}}, and \code{\link[qgcomp]{qgcomp.noboot}},
  #'  but is documented here for clarity. Generally, users will not need to call
  #'  this function directly. This function tries to determine whether there are 
  #'  non-linear terms in the underlying model, which helps infer whether the 
  #'  appropriate function is called, and whether more explicit function calls
  #'  are needed.
  #' @param terms model terms from attr(terms(modelfunction, data), "term.labels")
  nonlin <- ifelse(sum(grep("\\(|\\:|\\^", terms))>0, TRUE, FALSE)
  if(nonlin){
    return(FALSE)
  }else{
    return(TRUE)
  }
}


msm.fit <- function(f, 
                    qdata, 
                    intvals, 
                    expnms, 
                    rr=TRUE, 
                    main=TRUE, 
                    degree=1, 
                    id=NULL,
                    weights,
                    bayes=FALSE, 
                    MCsize=nrow(qdata), ...){
  #' @title Fitting marginal structural model (MSM) within quantile g-computation 
  #' @description This is an internal function called by \code{\link[qgcomp]{qgcomp}},
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
  #' should be represented via the \code{\link[base]{AsIs}} function
  #' @param qdata a data frame with quantized exposures
  #' @param expnms a character vector with the names of the columns in qdata that represent
  #' the exposures of interest (main terms only!)
  #' @param intvals sequence, the sequence of integer values that the joint exposure 
  #' is 'set' to for estimating the msm. For quantile g-computation, this is just 
  #' 0:(q-1), where q is the number of quantiles of exposure.
  #' @param rr logical, estimate log(risk ratio) (family='binomial' only)
  #' @param main logical, internal use: produce estimates of exposure effect (psi)
  #'  and expected outcomes under g-computation and the MSM
  #' @param degree polynomial bases for marginal model (e.g. degree = 2
  #'  allows that the relationship between the whole exposure mixture and the outcome
  #'  is quadratic. Default=1)
  #' @param id (optional) NULL, or variable name indexing individual units of 
  #' observation (only needed if analyzing data with multiple observations per 
  #' id/cluster)
  #' @param weights "case weights" - passed to the "weight" argument of 
  #' \code{\link[stats]{glm}} or \code{\link[arm]{bayesglm}}
  #' @param bayes use underlying Bayesian model (`arm` package defaults). Results
  #' in penalized parameter estimation that can help with very highly correlated 
  #' exposures. Note: this does not lead to fully Bayesian inference in general, 
  #' so results should be interpreted as frequentist.
  #' @param MCsize integer: sample size for simulation to approximate marginal 
  #'  zero inflated model parameters. This can be left small for testing, but should be as large
  #'  as needed to reduce simulation error to an acceptable magnitude (can compare psi coefficients for 
  #'  linear fits with qgcomp.zi.noboot to gain some intuition for the level of expected simulation 
  #'  error at a given value of MCsize)
  #' @param ... arguments to glm (e.g. family)
  #' @seealso \code{\link[qgcomp]{qgcomp.boot}}, and \code{\link[qgcomp]{qgcomp}}
  #' @concept variance mixtures
  #' @import stats arm
  #' @examples
  #' set.seed(50)
  #' dat <- data.frame(y=runif(200), x1=runif(200), x2=runif(200), z=runif(200))
  #' X <- c('x1', 'x2')
  #' qdat <- quantize(dat, X, q=4)$data
  #' mod <- qgcomp:::msm.fit(f=y ~ z + x1 + x2 + I(x1*x2), 
  #'         expnms = c('x1', 'x2'), qdata=qdat, intvals=1:4, bayes=FALSE)
  #' summary(mod$fit) # outcome regression model
  #' summary(mod$msmfit) # msm fit (variance not valid - must be obtained via bootstrap)
  
    newform <- terms(f, data = qdata)
    nobs = nrow(qdata)
    thecall <- match.call(expand.dots = FALSE)
    names(thecall) <- gsub("qdata", "data", names(thecall))
    names(thecall) <- gsub("f", "formula", names(thecall))
    m <- match(c("formula", "data", "weights", "offset"), names(thecall), 0L)
    hasweights = ifelse("weights" %in% names(thecall), TRUE, FALSE)
    thecall <- thecall[c(1L, m)]
    thecall$drop.unused.levels <- TRUE
    
    thecall[[1L]] <- quote(stats::model.frame)
    thecalle <- eval(thecall, parent.frame())
    if(hasweights){
      qdata$weights <- as.vector(model.weights(thecalle))
    } else qdata$weights = rep(1, nobs)
  
    if(is.null(id)) {
      id <- "id__"
      qdata$id__ <- seq_len(dim(qdata)[1])
    }
    # conditional outcome regression fit
    nidx = which(!(names(qdata) %in% id))
    if(!bayes) fit <- glm(newform, data = qdata,
                          weights=weights, 
                          ...)
    if(bayes){
      requireNamespace("arm")
      fit <- bayesglm(f, data = qdata[,nidx,drop=FALSE],
                      weights=weights, 
                      ...)
    } 
    if(fit$family$family %in% c("gaussian", "poisson")) rr=FALSE
    ### 
    # get predictions (set exposure to 0,1,...,q-1)
    if(is.null(intvals)){
      intvals <- (seq_len(length(table(qdata[expnms[1]])))) - 1
    }
    predit <- function(idx, newdata){
      #newdata <- qdata
      newdata[,expnms] <- idx
      suppressWarnings(predict(fit, newdata=newdata, type='response'))
    }
    if(MCsize==nrow(qdata)){
      newdata <- qdata
    }else{
      newids <- data.frame(temp=sort(sample(unique(qdata[,id, drop=TRUE]), MCsize, 
                                            #probs=weights, #bootstrap sampling with weights works with fixed weights, but not time-varying weights
                                            replace = TRUE
      )))
      names(newids) <- id
      newdata <- merge(qdata,newids, by=id, all.x=FALSE, all.y=TRUE)[seq_len(MCsize),]
    }
    predmat = lapply(intvals, predit, newdata=newdata)
    # fit MSM using g-computation estimates of expected outcomes under joint 
    #  intervention
    #nobs <- dim(qdata)[1]
    msmdat <- data.frame(
      cbind(
        Ya = unlist(predmat),
        psi = rep(intvals, each=MCsize),
        weights = rep(newdata$weights, times=length(intvals))
                      #times=length(table(qdata[expnms[1]])))
      )
      )
    # to do: allow functional form variations for the MSM via specifying the model formula
    if(bayes){
      if(!rr) suppressWarnings(msmfit <- bayesglm(Ya ~ poly(psi, degree=degree, raw=TRUE), data=msmdat,
                                                  weights=weights, x=TRUE,
                                                  ...))
      if(rr)  suppressWarnings(msmfit <- bayesglm(Ya ~ poly(psi, degree=degree, raw=TRUE), data=msmdat, 
                                                  family=binomial(link='log'), start=rep(-0.0001, degree+1),
                                                  weights=weights, x=TRUE))
    }
    if(!bayes){
      if(!rr) suppressWarnings(msmfit <- glm(Ya ~ poly(psi, degree=degree, raw=TRUE), data=msmdat,
                                             weights=weights, x=TRUE,
                                             ...))
      if(rr)  suppressWarnings(msmfit <- glm(Ya ~ poly(psi, degree=degree, raw=TRUE), data=msmdat, 
                                             family=binomial(link='log'), start=rep(-0.0001, degree+1),
                                             weights=weights, x=TRUE))
    }
    res <- list(fit=fit, msmfit=msmfit)
    if(main) {
      res$Ya <- msmdat$Ya   # expected outcome under joint exposure, by gcomp
      res$Yamsm <- predict(msmfit, type='response')
      res$Yamsml <- predict(msmfit, type="link")
      res$A <- msmdat$psi # joint exposure (0 = all exposures set category with 
       # upper cut-point as first quantile)
    }
    res
}


qgcomp.noboot <- function(f, 
                          data, 
                          expnms=NULL, 
                          q=4, 
                          breaks=NULL, 
                          id=NULL, 
                          weights, 
                          alpha=0.05, 
                          bayes=FALSE, 
                          ...){
  #' @title Quantile g-computation for continuous, binary, and count outcomes under linearity/additivity
  #'
  #' @description This function estimates a linear dose-response parameter representing a one quantile
  #' increase in a set of exposures of interest. This function is limited to linear and additive
  #' effects of individual components of the exposure. This model estimates the parameters of a marginal 
  #' structural model (MSM) based on g-computation with quantized exposures. Note: this function is  
  #' valid only under linear and additive effects of individual components of the exposure, but when
  #' these hold the model can be fit with very little computational burden.
  #' 
  #' @details For continuous outcomes, under a linear model with no 
  #' interaction terms, this is equivalent to g-computation of the effect of
  #' increasing every exposure by 1 quantile. For binary/count outcomes
  #' outcomes, this yields a conditional log odds/rate ratio(s) representing the 
  #' change in the expected conditional odds/rate (conditional on covariates)
  #' from increasing every exposure by 1 quantile. In general, the latter 
  #' quantity is not equivalent to g-computation estimates. Hypothesis test
  #' statistics and confidence intervals are based on using the delta
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
  #' id/cluster). Note that qgcomp.noboot will not produce cluster-appropriate
  #' standard errors (this parameter is essentially ignored in qgcomp.noboot). 
  #' Qgcomp.boot can be used for this, which will use bootstrap
  #' sampling of clusters/individuals to estimate cluster-appropriate standard
  #' errors via bootstrapping.
  #' @param weights "case weights" - passed to the "weight" argument of 
  #' \code{\link[stats]{glm}} or \code{\link[arm]{bayesglm}}
  #' @param alpha alpha level for confidence limit calculation
  #' @param bayes use underlying Bayesian model (`arm` package defaults). Results
  #' in penalized parameter estimation that can help with very highly correlated 
  #' exposures. Note: this does not lead to fully Bayesian inference in general, 
  #' so results should be interpreted as frequentist.
  #' @param ... arguments to glm (e.g. family)
  #' @seealso \code{\link[qgcomp]{qgcomp.boot}}, and \code{\link[qgcomp]{qgcomp}}
  #' @return a qgcompfit object, which contains information about the effect
  #'  measure of interest (psi) and associated variance (var.psi), as well
  #'  as information on the model fit (fit) and information on the 
  #'  weights/standardized coefficients in the positive (pos.weights) and 
  #'  negative (neg.weights) directions.
  #' @concept variance mixtures
  #' @import stats arm
  #' @export
  #' @examples
  #' set.seed(50)
  #' # linear model
  #' dat <- data.frame(y=runif(50), x1=runif(50), x2=runif(50), z=runif(50))
  #' qgcomp.noboot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, family=gaussian())
  #' # logistic model
  #' dat2 <- data.frame(y=rbinom(50, 1,0.5), x1=runif(50), x2=runif(50), z=runif(50))
  #' qgcomp.noboot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat2, q=2, family=binomial())
  #' # poisson model
  #' dat3 <- data.frame(y=rpois(50, .5), x1=runif(50), x2=runif(50), z=runif(50))
  #' qgcomp.noboot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat3, q=2, family=poisson())
  #' # weighted model
  #' N=5000
  #' dat4 <- data.frame(y=runif(N), x1=runif(N), x2=runif(N), z=runif(N))
  #' dat4$w=runif(N)*2
  #' qdata = quantize(dat4, expnms = c("x1", "x2"))$data
  #' (qgcfit <- qgcomp.noboot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat4, q=4, 
  #'                          family=gaussian(), weights=w))
  #' qgcfit$fit
  #' glm(y ~ z + x1 + x2, data = qdata, weights=w)

  newform <- terms(f, data = data)

  nobs = nrow(data)
  origcall <- thecall <- match.call(expand.dots = FALSE)
  names(thecall) <- gsub("f", "formula", names(thecall))
  m <- match(c("formula", "data", "weights", "offset"), names(thecall), 0L)
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
  
    #expnms <- attr(terms(f, data = data), "term.labels")
    expnms <- attr(newform, "term.labels")
  
    message("Including all model terms as exposures of interest\n")      
  }
  lin = checknames(expnms)
  if(!lin) stop("Model appears to be non-linear: use qgcomp.boot instead")
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
      qdata$id__ = seq_len(dim(qdata)[1])
    }
  
  if(!bayes) fit <- glm(newform, data = qdata,
                        weights=weights,
                        ...)
  if(bayes){
    requireNamespace("arm")
    fit <- bayesglm(newform, data = qdata,
                    weights=weights,
                    ...)
  }
  

    #if(!bayes) fit <- glm(f, data = qdata[,!(names(qdata) %in% id), drop=FALSE],
    #                      #weights=weights,
    #                      ...)
    #if(bayes){
    #  requireNamespace("arm")
    #  fit <- bayesglm(f, data = qdata[,!(names(qdata) %in% id), drop=FALSE],
    #                  #weights=weights,
    #                  ...)
    #}
    mod <- summary(fit)
    if(length(setdiff(expnms, rownames(mod$coefficients)))>0){
      stop("Model aliasing occurred, likely due to perfectly correlated quantized exposures. 
           Try one of the following:
             1) set 'bayes' to TRUE in the qgcomp function (recommended)
             2) set 'q' to a higher value in the qgcomp function (recommended)
             3) check correlation matrix of exposures, and drop all but one variable in each highly correlated set  (not recommended)
           ")
    }
    estb <- c(fit$coefficients[1], sum(mod$coefficients[expnms,1, drop=TRUE]))
    seb <- c(sqrt(mod$cov.scaled[1,1]), se_comb(expnms, covmat = mod$cov.scaled))
    tstat <- estb / seb
    df <- mod$df.null - length(expnms)
    pval <- 2 - 2 * pt(abs(tstat), df = df)
    pvalz <- 2 - 2 * pnorm(abs(tstat))
    ci <- cbind(estb + seb * qnorm(alpha / 2), estb + seb * qnorm(1 - alpha / 2))
    # 'weights'
    wcoef <- fit$coefficients[expnms]
    names(wcoef) <- gsub("_q", "", names(wcoef))
    pos.coef <- which(wcoef > 0)
    neg.coef <- which(wcoef <= 0)
    pos.weights <- abs(wcoef[pos.coef]) / sum(abs(wcoef[pos.coef]))
    neg.weights <- abs(wcoef[neg.coef]) / sum(abs(wcoef[neg.coef]))
    # 'post-hoc' positive and negative estimators 
    # similar to constrained gWQS
    pos.psi <- sum(wcoef[pos.coef])
    neg.psi <- sum(wcoef[neg.coef])
    #nmpos <- names(pos.weights)
    #nmneg <- names(neg.weights)
    #se.pos.psi <- se_comb(nmpos, covmat = mod$cov.scaled)
    #se.neg.psi <- se_comb(nmneg, covmat = mod$cov.scaled)
    qx <- qdata[, expnms]
    names(qx) <- paste0(names(qx), "_q")
    names(estb) <- c('(Intercept)', "psi1")
    res <- list(
      qx = qx, fit = fit, 
      psi = estb[-1], var.psi = seb[-1] ^ 2, covmat.psi=c('psi1' = seb[-1]^2), ci = ci[-1,],
      coef = estb, var.coef = seb ^ 2, 
      #covmat.coef=c('(Intercept)' = seb[1]^2, 'psi1' = seb[2]^2), 
      covmat.coef=vc_comb(aname="(Intercept)", expnms=expnms, covmat = mod$cov.scaled),
      ci.coef = ci,
      expnms=expnms, q=q, breaks=br, degree=1,
      pos.psi = pos.psi, neg.psi = neg.psi,
      pos.weights = sort(pos.weights, decreasing = TRUE),
      neg.weights = sort(neg.weights, decreasing = TRUE), 
      pos.size = sum(abs(wcoef[pos.coef])),
      neg.size = sum(abs(wcoef[neg.coef])),
      bootstrap=FALSE,
      cov.yhat=NULL,
      alpha=alpha,
      call=origcall
    )
      if(fit$family$family=='gaussian'){
        res$tstat <- tstat
        res$df <- df
        res$pval <- pval
      }
      if(fit$family$family %in% c('binomial', 'poisson')){
        res$zstat <- tstat
        res$pval <- pvalz
      }
    attr(res, "class") <- "qgcompfit"
    res
}

qgcomp.boot <- function(f, 
                        data, 
                        expnms=NULL, 
                        q=4, 
                        breaks=NULL, 
                        id=NULL, 
                        weights, 
                        alpha=0.05, 
                        B=200, 
                        rr=TRUE, 
                        degree=1, 
                        seed=NULL, 
                        bayes=FALSE, 
                        MCsize=nrow(data), 
                        parallel=FALSE, ...){
  #' @title Quantile g-computation for continuous and binary outcomes 
  #'  
  #' @description This function estimates a linear dose-response parameter representing a one quantile
  #' increase in a set of exposures of interest. This model estimates the parameters of a marginal 
  #' structural model (MSM) based on g-computation with quantized exposures. Note: this function  
  #' allows linear and non-additive effects of individual components of the exposure, as well as
  #' non-linear joint effects of the mixture via polynomial basis functions, which increase the
  #' computational computational burden due to the need for non-parametric bootstrapping.
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
  #' id/cluster). Note that qgcomp.noboot will not produce cluster-appropriate
  #' standard errors. Qgcomp.boot can be used for this, which will use bootstrap
  #' sampling of clusters/individuals to estimate cluster-appropriate standard
  #' errors via bootstrapping.
  #' @param weights "case weights" - passed to the "weight" argument of 
  #' \code{\link[stats]{glm}} or \code{\link[arm]{bayesglm}}
  #' @param alpha alpha level for confidence limit calculation
  #' @param B integer: number of bootstrap iterations (this should typically be >=200,
  #'  though it is set lower in examples to improve run-time).
  #' @param rr logical: if using binary outcome and rr=TRUE, qgcomp.boot will 
  #'   estimate risk ratio rather than odds ratio
  #' @param degree polynomial bases for marginal model (e.g. degree = 2
  #'  allows that the relationship between the whole exposure mixture and the outcome
  #'  is quadratic (default = 1).
  #' @param seed integer or NULL: random number seed for replicable bootstrap results
  #' @param bayes use underlying Bayesian model (`arm` package defaults). Results
  #' in penalized parameter estimation that can help with very highly correlated 
  #' exposures. Note: this does not lead to fully Bayesian inference in general, 
  #' so results should be interpreted as frequentist.
  #' @param MCsize integer: sample size for simulation to approximate marginal 
  #'  zero inflated model parameters. This can be left small for testing, but should be as large
  #'  as needed to reduce simulation error to an acceptable magnitude (can compare psi coefficients for 
  #'  linear fits with qgcomp.noboot to gain some intuition for the level of expected simulation 
  #'  error at a given value of MCsize). This likely won't matter much in linear models, but may 
  #'  be important with binary or count outcomes.
  #' @param parallel use (safe) parallel processing from the future and future.apply packages
  #' @param ... arguments to glm (e.g. family)
  #' @seealso \code{\link[qgcomp]{qgcomp.noboot}}, and \code{\link[qgcomp]{qgcomp}}
  #' @return a qgcompfit object, which contains information about the effect
  #'  measure of interest (psi) and associated variance (var.psi), as well
  #'  as information on the model fit (fit) and information on the 
  #'  marginal structural model (msmfit) used to estimate the final effect
  #'  estimates.
  #' @concept variance mixtures
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
  #'  \dontrun{
  #' qgcomp.boot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=4, 
  #'   family=gaussian(), B=200) # B should be at least 200 in actual examples
  #'   
  #'  # Note that these give different answers! In the first, the estimate is conditional on Z,
  #'  # but in the second, Z is marginalized over via standardization. The estimates
  #'  # can be made approximately the same by centering Z (for linear models), but
  #'  # the conditional estimate will typically have lower standard errors.
  #'  dat$z = dat$z - mean(dat$z)
  #'   
  #' # Conditional linear slope
  #' qgcomp.noboot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=4, family=gaussian())
  #' # Marginal linear slope (population average slope, for a purely linear, 
  #' #  additive model this will equal the conditional)
  #'  
  #' qgcomp.boot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=4, 
  #'   family=gaussian(), B=200) # B should be at least 200 in actual examples
  #'  
  #' # Population average mixture slope which accounts for non-linearity and interactions
  #' qgcomp.boot(y ~ z + x1 + x2 + I(x1^2) + I(x2*x1), family="gaussian", 
  #'  expnms = c('x1', 'x2'), data=dat, q=4, B=200)
  #'  
  #' # generally non-linear/non-addiive underlying models lead to non-linear mixture slopes
  #' qgcomp.boot(y ~ z + x1 + x2 + I(x1^2) + I(x2*x1), family="gaussian", 
  #'  expnms = c('x1', 'x2'), data=dat, q=4, B=200, deg=2)
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
  #'   data=dat, q=2, B=3, rr=FALSE)
  #'   
  #' # Population average mixture RR
  #' qgcomp.boot(y ~ z + x1 + x2, family="binomial", expnms = c('x1', 'x2'), 
  #'   data=dat, q=2, rr=TRUE, B=3)
  #'   
  #' # Population average mixture RR, indicator variable representation of x2
  #' # note that I(x==...) operates on the quantile-based category of x,
  #' # rather than the raw value
  #' res = qgcomp.boot(y ~ z + x1 + I(x2==1) + I(x2==2) + I(x2==3), 
  #'   family="binomial", expnms = c('x1', 'x2'), data=dat, q=4, rr=TRUE, B=200)
  #' res$fit  
  #' plot(res)
  #' 
  #' # now add in a non-linear MSM
  #' res2 = qgcomp.boot(y ~ z + x1 + I(x2==1) + I(x2==2) + I(x2==3), 
  #'   family="binomial", expnms = c('x1', 'x2'), data=dat, q=4, rr=TRUE, B=200, 
  #'   degree=2)
  #' res2$fit  
  #' res2$msmfit  # correct point estimates, incorrect standard errors
  #' res2  # correct point estimates, correct standard errors
  #' plot(res2)
  #' # Log risk ratio per one IQR change in all exposures (not on quantile basis)
  #' dat$x1iqr <- dat$x1/with(dat, diff(quantile(x1, c(.25, .75))))
  #' dat$x2iqr <- dat$x2/with(dat, diff(quantile(x2, c(.25, .75))))
  #' # note that I(x>...) now operates on the untransformed value of x,
  #' # rather than the quantized value
  #' res2 = qgcomp.boot(y ~ z + x1iqr + I(x2iqr>0.1) + I(x2>0.4) + I(x2>0.9), 
  #'   family="binomial", expnms = c('x1iqr', 'x2iqr'), data=dat, q=NULL, rr=TRUE, B=200, 
  #'   degree=2)
  #' res2
  #' # using parallel processing
  #'
  #' qgcomp.boot(y ~ z + x1iqr + I(x2iqr>0.1) + I(x2>0.4) + I(x2>0.9), 
  #'   family="binomial", expnms = c('x1iqr', 'x2iqr'), data=dat, q=NULL, rr=TRUE, B=200, 
  #'   degree=2, parallel=TRUE)
  #'
  #'
  #' # weighted model
  #' N=5000
  #' dat4 <- data.frame(id=seq_len(N), x1=runif(N), x2=runif(N), z=runif(N))
  #' dat4$y <- with(dat4, rnorm(N, x1*z + z, 1))
  #' dat4$w=runif(N) + dat4$z*5
  #' qdata = quantize(dat4, expnms = c("x1", "x2"), q=4)$data
  #' # first equivalent models with no covariates
  #' qgcomp.noboot(f=y ~ x1 + x2, expnms = c('x1', 'x2'), data=dat4, q=4, family=gaussian())
  #' qgcomp.noboot(f=y ~ x1 + x2, expnms = c('x1', 'x2'), data=dat4, q=4, family=gaussian(), 
  #'               weights=w)
  #' 
  #' set.seed(13)
  #' qgcomp.boot(f=y ~ x1 + x2, expnms = c('x1', 'x2'), data=dat4, q=4, family=gaussian(), 
  #'             weights=w)
  #' # using the correct model
  #' set.seed(13)
  #' qgcomp.boot(f=y ~ x1*z + x2, expnms = c('x1', 'x2'), data=dat4, q=4, family=gaussian(), 
  #'             weights=w, id="id")
  #' (qgcfit <- qgcomp.boot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat4, q=4, 
  #'                        family=gaussian(), weights=w))
  #' qgcfit$fit
  #' summary(glm(y ~ z + x1 + x2, data = qdata, weights=w))
  #' }
  # character names of exposure mixture components
    oldq = NULL
    if(is.null(seed)) seed = round(runif(1, min=0, max=1e8))
  
    newform <- terms(f, data = data)
    class(newform) <- "formula"
  
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
      # TODO: allow user specification of this
      nvals = length(table(unlist(data[,expnms])))
      if(nvals < 10){
        message("\nNote: using all possible values of exposure as the 
              intervention values\n")
        p = length(expnms)
        intvals <- as.numeric(names(table(unlist(data[,expnms]))))
        
        br <- lapply(seq_len(p), function(x) c(-1e16, intvals[2:nvals]-1e-16, 1e16))
      }else{
    message("\nNote: using quantiles of all exposures combined in order to set 
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
    ###
    msmfit <- msm.fit(newform, qdata, intvals, expnms, rr, main=TRUE,degree=degree, id=id, 
                      weights, 
                      bayes, 
                      MCsize=MCsize, 
                      ...)
    # main estimate  
    #estb <- as.numeric(msmfit$msmfit$coefficients[-1])
    estb <- as.numeric(msmfit$msmfit$coefficients)
    #bootstrap to get std. error
    nobs <- dim(qdata)[1]
    nids <- length(unique(qdata[,id, drop=TRUE]))
    starttime = Sys.time()
    psi.only <- function(i=1, f=f, qdata=qdata, intvals=intvals, expnms=expnms, rr=rr, degree=degree, 
                         nids=nids, id=id,
                         weights,MCsize=MCsize,
                         ...){
      if(i==2 & !parallel){
        timeiter = as.numeric(Sys.time() - starttime)
        if((timeiter*B/60)>0.5) message(paste0("Expected time to finish: ", round(B*timeiter/60, 2), " minutes \n"))
      }
      bootids <- data.frame(temp=sort(sample(unique(qdata[,id, drop=TRUE]), nids, 
                                             replace = TRUE
                                             )))
      names(bootids) <- id
      qdata_ <- merge(qdata,bootids, by=id, all.x=FALSE, all.y=TRUE)
      ft = msm.fit(f, qdata_, intvals, expnms, rr, main=FALSE, degree, id, weights=weights, bayes, MCsize=MCsize,
                   ...)
      yhatty = data.frame(yhat=predict(ft$msmfit), psi=ft$msmfit$data[,"psi"])
      as.numeric(
        # the yhat estimates will be identical across individuals due to this being a marginal model
        c(ft$msmfit$coefficients, with(yhatty, tapply(yhat, psi, mean)))
      )
    }
    set.seed(seed)
    if(parallel){
      Sys.setenv(R_FUTURE_SUPPORTSMULTICORE_UNSTABLE="quiet")
      future::plan(strategy = future::multisession)
      bootsamps <- future.apply::future_sapply(X=seq_len(B), FUN=psi.only,f=f, qdata=qdata, intvals=intvals, 
      #bootsamps <- future.apply::future_vapply(X=seq_len(B), FUN=psi.only,f=f, qdata=qdata, intvals=intvals, 
                          expnms=expnms, rr=rr, degree=degree, nids=nids, id=id,
                          weights=qdata$weights,MCsize=MCsize,
                          future.seed=TRUE,
                          ...)
      
      future::plan(strategy = future::transparent)
    }else{
      bootsamps <- sapply(X=seq_len(B), FUN=psi.only,f=f, qdata=qdata, intvals=intvals, 
      #bootsamps <- vapply(X=seq_len(B), FUN=psi.only,f=f, qdata=qdata, intvals=intvals, 
                          expnms=expnms, rr=rr, degree=degree, nids=nids, id=id,
                          weights=weights, MCsize=MCsize,
                          ...)
      
    }
    # these are the linear predictors
    hats = t(bootsamps[-c(seq_len(degree+1)),])
    # covariance of the linear predictors
    cov.yhat = cov(hats)
    bootsamps = bootsamps[seq_len(degree+1),]
    seb <- apply(bootsamps, 1, sd)
    covmat <- cov(t(bootsamps))
    colnames(covmat) <- rownames(covmat) <- names(estb) <- c("(intercept)", paste0("psi", seq_len(nrow(bootsamps)-1)))

    tstat <- estb / seb
    df <- nobs - length(attr(terms(f, data = data), "term.labels")) - 1 - degree # df based on obs - gcomp terms - msm terms
    pval <- 2 - 2 * pt(abs(tstat), df = df)
    pvalz <- 2 - 2 * pnorm(abs(tstat))
    ci <- cbind(estb + seb * qnorm(alpha / 2), estb + seb * qnorm(1 - alpha / 2))
    # outcome 'weights' not applicable in this setting, generally (i.e. if using this function for non-linearity, 
    #   then weights will vary with level of exposure)
    if (!is.null(oldq)){
      q = oldq
    }
    qx <- qdata[, expnms]
    res <- list(
      qx = qx, fit = msmfit$fit, msmfit = msmfit$msmfit, 
      psi = estb[-1], var.psi = seb[-1] ^ 2, covmat.psi=covmat[-1,-1, drop=FALSE], ci = ci[-1,],
      coef = estb, var.coef = seb ^ 2, covmat.coef=covmat, ci.coef = ci,
      expnms=expnms, q=q, breaks=br, degree=degree,
      pos.psi = NULL, neg.psi = NULL, 
      pos.weights = NULL, neg.weights = NULL, pos.size = NULL,neg.size = NULL, bootstrap=TRUE,
      y.expected=msmfit$Ya, y.expectedmsm=msmfit$Yamsm, index=msmfit$A,
      bootsamps = bootsamps,
      cov.yhat=cov.yhat,
      alpha=alpha,
      call=origcall
    )
      if(msmfit$fit$family$family=='gaussian'){
        res$tstat <- tstat
        res$df <- df
        res$pval <- pval
      }
      if(msmfit$fit$family$family %in% c('binomial', 'poisson')){
        res$zstat <- tstat
        res$pval <- pvalz
      }
    attr(res, "class") <- "qgcompfit"
    res
}


qgcomp <- function(f,data=data,family=gaussian(),rr=TRUE,...){
  #' @title Quantile g-computation for continuous, binary, count, and censored survival outcomes
  #' 
  #' @description This function automatically selects between qgcomp.noboot, qgcomp.boot,
  #'  qgcomp.cox.noboot, and qgcomp.cox.boot
  #'  for the most efficient approach to estimate the average expected 
  #'  change in the (log) outcome per quantile increase in the joint 
  #'  exposure to all exposures in `expnms', given the underlying model. For example,
  #'  if the underlying model (specified by the formula `f`) is a linear model with all
  #'  linear terms for exposure, then `qgcomp.noboot`` will be called to fit the model. Non-linear
  #'  terms or requesting the risk ratio for binomial outcomes will result in the `qgcomp.boot`
  #'  function being called. For a given linear model, boot and noboot versions will give identical
  #'  inference, though when using survival outcomes, the `boot` version uses simulation based 
  #'  inference, which can vary from the `noboot` version due to simulation error (which can 
  #'  be minimized via setting the MCsize parameter very large - see 
  #'  \code{\link[qgcomp]{qgcomp.cox.boot}} for details).
  #'
  #' @param f R style formula (may include survival outcome via \code{\link[survival]{Surv}})
  #' @param data data frame
  #' @param family \code{gaussian()}, \code{binomial()}, \code{cox()}, \code{poisson()} (works as
  #' argument to 'family' parameter in \code{\link[stats]{glm}}` or 'dist' parameter in 
  #' \code{\link[pscl]{zeroinfl}})
  #' @param rr logical: if using binary outcome and rr=TRUE, qgcomp.boot will 
  #' estimate risk ratio rather than odds ratio. Note, to get population average 
  #' effect estimates for a binary outcome, set rr=TRUE (default: ORs are generally not
  #' of interest as population average effects, so if rr=FALSE then a conditional
  #' OR will be estimated, which cannot be interpreted as a population average
  #' effect
  #' @param ... arguments to qgcomp.noboot or qgcomp.boot (e.g. q) or glm
  #' @seealso \code{\link[qgcomp]{qgcomp.noboot}}, \code{\link[qgcomp]{qgcomp.boot}}, 
  #'  \code{\link[qgcomp]{qgcomp.cox.noboot}}, \code{\link[qgcomp]{qgcomp.cox.boot}}
  #'  \code{\link[qgcomp]{qgcomp.zi.noboot}} and \code{\link[qgcomp]{qgcomp.zi.boot}}
  #'  (\code{\link[qgcomp]{qgcomp}} is just a wrapper for these functions)
  #' @return a qgcompfit object, which contains information about the effect
  #'  measure of interest (psi) and associated variance (var.psi), as well
  #'  as information on the model fit (fit) and possibly information on the 
  #'  marginal structural model (msmfit) used to estimate the final effect
  #'  estimates (qgcomp.boot, qgcomp.cox.boot only). 
  #'  If appropriate, weights are also reported, which represent the proportion
  #'  of a directional (positive/negative) effect that is accounted for by
  #'  each exposure.
  #' @concept variance mixtures
  #' @import stats
  #' @export
  #' @examples
  #' set.seed(50)
  #' dat <- data.frame(y=runif(50), x1=runif(50), x2=runif(50), z=runif(50))
  #' qgcomp.noboot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2)
  #' qgcomp.boot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, B=10, seed=125)
  #' # automatically selects appropriate method
  #' qgcomp(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2)
  #' # note for binary outcome this will choose the risk ratio (and bootstrap methods) by default
  #' dat <- data.frame(y=rbinom(100, 1, 0.5), x1=runif(100), x2=runif(100), z=runif(100))
  #' \dontrun{
  #' qgcomp.noboot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, family=binomial())
  #' set.seed(1231)
  #' qgcomp.boot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, family=binomial())
  #' set.seed(1231)
  #' qgcomp(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, family=binomial())
  #'   
  #' # automatically selects appropriate method when specifying rr or degree explicitly
  #' qgcomp(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, family=binomial(), rr=FALSE)
  #' qgcomp(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, family=binomial(), rr=TRUE)
  #' qgcomp(y ~ z + factor(x1) + factor(x2), degree=2, expnms = c('x1', 'x2'), data=dat, q=4, 
  #' family=binomial())
  # 
  #' 
  #' #survival objects
  #' set.seed(50)
  #' N=200
  #' dat <- data.frame(time=(tmg <- pmin(.1,rweibull(N, 10, 0.1))), 
  #'                 d=1.0*(tmg<0.1), x1=runif(N), x2=runif(N), z=runif(N))
  #' expnms=paste0("x", 1:2)
  #' f = survival::Surv(time, d)~x1 + x2
  #' qgcomp(f, expnms = expnms, data = dat)
  #' # note if B or MCsize are set but the model is linear, an error will result
  #' try(qgcomp(f, expnms = expnms, data = dat, B1=, MCsize))
  #' # note that in the survival models, MCsize should be set to a large number
  #' #  such that results are repeatable (within an error tolerance such as 2 significant digits)
  #' # if you run them under different  seed values
  #' f = survival::Surv(time, d)~x1 + x2 + x1:x2
  #' qgcomp(f, expnms = expnms, data = dat, B=10, MCsize=100)
  #' }
  requireNamespace("survival")
  iszip = (length(grep("\\|", as.character(f)))>0)
  issurv = survival::is.Surv(eval(attr(terms(f, data = data), "variables")[[2]], envir = data))
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    if(issurv){
      message("Survival model\n")
    }else{
      print(family)
      stop("'family' not recognized")
    }
  }
  if(!(family$family %in% c("binomial")) & rr) {
    #warning("'rr=TRUE' is for bimomial family only, setting rr=FALSE")
    rr = FALSE
  }
  terms <- attr(terms(f,data=data), 'term.labels')
  doboot <- !checknames(terms)
  if(issurv & doboot ){
    res <- qgcomp.cox.boot(f=f,data=data,...)
  } else if(issurv & !doboot ){
    res <- qgcomp.cox.noboot(f=f,data=data,...)
  } else if(!iszip & (rr | doboot)){
    res <- qgcomp.boot(f=f,data=data,family=family,rr=rr,...)
  }else if(!iszip){
    res <- qgcomp.noboot(f=f,data=data,family=family,...)
  }else if(iszip & doboot){
    res <- qgcomp.zi.boot(f=f,data=data,dist=family$family,rr=rr,...)
  }else if(iszip & !doboot){
    res <- qgcomp.zi.noboot(f=f,data=data,dist=family$family,...)
  } else{
    stop('Unable to parse the model type: try one of the specific functions in the qgcomp package
           e.g. qgcomp.noboot, qgcomp.boot, qgcomp.cox.noboot, qgcomp.cox.boot, qgcomp.zi.noboot, 
         qgcomp.zi.boot, ')
  }
  res
}






msm.predict <- function(object, newdata=NULL){
  #' @title Secondary prediction method for the (non-survival) qgcomp MSM.
  #' 
  #' @description this is an internal function called by 
  #'  \code{\link[qgcomp]{qgcomp.boot}},
  #'  but is documented here for clarity. Generally, users will not need to call
  #'  this function directly.
  #' @description Get predicted values from a qgcompfit object from
  #' \code{\link[qgcomp]{qgcomp.boot}}.
  #'
  #' @details (Not usually called by user) Makes predictions from the MSM 
  #' (rather than the conditional g-computation 
  #' fit) from a "qgcompfit" object. Generally, this should not be used in 
  #' favor of the default \code{\link[qgcomp]{predict.qgcompfit}} function. This
  #' function can only be used following the `qgcomp.boot` function. For the 
  #' `qgcomp.noboot` function, \code{\link[qgcomp]{predict.qgcompfit}} gives 
  #' identical inference to predicting from an MSM.
  #'
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

