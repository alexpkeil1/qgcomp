.devinstall <- function(...){
  .qgc.require("devtools")
  devtools::install_github("alexpkeil1/qgcomp",...)
}


.fold_list <- function (fold, set1, set2) {
  fold_list <- list(fold = fold, set1 = set1, set2 = set2)
  fold_list
}

.xfitfold_list_from_foldvec <- function (fold, folds, ordermat) {
  nfolds <- length(unique(folds))
  set1 <- which(folds %in% ordermat[1:(nfolds - 1), fold])
  set2 <- which(folds == ordermat[nfolds, fold])
  .fold_list(fold, set1, set2)
}

.make_xfitfolds_iid <- function (n, V = 5) {
  folds <- rep(seq_len(V), length = n)
  folds <- sample(folds)
  combinations <- utils::combn(V, V - 1)
  combinations <- rbind(combinations, apply(combinations, 2, 
                                            function(x) setdiff(1:V, x)))
  if (V > 1) 
    foldobj = lapply(1:V, .xfitfold_list_from_foldvec, 
                     folds = folds, ordermat = combinations)
  if (V == 1) 
    foldobj = list(.fold_list(fold = 1, set1 = 1:n, 
                                           set2 = 1:n))
  foldobj
}


.xfit_grab <- function(foldres, stat,initval=0,whichval=1){
  statmat = matrix(initval, nrow=length(foldres), ncol=2)
  for(i in seq_len(length(foldres))){
    res = foldres[[i]]
    if(!is.null(res$pos.fit)){
      statmat[i,1] = unlist(res$pos.fit[stat])[whichval]
    }
    if(!is.null(res$neg.fit)){
      statmat[i,2] = unlist(res$neg.fit[stat])[whichval]
    }
  }
  statmat
}

.xfit_procfolds <- function(foldres){
  list(
    intercepts = .xfit_grab(foldres, "coef"),
    psis = .xfit_grab(foldres, "coef",whichval=2),
    vars_intercept = .xfit_grab(foldres, "var.coef"),
    vars_psi = .xfit_grab(foldres, "var.coef",whichval=2)
  )
}

.xfit_proclist <- function(proclist,n){
  int_ests = apply(proclist$intercepts, 2, median)
  psi_ests = apply(proclist$psis, 2, median)
  int_resids = t(t(proclist$intercepts)-int_ests)
  psi_resids = t(t(proclist$psis)-psi_ests)
  int_vars = proclist$vars_intercept/n + int_resids^2
  psi_vars = proclist$vars_psi/n + psi_resids^2
  c(proclist, list(int_ests=int_ests,
                   psi_ests=psi_ests,
                   int_resids =int_resids ,
                   psi_resids =psi_resids ,
                   int_vars = apply(int_vars, 2, median),
                   psi_vars = apply(psi_vars, 2, median)
                   )
  )
}
  
.calcxfitci <- function(xdf, alpha=0.05){
  xdf =  cbind(xdf, xdf[,1] + qnorm(alpha/2)* xdf[,2])
  xdf =  cbind(xdf, xdf[,1] + qnorm(1-alpha/2)* xdf[,2])
  xdf =  cbind(xdf, xdf[,1]/xdf[,2])
  xdf =  cbind(xdf, 2 - 2 * pnorm(abs(xdf[,1]/xdf[,2])))
  colnames(xdf) <- c("Estimate", "Std. Err", "Lower CI", "Upper CI", "z value", "Pr(>|t|)")
  rownames(xdf) <- c("(Intercept)", "psi")
  xdf
}

.xfit_coeftab <- function(x){
  x$neg.coeftab = (rbind(
    c(x$int_ests[2],sqrt(x$int_vars[2])),
    c(x$psi_ests[2],sqrt(x$psi_vars[2]))
  ))
  x$pos.coeftab = (rbind(
    c(x$int_ests[1], sqrt(x$int_vars[1])),
    c(x$psi_ests[1], sqrt(x$psi_vars[1]))
  ))
  x$neg.coeftab = .calcxfitci(x$neg.coeftab)
  x$pos.coeftab = .calcxfitci(x$pos.coeftab)
  x
}


.qgcomp_xfitpartials <- function(
    data,
    fun = c("qgcomp.glm.noboot", "qgcomp.cox.noboot", "qgcomp.zi.noboot"),
    V=10,
    expnms=NULL,
    .fixbreaks=TRUE,
    ...
){
  n = nrow(data)
  folds <- .make_xfitfolds_iid(n=n, V=V)
  foldres = list()
  for(f in folds){
    traindata = data[f$set1,,drop=FALSE]
    validdata = data[f$set2,,drop=FALSE]
    foldres[[f$fold]] = qgcomp.partials(
      fun=fun,
      expnms = expnms, 
      traindata = traindata, 
      validdata = validdata,
      .fixbreaks = .fixbreaks,
      ...
    )
  }
  res = .xfit_procfolds(foldres)
  res = .xfit_proclist(res, n)
  res = .xfit_coeftab(res)
  alpha = foldres[[1]]$pos.fit$alpha
  res = c(res, list(foldres = foldres, n=n))
  class(res) <-  "qgcompmultixfit"
  res
}


print.qgcompmultixfit <- function(x,...){
  #' @export
  #cat(paste0("\nVariables with positive effect sizes in training data: ", paste(x$posmix, collapse = ", ")))
  #cat(paste0("\nVariables with negative effect sizes in training data: ", paste(x$negmix, collapse = ", ")))
  cat("\nPartial effect sizes estimated using V-fold cross-fitting\n")
  cat("Positive direction \n")
  if(!is.null(x$pos.coeftab)){
    printCoefmat(x$pos.coeftab, has.Pvalue = TRUE)
  } else{
    cat("\n No positive coefficients in model fit to training data ")
  }
  cat("\nNegative direction \n")
  if(!is.null(x$neg.coeftab)){
    printCoefmat(x$neg.coeftab, has.Pvalue = TRUE)
  } else{
    cat("\n No negative coefficients in model fit to training data ")
  } 
  
}


#' @title Quantile g-computation for survival outcomes in a case-cohort design under linearity/additivity
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
#' @param subcoh (From \code{\link[survival]{cch}} help) Vector of indicators for subjects sampled as part of the sub-cohort. Code 1 or TRUE for members of the sub-cohort, 0 or FALSE for others. If data is a data frame then subcoh may be a one-sided formula.
#' @param id (From \code{\link[survival]{cch}} help) Vector of unique identifiers, or formula specifying such a vector.
#' @param cohort.size (From \code{\link[survival]{cch}} help) Vector with size of each stratum original cohort from which subcohort was sampled
#' @param expnms character vector of exposures of interest
#' @param q NULL or number of quantiles used to create quantile indicator variables
#' representing the exposure variables. If NULL, then gcomp proceeds with un-transformed
#' version of exposures in the input datasets (useful if data are already transformed,
#' or for performing standard g-computation)
#' @param breaks (optional) NULL, or a list of (equal length) numeric vectors that 
#' characterize the minimum value of each category for which to 
#' break up the variables named in expnms. This is an alternative to using 'q'
#' to define cutpoints.
#' @param weights Not used here
#' @param cluster not yet implemented
#' @param alpha alpha level for confidence limit calculation
#' @param ... arguments to glm (e.g. family)
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
qgcomp.cch.noboot <- function(f, data, subcoh=NULL, id=NULL, cohort.size=NULL, expnms = NULL, q = 4, breaks = NULL, 
                              weights, cluster = NULL, alpha = 0.05, ...) 
{
  of <- f
  newform <- terms(f, data = data)
  class(newform) <- "formula"
  nobs = nrow(data)
  origcall <- thecall <- match.call(expand.dots = FALSE)
  names(thecall) <- gsub("f", "formula", names(thecall))
  m <- match(c("f", "data", "weights", "offset"), names(thecall), 
             0L)
  hasweights = ifelse("weights" %in% names(thecall), TRUE, 
                      FALSE)
  thecall <- thecall[c(1L, m)]
  thecall$drop.unused.levels <- TRUE
  thecall[[1L]] <- quote(stats::model.frame)
  thecalle <- eval(thecall, parent.frame())
  if (hasweights) {
    data$weights <- as.vector(model.weights(thecalle))
  }
  else data$weights = rep(1, nobs)
  if (is.null(expnms)) {
    message("Including all model terms as exposures of interest")
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
  environment(newform) <- list2env(list(qdata = qdata))
  fit <- cch(newform, data = qdata, subcoh = subcoh, id=id, cohort.size=cohort.size, ...)
  cchfam = list(family = "cch", link = "log", linkfun = log)
  class(cchfam) = "family"
  fit[["family"]] = cchfam
  mod <- summary(fit)
  estb <- sum(mod$coefficients[expnms, 1])
  covMat = fit$var
  colnames(covMat) <- names(mod$coefficients[, 1])
  seb <- se_comb(expnms, covmat = covMat)
  tstat <- estb/seb
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
  qx <- qdata[, expnms]
  names(qx) <- paste0(names(qx), "_q")
  names(estb) = "psi1"
  res <- .qgcomp_object(qx = qx, fit = fit, psi = estb, var.psi = seb^2, 
                                 covmat.psi = seb^2, ci = ci, coef = estb, var.coef = seb^2, 
                                 covmat.coef = seb^2, ci.coef = ci, expnms = expnms, q = q, 
                                 breaks = br, degree = 1, pos.psi = pos.psi, neg.psi = neg.psi, 
                                 pos.weights = sort(pos.weights, decreasing = TRUE), neg.weights = sort(neg.weights, 
                                                                                                        decreasing = TRUE), pos.size = sum(abs(wcoef[poscoef])), 
                                 neg.size = sum(abs(wcoef[negcoef])), zstat = tstat, pval = pvalz, 
                                 alpha = alpha, call = origcall, hasintercept = FALSE)
  attr(res, "class") <- c("survqgcompfit", attr(res, "class"))
  res
}

# in progress plotting function
.plot_boot_multinomial_base <- function(r, p, x, modelband, flexfit, modelfitline, pointwisebars, pointwiseref=1, alpha=0.05){
  #
  x$fit
  x$msmfit
  stop("Not yet implemented")
  p <- p #+ #scale_y_log10()
    p <- p + labs(x = "Joint exposure quantile", y = "P(Y=y)") + lims(x=c(0,1))
    if(modelband) p <- p + .plot.or.mod.bounds(x,alpha)
    if(flexfit)   p <- p + .plot.or.smooth.line(x)
    if(modelfitline) p <- p + .plot.logitlin.line(x)
    if(pointwisebars) p <- p + .plot.or.pw.boot(x,alpha,pointwiseref)
  p
}

