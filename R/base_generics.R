# generics from other packages

coef.qgcompfit <- function(object, ...){
  #' @exportS3Method stats::coef
  object$coef
}

coef.eefit <- function(object,...){
  #' @exportS3Method stats::coef
  object$est
}

deviance.qgcompfit <- function(object, ...){
  #' @exportS3Method stats::deviance
  object$fit$deviance
}


df.residual.qgcompfit <- function(object, ...){
  #' @exportS3Method stats::df.residual
  df.residual(object$fit)
}


model.matrix.eeqgcompfit <- function(object, ...) {
  #' @exportS3Method stats::model.matrix
  object$fit$X
}

formula.eeqgcompfit <- function(x, ...) {
  #' @exportS3Method stats::formula
  x$call$f
}

vcov.qgcompfit <- function(object, ...){
  #' @exportS3Method stats::vcov
  object$covmat.coef
}

vcov.eefit <- function(object,...){
  #' @exportS3Method stats::vcov
  # a fit or msmfit object from ee methods
  object$vcov
}

AIC.qgcompfit <- function(object, ...){
  #' @exportS3Method stats::AIC
  if(!is.null(object$bread)) stop("Calculation not available for estimating equation methods")
  AIC(object$fit)
}

BIC.qgcompfit <- function(object, ...){
  #' @exportS3Method stats::BIC
  if(!is.null(object$bread)) stop("Calculation not available for estimating equation methods")
  BIC(object$fit)
}

logLik.qgcompfit <- function(object, ...){
  #' @exportS3Method stats::logLik
  if(!is.null(object$bread)) stop("Calculation not available for estimating equation methods")
  logLik(object$fit)
}

anova.qgcompfit <- function(object, ...){
  #' @exportS3Method stats::anova
  if(!is.null(object$bread)) stop("Calculation not available for estimating equation methods")
  anova(object$fit)
}

anova.eeqgcompfit = function(object, ..., dispersion = NULL, test = NULL)
{
  #' @exportS3Method stats::anova
  # based on geepack:::anova.geeglm
  dotargs <- list(...)
  named <- if (is.null(names(dotargs)))
    rep(FALSE, length(dotargs))
  else (names(dotargs) != "")
  if (any(named))
    warning("The following arguments to anova(..) are invalid and dropped: ",
            paste(deparse(dotargs[named]), collapse = ", "))
  dotargs <- dotargs[!named]
  is.eefit <- unlist(lapply(dotargs, function(x) inherits(x, "eeqgcompfit")))
  dotargs <- dotargs[is.eefit]
  if (length(dotargs) > 0)
    return(anova.eeqgcompfitlist(c(list(object), dotargs), dispersion = dispersion,
                                 test = test))
  else{
    stop("Two model fits are needed")
  }
}

anova.eeqgcompfitlist <- function(object, ..., dispersion = NULL, test = NULL)
{
  #' @exportS3Method stats::anova
  responses <- as.character(lapply(object, function(x) {
    deparse(formula(x)[[2]])
  }))
  sameresp <- responses == responses[1]
  if (!all(sameresp)) {
    object <- object[sameresp]
    warning("Models with response ", deparse(responses[!sameresp]),
            " removed because response differs from ", "model 1")
  }
  ns <- sapply(object, function(x) length(x$y.expected))
  if (any(ns != ns[1]))
    stop("models were not all fitted to the same size of dataset")
  objects <- list(object, ...)
  m1 <- objects[[1]][[1]]
  if (length(objects[[1]]) > 1)
    m2 <- objects[[1]][[2]]
  else m2 <- NULL
  value <- anovaqgcompgee(m1, m2)
  return(value)
}

anovaqgcompgee <- function(m1, m2, ...) {
  mm1 <- model.matrix(m1)
  mm2 <- model.matrix(m2)
  P1 <- mm1 %*% solve(t(mm1) %*% mm1) %*% t(mm1)
  P2 <- mm2 %*% solve(t(mm2) %*% mm2) %*% t(mm2)
  e2 <- mm2 - P1 %*% mm2
  e1 <- mm1 - P2 %*% mm1
  m2inm1 <- all(apply(e2, 2, var) < 1e-10)
  m1inm2 <- all(apply(e1, 2, var) < 1e-10)
  if (!any(c(m2inm1, m1inm2)))
    cat("Models not nested\n")
  else if (all(c(m2inm1, m1inm2)))
    cat("Models are identical\n")
  else {
    if (m1inm2) {
      tmp <- m1
      m1 <- m2
      m2 <- tmp
    }
    mm1 <- model.matrix(m1)
    mm2 <- model.matrix(m2)
    m1emm = m1$call$emmvar
    m2emm = m2$call$emmvar
    mf1 <- paste(paste(formula(m1))[c(2, 1, 3)], collapse = " ")
    mf2 <- paste(paste(formula(m2))[c(2, 1, 3)], collapse = " ")
    if (!any(is.null(m1$expnms))) mf1 = paste0(mf1, ", expnms: ", paste(m1$expnms, collapse=", "))
    if (!any(is.null(m2$expnms))) mf2 = paste0(mf2, ", expnms: ", paste(m2$expnms, collapse=", "))
    if (!any(is.null(m1emm))) mf1 = paste0(mf1, ", EMM: ", m1emm)
    if (!any(is.null(m2emm))) mf2 = paste0(mf2, ", EMM: ", m2emm)
    
    
    mm <- cbind(mm2, mm1)
    qmm <- qr(mm)
    qmmq <- qr.Q(qmm)
    nymm1 <- as.data.frame(qmmq[, 1:qmm$rank])
    colnames(nymm1) <- paste("parm", 1:ncol(nymm1), sep = ".")
    nymm2 <- nymm1[, 1:ncol(mm2), drop = FALSE]
    formula1 <- formula(paste(formula(m1)[[2]], formula(m1)[[1]],
                              paste(c("-1", colnames(nymm1)), collapse = "+"),
                              collapse = ""))
    beta = coef(m1$fit)
    vbeta = vcov(m1$fit)
    df <- dim(mm1)[2] - dim(mm2)[2]
    rbeta <- rep(1, length(beta))
    rbeta[1:df] <- 0
    beta0 <- rev(rbeta)
    zeroidx <- beta0 == 0
    V0 <- vbeta[zeroidx, zeroidx, drop = FALSE]
    b0 <- beta[zeroidx]
    #X2 <- as.numeric(t(b0) %*% ginv(V0) %*% b0) # MASS::ginv is not in the dependencies, reverting to solve
    X2 <- as.numeric(t(b0) %*% solve(V0) %*% b0)
    ev <- eigen(V0, only.values = TRUE)$values
    df.real <- sum(ev > 1e-12)
    topnote <- paste("Model 1", mf1, "\nModel 2", mf2)
    title <- "Analysis of 'Wald statistic' Table\n"
    table <- data.frame(Df = df.real, X2 = X2, p = 1 - pchisq(X2,
                                                              df.real))
    dimnames(table) <- list("1", c("Df", "X2", "P(>|Chi|)"))
    val <- structure(table, heading = c(title, topnote),
                     class = c("anova", "data.frame"))
    return(val)
  }
}

confint.qgcompfit <- function(object, parm, level=0.95, ...){
  #' @exportS3Method stats::confint
  est = coef(object)
  se = sqrt(diag(vcov(object)))
  cbind(lower = est+se*qnorm((1-level)/2), upper = est+se*qnorm(1-(1-level)/2))
}

confint.eefit <- function(object, parm, level=0.95, ...){
  #' @exportS3Method stats::confint
  est = coef(object)
  std = sqrt(diag(vcov(object)))
  critval = qnorm(c((1-level)/2, 1-(1-level)/2))
  res = cbind(est + critval[1]*std, est + critval[2]*std)
  colnames(res) <- paste(100*c((1-level)/2, 1-(1-level)/2), "%")
  res
}


print.qgcompfit <- function(x, showweights=TRUE, ...){
  #' @title Default printing method for a qgcompfit object
  #' 
  #' @description Gives variable output depending on whether `qgcomp.glm.noboot` or `qgcomp.glm.boot`
  #' is called. For `qgcomp.glm.noboot` will output final estimate of joint exposure
  #' effect (similar to the 'index' effect in weighted quantile sums), as well
  #' as estimates of the 'weights' (standardized coefficients). For `qgcomp.glm.boot`,
  #' the marginal effect is given, but no weights are reported since this approach
  #' generally incorporates non-linear models with interaction terms among exposures,
  #' which preclude weights with any useful interpretation.
  #' 
  #' @param x "qgcompfit" object from `qgcomp`, `qgcomp.glm.noboot` or `qgcomp.glm.boot` 
  #' function
  #' @param showweights logical: should weights be printed, if estimated?
  #' @param ... unused
  #' @seealso \code{\link[qgcomp]{qgcomp.glm.noboot}}, \code{\link[qgcomp]{qgcomp.glm.boot}}, and \code{\link[qgcomp]{qgcomp}}
  #' @concept variance mixtures
  #' @exportS3Method base::print
  #' @examples
  #' set.seed(50)
  #' dat <- data.frame(y=runif(50), x1=runif(50), x2=runif(50), z=runif(50))
  #' obj1 <- qgcomp.glm.noboot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2)
  #' obj2 <- qgcomp.glm.boot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, B=10, seed=125)
  #' # does not need to be explicitly called, but included here for clarity
  #' print(obj1)
  #' print(obj2)
  fam <- x$fit$family$family
  rnm =  c(paste0('psi',seq_len(max(1, length(coef(x))))))
  if(x$hasintercept){
    rnm = c("(Intercept)",rnm[-length(rnm)])
  }
  if(inherits(x, "ziqgcompfit")){
    printZI(x, showweights=showweights, ...)
    return(invisible(x))
  }
  if(!is.null(x$pos.size) & showweights) {
    cat(paste0("Scaled effect size (positive direction, sum of positive coefficients = ", signif(x$pos.size, 3) , ")\n"))
    if (length(x$pos.weights) > 0) {
      print(x$pos.weights, digits = 3)
    } else cat("None\n")
    cat("\n")
  }
  if(!is.null(x$neg.size) & showweights) {
    cat(paste0("Scaled effect size (negative direction, sum of negative coefficients = ", signif(-x$neg.size, 3) , ")\n"))
    if (length(x$neg.weights) > 0) {
      print(x$neg.weights, digits = 3)
    } else cat("None\n")
    cat("\n")
  }
  if (fam %in% c("binomial", "quasibinomial")){
    estimand <- 'OR'
    if(x$bootstrap && x$msmfit$family$link=='log') estimand = 'RR'
    cat(paste0("Mixture log(",estimand,")", ifelse(x$bootstrap, " (bootstrap CI)", ifelse(is.null(x$covmat.all_robust), " (delta method CI)", " (robust CI)")), ":\n\n"))
    testtype = "Z"
  }
  if (fam %in% c("poisson", "quasipoisson")){
    #message("Poisson family still experimental: use with caution")
    estimand <- 'RR'
    cat(paste0("Mixture log(",estimand,")", ifelse(x$bootstrap, " (bootstrap CI)", ifelse(is.null(x$covmat.all_robust), " (delta method CI)", " (robust CI)")), ":\n\n"))
    testtype = "Z"
  }
  if (fam %in% c("gaussian", "inverse.gaussian", "tobit")){
    cat(paste0("Mixture slope parameters", ifelse(x$bootstrap, " (bootstrap CI)", ifelse(is.null(x$covmat.all_robust), " (delta method CI)", " (robust CI)")), ":\n\n"))
    testtype = "t"
    x$zstat = x$tstat
  }
  if (fam %in% c("cox", "cch")){
    cat(paste0("Mixture log(hazard ratio)", ifelse(x$bootstrap, " (bootstrap CI)", ifelse(is.null(x$covmat.all_robust), " (delta method CI)", " (robust CI)")), ":\n\n"))
    testtype = "Z"
  }
  if (!(fam %in% c("poisson", "quasipoisson", "binomial", "cox", "cch", "gaussian", "tobit"))){
    testtype = "Z"
    rnm = c(paste0('psi',seq_len(max(1, length(coef(x))))))
    warning(paste0("The ", fam, " distribution has not been tested with qgcomp! Please use with extreme caution
                   and check results thoroughly with simulated data to ensure it works."))
  }
  plab = ifelse(testtype=="Z", "Pr(>|z|)", "Pr(>|t|)")
  if(is.null(dim(x$ci.coef))){
    pdat <- cbind(Estimate=coef(x), "Std. Error"=sqrt(x$var.coef), "Lower CI"=x$ci.coef[1], "Upper CI"=x$ci.coef[2], "test"=x$zstat, "pval"=x$pval)
  } else{
    pdat <- cbind(Estimate=coef(x), "Std. Error"=sqrt(x$var.coef), "Lower CI"=x$ci.coef[,1], "Upper CI"=x$ci.coef[,2], "test"=x$zstat, "pval"=x$pval)
  }
  colnames(pdat)[which(colnames(pdat)=="test")] = eval(paste(testtype, "value"))
  colnames(pdat)[which(colnames(pdat)=="pval")] = eval(paste(plab))
  rownames(pdat) <- rnm
  printCoefmat(pdat,has.Pvalue=TRUE,tst.ind=5L,signif.stars=FALSE, cs.ind=1L:2)
  invisible(x)
}

print.eefit <- function(x, ...){
  #' @exportS3Method base::print
  ci = confint(x, ...)
  zp = ztest.eefit(x, ...)
  pdat <- cbind(Estimate=coef(x), "Std. Error"=sqrt(diag(vcov(x))), "Lower CI"=ci[,1], "Upper CI"=ci[,2], "test"=zp[,1], "pval"=zp[,2])
  testtype = "Z"
  plab = ifelse(testtype=="Z", "Pr(>|z|)", "Pr(>|t|)")
  colnames(pdat)[which(colnames(pdat)=="test")] = eval(paste(testtype, "value"))
  colnames(pdat)[which(colnames(pdat)=="pval")] = eval(paste(plab))
  rownames(pdat) <- names(coef(x))
  printCoefmat(pdat,has.Pvalue=TRUE,tst.ind=5L,signif.stars=FALSE, cs.ind=1L:2)
  invisible(x)
}



summary.qgcompfit <- function(object, ...){
  #' @exportS3Method base::summary
  fam <- object$fit$family$family
  if(inherits(object, "ziqgcompfit")){
    res = summaryZI(object)
    return(res)
  }
  if (fam %in% c("binomial", "quasibinomial")){
    estimand <- 'OR'
    if(object$bootstrap && object$msmfit$family$link=='log') estimand = 'RR'
    cat(paste0("Mixture log(",estimand,")", ifelse(object$bootstrap, " (bootstrap CI)", ifelse(is.null(object$covmat.all_robust), " (delta method CI)", " (robust CI)")), ":\n\n"))
    testtype = "Z"
    rnm = c("(Intercept)", c(paste0('psi',seq_len(max(1, length(coef(object))-1)))))
  }
  if (fam %in% c("poisson", "quasipoisson")){
    #message("Poisson family still experimental: use with caution")
    estimand <- 'RR'
    cat(paste0("Mixture log(",estimand,")", ifelse(object$bootstrap, " (bootstrap CI)", ifelse(is.null(object$covmat.all_robust), " (delta method CI)", " (robust CI)")), ":\n\n"))
    testtype = "Z"
    rnm = c("(Intercept)", c(paste0('psi',seq_len(max(1, length(coef(object))-1)))))
  }
  if (fam %in% c("gaussian", "inverse.gaussian", "tobit")){
    cat(paste0("Mixture slope parameters", ifelse(object$bootstrap, " (bootstrap CI)", ifelse(is.null(object$covmat.all_robust), " (delta method CI)", " (robust CI)")), ":\n\n"))
    testtype = "t"
    rnm = c("(Intercept)", c(paste0('psi',seq_len(max(1, length(coef(object))-1)))))
  }
  if (fam %in% c("cox", "cch")){
    cat(paste0("Mixture log(hazard ratio)", ifelse(object$bootstrap, " (bootstrap CI)", ifelse(is.null(object$covmat.all_robust), " (delta method CI)", " (robust CI)")), ":\n\n"))
    testtype = "Z"
    rnm = c(paste0('psi',seq_len(max(1, length(coef(object))))))
  }
  if (!(fam %in% c("poisson", "quasipoisson", "binomial", "cox", "cch", "gaussian", "tobit"))){
    # in development: quasipoisson, Gamma, quasi, quasibinomial, inverse.gaussian
    testtype = "Z"
    rnm = c(paste0('psi',seq_len(max(1, length(coef(object))))))
    warning(paste0("The ", fam, " distribution has not been tested with qgcomp! Please use with extreme caution
                   and check results thoroughly with simulated data to ensure it works."))
  }
  plab = ifelse(testtype=="Z", "Pr(>|z|)", "Pr(>|t|)")
  if(is.null(dim(object$ci.coef))){
    pdat <- cbind(Estimate=coef(object), "Std. Error"=sqrt(object$var.coef), "Lower CI"=object$ci.coef[1], "Upper CI"=object$ci.coef[2], "test"=object$zstat, "pval"=object$pval)
  } else{
    pdat <- cbind(Estimate=coef(object), "Std. Error"=sqrt(object$var.coef), "Lower CI"=object$ci.coef[,1], "Upper CI"=object$ci.coef[,2], "test"=object$zstat, "pval"=object$pval)
  }
  colnames(pdat)[which(colnames(pdat)=="test")] = eval(paste(testtype, "value"))
  colnames(pdat)[which(colnames(pdat)=="pval")] = eval(paste(plab))
  rownames(pdat) <- rnm
  list(coefficients=pdat)
}

family.qgcompfit <- function(object, ...){
  #' @exportS3Method stats::family
  object$fit$family
}



predict.qgcompfit <- function(object, expnms=NULL, newdata=NULL, type="response", ...){
  #' @title Default prediction method for a qgcompfit object (non-survival 
  #' outcomes only)
  #'
  #' @description get predicted values from a qgcompfit object, or make predictions
  #' in a new set of data based on the qgcompfit object. Note that when making predictions
  #' from an object from qgcomp.glm.boot, the predictions are made from the (conditional) g-computation
  #' model rather than the marginal structural model. Predictions from the marginal
  #' structural model can be obtained via \code{\link[qgcomp]{msm.predict}}. Note
  #' that this function accepts non-quantized exposures in "newdata" and automatically
  #' quantizes them according to the quantile cutpoints in the original fit.
  #' 
  #' @param object "qgcompfit" object from `qgcomp.glm.noboot`, `qgcomp.glm.boot`, `qgcomp.zi.noboot`, 
  #' or `qgcomp.zi.boot`functions
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
  #' @exportS3Method stats::predict
  #' @examples
  #' set.seed(50)
  #' dat <- data.frame(y=runif(50), x1=runif(50), x2=runif(50), z=runif(50))
  #' obj1 <- qgcomp.glm.noboot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2)
  #' obj2 <- qgcomp.glm.boot(y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, B=10, seed=125)
  #' set.seed(52)
  #' dat2 <- data.frame(y=runif(50), x1=runif(50), x2=runif(50), z=runif(50))
  #' summary(predict(obj1, expnms = c('x1', 'x2'), newdata=dat2))
  #' summary(predict(obj2, expnms = c('x1', 'x2'), newdata=dat2))
  if(is.null(newdata)){
    pred <- predict(object$fit, type=type, ...) 
  }
  if(!is.null(newdata)){
    if(is.null(expnms[1])) expnms = object$expnms # testing
    newqdata <- quantize(newdata, expnms, q=NULL, object$breaks)$data
    pred <- predict(object$fit, newdata=newqdata, type=type, ...) 
  }
  return(pred)
}

predict.eefit <- function(object, expnms=NULL, newdata=NULL, type=c("response", "link"), ...){
  #' @exportS3Method stats::predict
  if(is.null(object$X)){
    stop("IncludeX must be TRUE in the call to qgcomp.glm.ee")
  }
  if(type[1]=="response"){
    afun = object$family$linkinv
  } else{
    afun = function(a) a
  }
  if(is.null(newdata)){
    pred <- afun(object$X %*% object$est)
  }
  if(!is.null(newdata)){
    if(is.null(expnms[1])) expnms = object$expnms # testing
    newqdata <- quantize(newdata, expnms, q=NULL, object$breaks)$data
    X = model.matrix(object$formula, newqdata)
    pred <- afun(object$X %*% object$est)
  }
  return(pred)
}



ztest.eefit <- function(object, parm, level=0.95, ...){
  est = coef(object)
  std = sqrt(diag(vcov(object)))
  critval = qnorm(1-(1-level)/2)
  z = (est/std)
  p = 2*(1-pnorm(abs(z), lower.tail = TRUE))
  res = cbind(z,p)
  colnames(res) <- c("Z", "p(>Z)")
  res
}





