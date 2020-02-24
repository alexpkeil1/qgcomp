# experimental functions that may make it into permanent base files

# allow dependencies that don't require installation at initial installation of qgcomp
.qgc.require <- function (package, message = paste("loading required package (", 
                                   package, ") failed", sep = "")){
  if (!requireNamespace(package, quietly = FALSE)) {
    stop(message, call. = FALSE)
  }
  invisible(TRUE)
}


df.residual.qgcompfit <- function(object, ...){
  #' @export
  object$fit$df.residual
}


glance.qgcompfit <- function(x){
  #' @export
  #' @importFrom generics glance
  # @importFrom broom unrowname finish_glance glance
  .qgc.require("broom")
  s <- summary(x$fit)
  # taken from broom package
  unrowname <- function (x) {
    rownames(x) <- NULL
    x
  }
  ret <- unrowname(as.data.frame(s[c("null.deviance", "df.null", "df.residual")]))
  broom::finish_glance(ret, x)
}


tidy.qgcompfit <- function (x, conf.int = FALSE, conf.level = 1-x$alpha, exponentiate = FALSE, 
                            quick = FALSE, ...) {
  #' @export
  #' @importFrom generics tidy
  #' @importFrom tibble as_tibble
  co <- coef(x)
  nms = c("(Intercept)", paste0("psi", 1:(length(names(co))-1)))
  ret <- data.frame(term = nms, estimate = unname(co), 
                    stringsAsFactors = FALSE)
  ll = x$coef + qnorm((1-conf.level)/2)*sqrt(x$var.coef)
  ul = x$coef + qnorm(1-(1-conf.level)/2)*sqrt(x$var.coef)
  z = x$coef/sqrt(x$var.coef)
  ret <- within(ret, {
    "std.error"=sqrt(x$var.coef)
    "Lower CI"=ll
    "Upper CI"=ul
    "test"= z
    "Pr(>|z|)"= x$pval
  }) 
  tibble::as_tibble(ret)
  ret
}


mice.impute.leftcenslognorm <- function(y, ry, x, wy = NULL, lod = NULL, debug=FALSE, ...){
  #' @title imputation for limits of detection problems
  #' 
  #' @description This function integrates with \code{\link[mice]{mice}} to impute values below the LOD using a left
  #' censored log-normal distribution. 
  #' 
  #' @details While this function has utility far beyond qgcomp,
  #' it is included in the qgcomp package because it will be useful for a variety of
  #' settings in which qgcomp is useful. Note that LOD problems where the LOD is small,
  #' and the \code{q} parameter from \code{\link[qgcomp]{qgcomp.noboot}} or 
  #' \code{\link[qgcomp]{qgcomp.boot}} is not large, the LOD may be below the lowest
  #' quantile cutpoint which will yield identical datasets from the MICE procedure in terms
  #' of quantized exposure data. If only exposures are missing, and they have low LODs, then
  #' there will be no benefit in qgcomp from using MICE rather than imputing some small value 
  #' below the LOD.
  #' @param y Vector to be imputed
  #' @param ry Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted. The ry generally distinguishes the observed (TRUE) and missing values (FALSE) in y.
  #' @param x Numeric design matrix with length(y) rows with predictors for y. Matrix x may have no missing values.
  #' @param wy Logical vector of length length(y). A TRUE value indicates locations in y for which imputations are created.
  #' @param lod numeric vector of limits of detection
  #' @param debug logical, print extra info
  #' @param ...	arguments to survreg
  #' @return Vector with imputed data, same type as y, and of length sum(wy)
  #' @import survival
  #' @export
  #' 
  #' @examples
  #' N = 100
  #' set.seed(123)
  #' dat <- data.frame(y=runif(N), x1=runif(N), x2=runif(N), z=runif(N))
  #' true = qgcomp.noboot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, family=gaussian())
  #' mdat <- dat
  #' mdat$x1 = ifelse(mdat$x1>0.5, mdat$x1, NA)
  #' mdat$x2 = ifelse(mdat$x2>0.75, mdat$x2, NA)
  #' cc <- qgcomp.noboot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=mdat[complete.cases(mdat),], q=2, family=gaussian())
  #' 
  #' \donttest{
  #' library("mice")
  #' library("survival")
  #' impdat = mice(data = mdat, 
  #'   method = c("", "leftcenslognorm", "leftcenslognorm", ""),
  #'   lod=c(0.5, 0.75), debug=FALSE)
  #' qc.fit.imp <- list(
  #'   call = call("qgcomp.noboot(y~., expnms = c('x1', 'x2'), family=gaussian())"),
  #'   call1 = impdat$call,
  #'   nmis = impdat$nmis,
  #'   analyses = lapply(1:5, function(x) qgcomp.noboot(y~., expnms = c("x1", "x2"),
  #'     data=complete(impdat, x), family=gaussian(), bayes=TRUE))
  #' )
  #' obj = pool(as.mira(qc.fit.imp))
  #' # true values
  #' true
  #' # complete case analysis
  #' cc
  #' # MI based analysis
  #' summary(obj)
  #' }
  .qgc.require("mvtnorm")
  whichvar = eval(as.name("j"), parent.frame(n = 2))
  nms = names(eval(as.name("data"), parent.frame(n = 3)))
  j = which(nms == whichvar)
  # start
  N = length(y)
  if(is.null(lod)) {
    lod = min(y[ry])
  } else{
    lod = lod[j]
  }
  
  if (is.null(wy)) 
    wy <- !ry
  nmiss = sum(wy)
  x <- cbind(1, as.matrix(x))
  LOD = rep(lod, N)
  ylod = y
  ylod[wy] = LOD[wy]
  fit <-  survreg(
    Surv(time=ylod, event=ry, type='left') ~ x[,-1, drop=FALSE],
    dist='lognormal', ...
  )
  # take a draw from the coefficients
  draw = mvtnorm::rmvnorm(1, c(fit$coefficients, `Log(Scale)`=log(fit$scale)), fit$var)
  #  get predictions under new draws
  fit2 <-  survreg(
    Surv(time=ylod, event=ry, type='left') ~ x[,-1, drop=FALSE],
    dist='lognormal', init=draw, control=survreg.control(maxiter=0)
  )
  fit2$linear.predictors[wy] = ifelse(is.na(fit2$linear.predictors[wy]), -20, fit2$linear.predictors[wy])
  fub = plnorm(LOD[wy], fit2$linear.predictors[wy], fit2$scale)
  # note plnorm(mu,0,1) = pnorm(log(mu), 0, 1)
  u <- runif(nmiss)*fub # random draw from cdf
  #convert cdf back to normal (need to fix with highly negative predictions)
  returny <- qlnorm(u, fit2$linear.predictors[wy], fit2$scale)
  if(any(is.na(returny))){
    ck = data.frame(returny, u,  fit2$linear.predictors[wy], scale = rep(fit2$scale, length(u)))
    print(ck[which(is.na(returny)),])
  }
  if(debug) {
    dlist <- c(
      j = j,
      nmissing = nmiss,
      totalN = N,
      min_imp = min(returny),
      max_imp = max(returny),
      lod = LOD[1]
    )
    cat("\n")
    print(dlist)
  }
  returny
}
