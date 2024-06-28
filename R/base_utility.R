# utility functions that are not part of the base functions but serve miscellaneous purposes


# allow dependencies that don't require installation at initial installation of qgcomp
.qgc.require <- function (package, message = paste("loading required package (",
                                                   package, ") failed", sep = "")){
  if (!requireNamespace(package, quietly = FALSE)) {
    stop(message, call. = FALSE)
  }
  invisible(TRUE)
}

.rmvnorm <- function(n,mu,Sigma){
  # draw from multivariate normal distribution - this is a thin
  #  wrapper for either mvtnorm::rmvnorm or MASS::mvrnorm,
  #  depending on the moods of those package developers
  #.qgc.require("mvtnorm")
  #draw = mvtnorm::rmvnorm(n, mu, Sigma)
  .qgc.require("MASS")
  draw = MASS::mvrnorm(n, mu, Sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  draw
}

construction <- function(i="", j=""){
  msg ="This function is not yet fully implemented"
  if(j != "") msg = paste0(msg, ": ", j)
  if(i %in% c("msg", "message", "warning", "wrn", "warn")) warning(msg)
  else stop(msg)
}

# fake cox family function
cox <- function(){
  obj = binomial(link="log")
  obj$family="cox"
  obj
}

# fake zi family function
zi <- function(){
  obj = binomial(link="log")
  obj$family="zi"
  obj
}

# fake multinomial family function
multinom_family <- function(){
  obj = binomial(link="log")
  obj$family="multinomial"
  obj
}


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
  #' grad_poly procedure, in the case of coming up with proper weights when the combination
  #' of expnms derives from a polynomial function (as in qgcomp.glm.boot with degree>1).
  #'
  #' @export
  #' @examples
  #' vcov = rbind(c(1.2, .9),c(.9, 2.0))
  #' colnames(vcov) <- rownames(vcov) <- expnms <- c("x1", "x2")
  #' se_comb(expnms, vcov, c(1, 0))^2 # returns the given variance
  #' se_comb(expnms, vcov, c(1, 1)) # default linear MSM fit: all exposures
  #' # have equal weight
  #' se_comb(expnms, vcov, c(.3, .1)) # used when one exposure contributes
  #'   # to the overall fit more than others  = d(msmeffect)/dx
  
  if(!is.matrix(covmat)) {
    nm <- names(covmat)
    covmat = matrix(covmat)
    colnames(covmat) <- nm
  }
  weightvec <- rep(0, dim(covmat)[1])
  #if(is.null(grad)) weightvec[which(colnames(as.matrix(covmat)) %in% expnms)] <- 1
  #if(!is.null(grad)) weightvec[which(colnames(as.matrix(covmat)) %in% expnms)] <- grad
  if (is.null(grad)){
    weightvec[which(colnames(as.matrix(covmat)) %in% expnms)] <- 1
  } else if (!is.null(grad) && length(grad)==1){
    weightvec[which(colnames(as.matrix(covmat)) %in% expnms)] <- grad
  } else if (!is.null(grad) && length(grad)==length(weightvec)){
    weightvec <- grad
  }
  var <- weightvec %*% covmat %*% weightvec
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
  #' @export
  #' @examples
  #' vcov = rbind(c(0.010051348, -0.0039332248, -0.0036965571),
  #'              c(-0.003933225,  0.0051807876,  0.0007706792),
  #'              c(-0.003696557,  0.0007706792,  0.0050996587))
  #' colnames(vcov) <- rownames(vcov) <- c("(Intercept)", "x1", "x2")
  #' expnms <- rownames(vcov)[2:3]
  #' aname = rownames(vcov)[1]
  #' vc_comb(aname, expnms, vcov) # returns the given covariance matrix
  
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
  #if(!is.null(grad)) grad = NULL # not yet used
  #if(is.null(grad)) weightvec[which(colnames(as.matrix(covmat)) %in% expnms)] <- 1
  if (is.null(grad)){
    weightvec[which(colnames(as.matrix(covmat)) %in% expnms)] <- 1
  } else if (!is.null(grad) && length(grad)==1){
    weightvec[which(colnames(as.matrix(covmat)) %in% expnms)] <- grad # need more testing with non-unity weights
  } else if (!is.null(grad) && length(grad)==length(weightvec)){
    weightvec <- grad
  }
  outcov = matrix(NA, nrow=2, ncol=2)
  acol = which(colnames(as.matrix(covmat)) %in% aname)
  bcol = which(colnames(as.matrix(covmat)) %in% expnms)
  outcov[1,1] <- covmat[acol,acol]
  outcov[1,2] <- outcov[2,1] <- sum(covmat[acol, bcol])
  outcov[2,2] <- weightvec %*% covmat %*% weightvec # delta method
  outcov
}


grad_poly <- function(intvals, degree){
  #' @import stats
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
  #'  \code{\link[qgcomp]{qgcomp.glm.boot}}, and \code{\link[qgcomp]{qgcomp.glm.noboot}},
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

.qgcomp_object <- function(...){
  res = list(...)
  nms = names(res)
  if(is.na(match("hasintercept", nms))) res$hasintercept = TRUE
  if(is.na(match("bootstrap", nms))) res$bootstrap=FALSE
  if(is.na(match("cov.yhat", nms))) res$cov.yhat=NULL
  if(is.na(match("degree", nms))) res$degree=1
  if(is.na(match("pos.psi", nms))) res$pos.psi = NULL
  if(is.na(match("neg.psi", nms))) res$neg.psi = NULL
  if(is.na(match("pos.weights", nms))) res$pos.weights = NULL
  if(is.na(match("neg.weights", nms))) res$neg.weights = NULL
  if(is.na(match("pos.size", nms))) res$pos.size = NULL
  if(is.na(match("neg.size", nms))) res$neg.size = NULL
  if(is.na(match("df", nms))) res$df = NULL
  if(is.na(match("qx", nms))) res$qx = NULL
  attr(res, "class") <- c("qgcompfit", "list")
  res
}

.qgcompmult_object <- function(...){
  res = .qgcomp_object(...)
  attr(res, "class") <- c("qgcompmultfit", attr(res, "class"))
  res
}
