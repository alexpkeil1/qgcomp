.devinstall <- function(...){
  .qgc.require("devtools")
  devtools::install_github("alexpkeil1/qgcomp",...)
}

# experimental functions that may make it into permanent base files

.split.iid.data <- function(data, prop.train=0.4){
  trainidx <- sort(sample(1:nrow(data), round(nrow(data)*prop.train)))
  valididx <- setdiff(1:nrow(data),trainidx)
  traindata <- data[trainidx,]
  validdata <- data[valididx,]
  list(
    trainidx = trainidx,
    valididx = valididx,
    traindata = traindata,
    validdata = validdata
  )
}

.split.cluster.data <- function(data, cluster="id", prop.train=0.4){
 ids = sort(unique(data[[cluster]]))
 trididx = sort(sample(1:length(ids), round(length(ids) * 
                                               prop.train)))
 trainids = ids[trididx]
 validids = setdiff(ids, trainids)
 trainidx <- which(data[[cluster]] %in% trainids)
 valididx <- which(data[[cluster]] %in% validids)
 traindata <- data[trainidx, ]
 validdata <- data[valididx, ]
 list(
    trainidx = trainidx,
    valididx = valididx,
    traindata = traindata,
    validdata = validdata
 )
}

qgcomp.partials <- function(
  fun = c("qgcomp.noboot", "qgcomp.cox.noboot", "qgcomp.zi.noboot"),
  traindata=NULL,
  validdata=NULL,
  expnms=NULL,
  ...
){
  #' @title Partial effect sizes, confidence intervals, hypothesis tests
  #' @description Obtain effect estimates for "partial positive" and "partial
  #' negative" effects using quantile g-computation. This approach uses sample
  #' splitting to evaluate the overall impact of a set of variables with 
  #' effects in a single direction, where, using training data, all variables
  #' with effects in the same direction are grouped.  
  #' 
  #' @details In the basic (non bootstrapped) `qgcomp` functions, the positive and 
  #' negative "sums
  #' of coefficients" or "partial effect sizes" are given, which equal the sum
  #' of the negative and positive coefficients in the underlying model. Unfortunately,
  #' these partial effects don't constitute variables for which we can derive confidence
  #' intervals or hypothesis tests, so they are mainly for exploratory purposes. By employing
  #' sample splitting, however, we can obtain better estimates of these partial effects.
  #' 
  #' Sample splitting proceeds by partitioning the data into two samples (40/60 training/validtion
  #' split seems acceptable in many circumstances). The "overall mixture effect" is then
  #' estimated in the training data, and the mixture variables with positive and negative coefficients
  #' are split into separate groups. These two different groups are then used as 
  #' "the mixture of interest" in two additional qgcomp fits, where the mixture of interest
  #' is adjusted for the other exposure variables. For example, if the "positive partial effect"
  #' is of interest, then this effect is equal to the sum of the coefficients in the 
  #' qgcomp model fit to the validation data, with the mixture of interest selected by the
  #' original fit to the training data (note that some of these coefficients may be negative
  #' in the fit to the validation data - this is expected and necessary for valid hypothesis tests).
  #' 
  #' The positive/negative partial effects are necessarily exploratory, but sample splitting preserves
  #' the statistical properties at the expense of wider confidence intervals and larger variances. The
  #' two resulting mixture groups groups should be inspected for 
  #' @param fun character variable in the set "qgcomp.noboot" (binary, count, continuous outcomes), 
  #' "qgcomp.cox.noboot" (survival outcomes), 
  #' "qgcomp.zi.noboot" (zero inflated outcomes). This describes which qgcomp package
  #' function is used to fit the model. (default = "qgcomp.noboot")
  #' @return A 'qgcompmultifit' object, which inherits from \code{\link[base]{list}}, which contains
  #' \describe{
  #' \item{posmix}{character vector of variable names with positive coefficients in the qgcomp model 
  #'   fit to the training data} 
  #' \item{negmix}{character vector of variable names with negative coefficients in the qgcomp model 
  #'   fit to the training data} 
  #' \item{pos.fit}{a qgcompfit object fit to the validation data, in which the exposures of
  #'  interest are contained in 'posmix'} 
  #' \item{neg.fit}{a qgcompfit object fit to the validation data, in which the exposures of
  #'  interest are contained in 'negmix'}
  #'  }
  #' @param traindata Data frame with training data
  #' @param validdata Data frame with validation data
  #' @param expnms Exposure mixture of interest
  #' @param ... Arguments to \code{\link[qgcomp]{qgcomp.noboot}}, 
  #'    \code{\link[qgcomp]{qgcomp.cox.noboot}}, or 
  #'    \code{\link[qgcomp]{qgcomp.zi.noboot}}
  #' @export
  #' @examples 
  #' set.seed(123223)
  #' dat = qgcomp:::.dgm_quantized(N=1000, coef=c(0.25,-0.25,0,0), ncor=1)
  #' cor(dat)
  #' # overall fit (more or less null due to counteracting exposures)
  #' (overall <- qgcomp.noboot(f=y~., q=NULL, expnms=c("x1", "x2", "x3", "x4"), data=dat))
  #' 
  #' # partial effects using 40% training/60% validation split
  #' trainidx <- sample(1:nrow(dat), round(nrow(dat)*0.4))
  #' valididx <- setdiff(1:nrow(dat),trainidx)
  #' traindata = dat[trainidx,]
  #' validdata = dat[valididx,]
  #' splitres <- qgcomp:::qgcomp.partials(fun="qgcomp.noboot", f=y~., q=NULL, 
  #'     traindata=traindata,validdata=validdata, expnms=c("x1", "x2", "x3", "x4"))
  #' splitres
  #' \donttest{
  #' # under the null, both should give null results
  #' set.seed(123223)
  #' dat <- qgcomp:::.dgm_quantized(N=1000, coef=c(0,0,0,0), ncor=1)
  #' # 40% training/60% validation
  #' trainidx2 <- sample(1:nrow(dat), round(nrow(dat)*0.4))
  #' valididx2 <- setdiff(1:nrow(dat),trainidx2)
  #' traindata2 <- dat[trainidx2,]
  #' validdata2 <- dat[valididx2,]
  #' splitres2 <- qgcomp:::qgcomp.partials(fun="qgcomp.noboot", f=y~., 
  #'    q=NULL, traindata=traindata2,validdata=validdata2, expnms=c("x1", "x2", "x3", "x4"))
  #' splitres2
  #' 
  #' # 60% training/40% validation
  #' trainidx3 <- sample(1:nrow(dat), round(nrow(dat)*0.6))
  #' valididx3 <- setdiff(1:nrow(dat),trainidx3)
  #' traindata3 <- dat[trainidx3,]
  #' validdata3 <- dat[valididx3,]
  #' splitres3 <- qgcomp:::qgcomp.partials(fun="qgcomp.noboot", f=y~., q=NULL, 
  #'     traindata=traindata3,validdata=validdata3, expnms=c("x1", "x2", "x3", "x4"))
  #' splitres3
  #' 
  #' # survival outcome
  #' set.seed(50)
  #' N=1000
  #' dat <- data.frame(time=(tmg <- pmin(.1,rweibull(N, 10, 0.1))), 
  #'      d=1.0*(tmg<0.1), x1=runif(N)+(tmg<0.1)*0.1, x2=runif(N)-(tmg<0.1)*0.1, x3=runif(N),
  #'                       x4=runif(N), x5=runif(N) , z=runif(N))
  #' trainidx4 <- sample(1:nrow(dat), round(nrow(dat)*0.6))
  #' valididx4 <- setdiff(1:nrow(dat),trainidx4)
  #' traindata4 <- dat[trainidx4,]
  #' validdata4 <- dat[valididx4,]
  #' expnms=paste0("x", 1:5)
  #' f = survival::Surv(time, d)~x1 + x2 + x3 + x4 + x5 + z
  #' (fit1 <- survival::coxph(f, data = dat))
  #' (overall <- qgcomp.cox.noboot(f, expnms = expnms, data = dat))
  #' (splitres4 <- qgcomp:::qgcomp.partials(fun="qgcomp.cox.noboot", f=f, q=4,
  #'                traindata=traindata4,validdata=validdata4,
  #'                 expnms=expnms))
  #'  
  #'  # zero inflated count outcome
  #'  set.seed(50)
  #'  n=1000
  #'  dat <- data.frame(y= (yany <- rbinom(n, 1, 0.5))*(ycnt <- rpois(n, 1.2)), x1=runif(n)+ycnt*0.2, 
  #'                        x2=runif(n)-ycnt*0.2, x3=runif(n),
  #'                       x4=runif(n) , z=runif(n))
  #'  # poisson count model, mixture in both portions, but note that the qgcomp.partials
  #'  # function defines the "positive" variables only  by the count portion of the model
  #'  (overall5 <- qgcomp.zi.noboot(f=y ~ z + x1 + x2 + x3 + x4 | x1 + x2 + x3 + x4 + z, 
  #'                  expnms = c("x1", "x2", "x3", "x4"), 
  #'                   data=dat, q=4, dist="poisson"))
  #'
  #' trainidx5 <- sample(1:nrow(dat), round(nrow(dat)*0.6))
  #' valididx5 <- setdiff(1:nrow(dat),trainidx5)
  #' traindata5 <- dat[trainidx5,]
  #' validdata5 <- dat[valididx5,]
  #' splitres5 <- qgcomp.partials(fun="qgcomp.zi.noboot", 
  #'     f=y ~ x1 + x2 + x3 + x4 + z | x1 + x2 + x3 + x4 + z, q=4, 
  #'     traindata=traindata5, validdata=validdata5, 
  #'     expnms=c("x1", "x2", "x3", "x4"))
  #' splitres5
  #'                 
  #' }
  if(is.null(traindata) | is.null(validdata))
    stop("traindata and validdata must both be specified")
  whichfun = fun[1]
  train.fit = switch(whichfun, 
                     qgcomp.noboot = qgcomp.noboot(data=traindata, expnms=expnms, ...),
                     qgcomp.cox.noboot = qgcomp.cox.noboot(data=traindata, expnms=expnms, ...),
                     qgcomp.zi.noboot = qgcomp.zi.noboot(data=traindata, expnms=expnms, ...)
  )
  posnms = names(train.fit$pos.weights)
  negnms = names(train.fit$neg.weights)
  if(length(posnms)==1 && all(posnms==c("count", "zero"))){
    posnms = names(train.fit$pos.weights$count)
    negnms = names(train.fit$neg.weights$count)
  }
  res = list(train.fit=train.fit)
  res$negmix <- res$posmix <- "none"
  if(length(posnms)>0){
    res$posmix = posnms
    res$pos.fit <- switch(whichfun, 
                      qgcomp.noboot = qgcomp.noboot(data=validdata, expnms=posnms, ...),
                      qgcomp.cox.noboot = qgcomp.cox.noboot(data=validdata, expnms=posnms, ...),
                      qgcomp.zi.noboot = qgcomp.zi.noboot(data=validdata, expnms=posnms, ...)
    )
  }
  if(length(negnms)>0){
    res$negmix = negnms
    res$neg.fit <- switch(whichfun, 
                      qgcomp.noboot = qgcomp.noboot(data=validdata, expnms=negnms, ...),
                      qgcomp.cox.noboot = qgcomp.cox.noboot(data=validdata, expnms=negnms, ...),
                      qgcomp.zi.noboot = qgcomp.zi.noboot(data=validdata, expnms=negnms, ...)
    )
    
  }
  class(res) <- "qgcompmultifit"
  res
}



print.qgcompmultifit <- function(x, ...){
  #' @export
  cat(paste0("\nVariables with positive effect sizes in training data: ", paste(x$posmix, collapse = ", ")))
  cat(paste0("\nVariables with negative effect sizes in training data: ", paste(x$negmix, collapse = ", ")))
  cat("\nPartial effect sizes estimated in validation data\n")
  cat("Positive direction ")
  if(!is.null(x$pos.fit)){
    print(x$pos.fit, showweights=FALSE)
  } else{
    cat("\n No positive coefficients in model fit to training data ")
  }
  cat("\nNegative direction ")
  if(!is.null(x$neg.fit)){
    print(x$neg.fit, showweights=FALSE)
  } else{
    cat("\n No negative coefficients in model fit to training data ")
  } 
}
