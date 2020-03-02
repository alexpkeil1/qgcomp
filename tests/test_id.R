cat("# testing use of 'id' to get appropriate standard errors \n")
library("qgcomp")
#install.packages('sandwich', repos="https://cran.rstudio.com/", dependencies = c("Depends"))
2#library("sandwich")
set.seed(2123)
N = 250
t = 4
dat <- data.frame(row.names = 1:(N*t))
dat <- within(dat, {
  id = do.call("c", lapply(1:N, function(x) rep(x, t)))
  u =  do.call("c", lapply(1:N, function(x) rep(runif(1), t)))
  x1 = rnorm(N, u)
  x2 = rnorm(N, u)
  y = rnorm(N) + u + x1 - x2*x1
})


# from sandwich package 2.5.1 (to avoid having to install zoo dependency)
bread.glm <- function (x, ...) {
  if (!is.null(x$na.action)) 
    class(x$na.action) <- "omit"
  sx <- summary(x)
  wres <- as.vector(residuals(x, "working")) * weights(x, "working")
  dispersion <- if (substr(x$family$family, 1, 17) %in% c("poisson", 
                                                          "binomial", "Negative Binomial")) 
    1
  else sum(wres^2)/sum(weights(x, "working"))
  return(sx$cov.unscaled * as.vector(sum(sx$df[1:2])) * dispersion)
}

bread <- function (x, ...) 
{
  UseMethod("bread")
}


sandwich <- function (x, bread. = bread, meat. = meat, ...) {
  if (is.list(x) && !is.null(x$na.action)) 
    class(x$na.action) <- "omit"
  if (is.function(bread.)) 
    bread. <- bread.(x)
  if (is.function(meat.)) 
    meat. <- meat.(x, ...)
  n <- NROW(estfun(x))
  return(1/n * (bread. %*% meat. %*% bread.))
}

meatCL <- function (x, cluster = NULL, type = NULL, cadjust = TRUE, multi0 = FALSE, 
          ...) 
{
  if (is.list(x) && !is.null(x$na.action)) 
    class(x$na.action) <- "omit"
  ef <- estfun(x, ...)
  k <- NCOL(ef)
  n <- NROW(ef)
  rval <- matrix(0, nrow = k, ncol = k, dimnames = list(colnames(ef), 
                                                        colnames(ef)))
  if (is.null(cluster)) 
    cluster <- attr(x, "cluster")
  if (is.null(cluster)) 
    cluster <- 1L:n
  if (inherits(cluster, "formula")) {
    cluster_tmp <- expand.model.frame(x, cluster, na.expand = FALSE)
    cluster <- model.frame(cluster, cluster_tmp, na.action = na.pass)
  }
  else {
    cluster <- as.data.frame(cluster)
  }
  if ((n != NROW(cluster)) && !is.null(x$na.action) && (class(x$na.action) %in% 
                                                        c("exclude", "omit"))) {
    cluster <- cluster[-x$na.action, , drop = FALSE]
  }
  if (NROW(cluster) != n) 
    stop("number of observations in 'cluster' and 'estfun()' do not match")
  p <- NCOL(cluster)
  if (p > 1L) {
    cl <- lapply(1L:p, function(i) combn(1L:p, i, simplify = FALSE))
    cl <- unlist(cl, recursive = FALSE)
    sign <- sapply(cl, function(i) (-1L)^(length(i) + 1L))
    paste_ <- function(...) paste(..., sep = "_")
    for (i in (p + 1L):length(cl)) {
      cluster <- cbind(cluster, Reduce(paste_, unclass(cluster[, 
                                                               cl[[i]]])))
    }
    if (multi0) 
      cluster[[length(cl)]] <- 1L:n
  }
  else {
    cl <- list(1)
    sign <- 1
  }
  g <- sapply(1L:length(cl), function(i) {
    if (is.factor(cluster[[i]])) {
      length(levels(cluster[[i]]))
    }
    else {
      length(unique(cluster[[i]]))
    }
  })
  if (is.null(type)) {
    type <- if (class(x)[1L] == "lm") 
      "HC1"
    else "HC0"
  }
  type <- match.arg(type, c("HC", "HC0", "HC1", "HC2", "HC3"))
  if (type == "HC") 
    type <- "HC0"
  if (type %in% c("HC2", "HC3")) {
    if (any(g == n)) 
      h <- hatvalues(x)
    if (!all(g == n)) {
      if (!(class(x)[1L] %in% c("lm", "glm"))) 
        warning("clustered HC2/HC3 are only applicable to (generalized) linear regression models")
      X <- model.matrix(x)
      if (any(alias <- is.na(coef(x)))) 
        X <- X[, !alias, drop = FALSE]
      attr(X, "assign") <- NULL
      w <- weights(x, "working")
      XX1 <- if (is.null(w)) 
        chol2inv(qr.R(qr(X)))
      else chol2inv(qr.R(qr(X * sqrt(w))))
      res <- rowMeans(ef/X, na.rm = TRUE)
      res[apply(abs(ef) < .Machine$double.eps, 1L, all)] <- 0
      matpower <- function(X, p) {
        if ((ncol(X) == 1L) && (nrow(X) == 1L)) 
          return(X^p)
        Xeig <- eigen(X, symmetric = TRUE)
        if (any(Xeig$values < 0)) 
          stop("matrix is not positive semidefinite")
        sqomega <- diag(Xeig$values^p)
        return(Xeig$vectors %*% sqomega %*% t(Xeig$vectors))
      }
    }
  }
  for (i in 1L:length(cl)) {
    efi <- ef
    adj <- if (multi0 & (i == length(cl))) {
      if (type == "HC1") 
        (n - k)/(n - 1L)
      else 1
    }
    else {
      if (cadjust) 
        g[i]/(g[i] - 1L)
      else 1
    }
    if (type %in% c("HC2", "HC3")) {
      if (g[i] == n) {
        efi <- if (type == "HC2") {
          efi/sqrt(1 - h)
        }
        else {
          efi/(1 - hatvalues(x))
        }
      }
      else {
        for (j in unique(cluster[[i]])) {
          ij <- which(cluster[[i]] == j)
          Hij <- if (is.null(w)) {
            X[ij, , drop = FALSE] %*% XX1 %*% t(X[ij, 
                                                  , drop = FALSE])
          }
          else {
            X[ij, , drop = FALSE] %*% XX1 %*% t(X[ij, 
                                                  , drop = FALSE]) %*% diag(w[ij], nrow = length(ij), 
                                                                            ncol = length(ij))
          }
          Hij <- if (type == "HC2") {
            matpower(diag(length(ij)) - Hij, -0.5)
          }
          else {
            solve(diag(length(ij)) - Hij)
          }
          efi[ij, ] <- drop(Hij %*% res[ij]) * X[ij, 
                                                 , drop = FALSE]
        }
      }
      efi <- sqrt((g[i] - 1L)/g[i]) * efi
    }
    efi <- if (g[i] < n) 
      apply(efi, 2L, rowsum, cluster[[i]])
    else efi
    rval <- rval + sign[i] * adj * crossprod(efi)/n
  }
  if (type == "HC1") 
    rval <- (n - 1L)/(n - k) * rval
  return(rval)
}

estfun.glm <- function (x, ...) {
  xmat <- model.matrix(x)
  xmat <- naresid(x$na.action, xmat)
  if (any(alias <- is.na(coef(x)))) 
    xmat <- xmat[, !alias, drop = FALSE]
  wres <- as.vector(residuals(x, "working")) * weights(x, "working")
  dispersion <- if (substr(x$family$family, 1, 17) %in% c("poisson", 
                                                          "binomial", "Negative Binomial")) 
    1
  else sum(wres^2, na.rm = TRUE)/sum(weights(x, "working"), 
                                     na.rm = TRUE)
  rval <- wres * xmat/dispersion
  attr(rval, "assign") <- NULL
  attr(rval, "contrasts") <- NULL
  #res <- residuals(x, type = "pearson")
  #if (is.ts(res)) 
  #  rval <- ts(rval, start = start(res), frequency = frequency(res))
  #if (is.zoo(res)) 
  #  rval <- zoo(rval, index(res), attr(res, "frequency"))
  return(rval)
}

estfun <- function (x, ...) 
{
  UseMethod("estfun")
}

vcovCL <- function (x, cluster = NULL, type = NULL, sandwich = TRUE, fix = FALSE, 
          ...) 
{
  rval <- meatCL(x, cluster = cluster, type = type, ...)
  if (sandwich) 
    rval <- sandwich(x, meat. = rval)
  if (fix && any((eig <- eigen(rval, symmetric = TRUE))$values < 
                 0)) {
    eig$values <- pmax(eig$values, 0)
    rval[] <- crossprod(sqrt(eig$values) * t(eig$vectors))
  }
  return(rval)
}
################################################################################
#
#
#
################################################################################
# pre quantize
expnms = c("x1")
datl = quantize(dat, expnms = expnms)


#' \donttest{
#' 
#' # delta method/bootstrap variance ignoring clustering
#' noclust = qgcomp.noboot(y~ x1, data=datl$dat, id="id", family=gaussian(), q = NULL)
#' noclust.b = qgcomp.boot(y~ x1, data=datl$dat, family=gaussian(), q = NULL, MCsize=1000)
#' 
#' # bootstrap variance with sampling by cluster
#' clust.b = qgcomp.boot(y~ x1, data=datl$dat, id="id", family=gaussian(), q = NULL, MCsize=5000, B = 500)
#' #clust.g = summary(geeglm(y~x1, data=datl$dat, id=id, corstr = "independence"))
#' fitglm = glm(y~x1, data=datl$dat)
#' # cluster robust variance
#' sw.cov = vcovCL(fitglm, cluster=~id, type = "HC0")[2,2]
#' 
#' stopifnot(all.equal(clust.b$var.psi, sw.cov, tolerance = 0.005))
#' 
#' # change in variance should be the same in gee and bootstrap
#' stopifnot( 
#'   (clust.b$var.psi > noclust.b$var.psi) & (sw.cov > noclust.b$var.psi) |
#'   (clust.b$var.psi < noclust.b$var.psi) & (sw.cov < noclust.b$var.psi) 
#'   )
#' stopifnot( 
#'   (clust.b$var.psi > noclust$var.psi) & (sw.cov > noclust$var.psi) |
#'   (clust.b$var.psi < noclust$var.psi) & (sw.cov < noclust$var.psi) 
#' )
#' }
cat("done")
