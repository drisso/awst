#' Asymmetric Within-Sample Transformation
#'
#' This function implements the asymmetric within-sample transformation
#' described in Risso and Pagnotta (2019). The function includes two steps: a
#' standardization step and a asymmetric winsorization step. See details.
#'
#' @details The standardization step is based on a log-normal distribution of
#'   the high-intensity genes. Optionally, only positive counts can be used in
#'   this step (this option is especially useful for single-cell data). The
#'   winsorization step is controlled by two parameters, sigma0 and lambda,
#'   which control the growth rate of the winsorization function.
#'
#' @references Risso and Pagnotta (2019). Within-sample standardization and
#'   asymmetric winsorization lead to accurate classification of RNA-seq
#'   expression profiles. Manuscript in preparation.
#'
#' @param x a matrix of (possibly normalized) RNA-seq read counts.
#' @param poscount a logical value indicating whether positive counts only
#'   should be used for the standardization step.
#' @param full_quantile a logical value indicating whether the data have been
#'   normalized with the full-quantile normalization. In this case, computations
#'   can be sped up.
#' @param sigma0 a multiplicative constant to be applied to the smoothing
#'   function.
#' @param lambda a parameter that controls the growth rate of the smoothing
#'   function.
#'
#' @return a matrix of transformed values, with genes in columns and samples in
#'   row, ready to be used in distance functions.
#'
#' @examples
#' x <- matrix(data = rpois(100, lambda=5), ncol=10, nrow=10)
#' awst(x)
#'
#' @export
awst <- function(x, poscount = FALSE, full_quantile = FALSE, sigma0 = 0.075,
                 lambda = 13) {

  zcount <- score(x, poscount = poscount, full_quantile = full_quantile)
  retval <- ssmooth(zcount, sigma0 = sigma0, lambda = lambda)

  return(retval)
}

#' @importFrom stats approxfun cov density dnorm integrate pnorm qnorm
#'   quantile sd var
score <- function(x, poscount = FALSE, full_quantile = FALSE) {

  if(full_quantile) {

    ccenter <- .compute_center(x[,1], poscount = poscount)
    tau <- .compute_tau(x[,1], ccenter, poscount = poscount)
    retval <- t((x - exp(ccenter))/tau)

  } else {

    ccenter <- apply(x, 2, .compute_center, poscount = poscount)
    tau <- lapply(1:NCOL(x), function(i) {
      .compute_tau(x[,i], ccenter[i], poscount = poscount)
    })
    tau <- simplify2array(tau)

    retval <- (t(x) - exp(ccenter))/tau

  }

  return(retval)
}

ssmooth <- function(zcount, sigma0 = 0.075, lambda = 13) {
  ### distribution
  ssigma <- function(z) return(1 + ifelse(z > 0, 2*lambda*(pnorm(z)-.5), 0))
  foo <- function(z, sigma = sigma0)  return(dnorm(z, sd = sigma * ssigma(z)))
  (tmp <- integrate(foo, -Inf, Inf)$value)
  xx <- seq(from = -1, to = 10, length.out = 10000)
  yy <- foo(xx)
  YY <- cumsum(yy[1:(length(xx)-1)]*(xx[2:length(xx)]-xx[1:(length(xx)-1)]))
  YY <- YY/tmp
  YY[YY > 1] <- 1

  ### exprData
  f <- approxfun(xx[-1], YY, yleft = 0, yright = 1)
  ans <- apply(zcount, 2, f)
  (tmp <- YY[which(xx > 0)[1]])
  ans <- 2 * (ans - tmp)/tmp
  rownames(ans) <- rownames(zcount)

  ### finalizing
  attr(ans, "lambda") <- lambda
  return(ans)
}

.compute_center <- function(x, poscount = FALSE, n = 1000) {

  if(poscount) {
    wd <- log1p(x[x>0])
  } else {
    wd <- log1p(x)
  }

  d <- density(wd, n = n)
  ccenter <- (d$x[d$x > 1])[which.max(d$y[d$x > 1])]
  return(ccenter)
}

.compute_tau <- function(x, ccenter, poscount = FALSE) {

  if(poscount) {
    wd <- log1p(x[x>0])
  } else {
    wd <- log1p(x)
  }

  wd <- wd - ccenter
  tmp <- wd[wd > 0]
  tmp <- c(tmp, -tmp)
  sigma2 <- var(tmp)
  tau <- sqrt((exp(sigma2) - 1) * exp(2*ccenter + sigma2))
}
