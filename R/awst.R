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
#' @return a matrix of transformed values.
#'
#' @examples
#' x <- matrix(data = rpois(100, lambda=5), ncol=10, nrow=10)
#' awst(x)
#'
#' @export
awst <- function(x, poscount = FALSE, full_quantile = FALSE, sigma0 = 0.075,
                 lambda = 5) {
  zcount <- score(x, poscount = poscount, full_quantile = full_quantile)

  retval <- ssmooth(zcount, sigma0 = sigma0, lambda = lambda)

  return(retval)
}

#' @importFrom stats approxfun cov density dnorm fivenum integrate pnorm qnorm
#'   quantile sd var
score <- function(x, poscount = FALSE, full_quantile = FALSE) {
  nr <- nrow(x)
  nc <- ncol(x)

  pp <- c(0.001, 0.005, 0.01, 0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95,
          0.975, 0.99, 0.999)
  qq <- qnorm(pp)
  percentiles <- matrix(0, ncol = length(pp), nrow = nr)
  colnames(percentiles) <- paste0(round(100*pp, 1), "%")
  rownames(percentiles) <- rownames(x)

  stats <- cbind(matrix(0, ncol = 3, nrow = nr), percentiles)
  colnames(stats)[1:3] <- c("totalSum", "location", "scale")

  for(i in 1:nr) {

    wd <- x[i,]
    if(poscount) {
      wd <- log(1 + wd[wd > 0])
    } else {
      wd <- log(1 + wd)
    }

    d <- density(wd, n = 1000)
    ccenter <- (d$x[d$x > 1])[which.max(d$y[d$x > 1])]
    wd <- wd - ccenter
    tmp <- wd[wd > 0]
    tmp <- c(tmp, -tmp)
    sigma2 <- var(tmp)
    tau <- sqrt((exp(sigma2) - 1) * exp(2*ccenter + sigma2))
    stats[i, 1:3] <- floor(c(sum(x[i,]), exp(ccenter), tau))
    x[i,] <- (x[i,] - exp(ccenter))/tau

    pcts <- exp(ccenter + sd(tmp) * qq)
    stats[i, 4:ncol(stats)] <- floor(pcts)
    percentiles[i,] <- (pcts - exp(ccenter))/tau
  }

  attr(x, "stats") <- stats
  attr(x, "percentiles") <- percentiles
  return(x)
}

ssmooth <- function(zcount, sigma0 = 0.075, lambda = 5) {
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
  (fivenum(ans))
  (tmp <- YY[which(xx > 0)[1]])
  ans <- 2 * (ans - tmp)/tmp
  rownames(ans) <- rownames(zcount)

  ### finalizing
  attr(ans, "stats") <- attr(zcount, "stats")
  attr(ans, "percentiles") <- attr(zcount, "percentiles")
  attr(ans, "lambda") <- lambda
  return(ans)
}

## unificare con score()
score_TPM <- function(x) {
  nr <- nrow(x)
  nc <- ncol(x)
  pp <- c(0.001, 0.005, 0.01, 0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95,
          0.975, 0.99, 0.999)
  qq <- qnorm(pp)
  percentiles <- matrix(0, ncol = length(pp), nrow = nr)
  colnames(percentiles) <- paste0(round(100*pp, 1), "%")
  rownames(percentiles) <- rownames(x)
  stats <- cbind(matrix(0, ncol = 3, nrow = nr), percentiles)
  colnames(stats)[1:3] <- c("totalSum", "location", "scale")

  for(i in 1:nr) {
    wd <- log(x[i, which(x[i,] > 0)])
    d <- density(wd[wd>0], n = 1000)

    ccenter <- (d$x)[which.max(d$y)]
    wd <- wd - ccenter
    tmp <- wd[wd > 0]
    tmp <- c(tmp, -tmp)
    ssd <- sqrt((exp(sd(tmp)^2)-1)*exp(2*ccenter))
    stats[i, 1:3] <- floor(c(sum(x[i,]), exp(ccenter), ssd))
    x[i,] <- (x[i,] - exp(ccenter))/ssd

    pcts <- exp(ccenter + sd(tmp) * qq)
    stats[i, 4:ncol(stats)] <- floor(pcts)
    percentiles[i,] <- (pcts - exp(ccenter))/ssd
  }

  attr(x, "stats") <- stats
  attr(x, "percentiles") <- percentiles
  return(x)
}
