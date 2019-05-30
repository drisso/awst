#' Gene filtering based on heterogeneity
#'
#' This function filters out genes that show a low heterogeneity, as measured by
#' Shannon's entropy.
#'
#' @details Shannon's entropy is computed on the categorized data after AWST
#'   transformation. Those genes that show a lower entropy than the predefined
#'   threshold are deemed to carry too low information to be useful for the
#'   classification of the samples, and are hence removed.
#'
#' @references Risso and Pagnotta (2019). Within-sample standardization and
#'   asymmetric winsorization lead to accurate classification of RNA-seq
#'   expression profiles. Manuscript in preparation.
#'
#' @param x a matrix of transformed gene expression counts (typically the
#'   results of \code{\link{awst}}).
#' @param from the minimum value from which to start binning data.
#' @param to the maximum value for the binning of the data.
#' @param nBins the number of bins.
#' @param heterogeneity_threshold the trheshold used for the filtering.
#'
#' @return a matrix of transformed values, with genes in columns and samples in
#'   row, ready to be used in distance functions.
#'
#' @examples
#' set.seed(222)
#' x <- matrix(rpois(75, lambda=5), ncol=5, nrow=15)
#' a <- awst(x)
#' gene_filter(a)
#'
#' @export
gene_filter <- function(x, from = min(x), to = max(x),
                        nBins = 20, heterogeneity_threshold = 0.1) {

  bbreaks <- seq(from = from, to = to, length.out = nBins + 1)
  classes <- levels(cut(x[1, ], breaks = bbreaks, include.lowest = TRUE))

  ddata_cut <- apply(x, 1, cut, breaks = bbreaks, include.lowest = TRUE)
  rownames(ddata_cut) <- colnames(x)

  ttable <- matrix(0, nrow = ncol(x), ncol = nBins)
  colnames(ttable) <- classes
  rownames(ttable) <- colnames(x)

  for(k in seq_len(nrow(ddata_cut))) {
    tt <- table(ddata_cut[k,])
    ttable[k, names(tt)] <- tt
  }

  tmp <- apply(ttable, 1, heterogeneity, nBins = nBins)
  noisy_features <- which(tmp < heterogeneity_threshold)
  ans <- x[, -noisy_features]
  attr(ans, "breaks") <- bbreaks
  attr(ans, "noisy_features") <- names(noisy_features)
  return(ans)
}

heterogeneity <- function(empirical_probabilities, nBins = nBins) {
  # Shannon's entropy
  empirical_probabilities <- empirical_probabilities/sum(empirical_probabilities)
  empirical_probabilities[empirical_probabilities == 0] <- 1
  ans <- -sum(empirical_probabilities * log2(empirical_probabilities))
  ans <- ans/log2(nBins)
  return(ans)
}
