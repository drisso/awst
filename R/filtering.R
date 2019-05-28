#' @export
gene_filter <- function(exprData, from = min(exprData), to = max(exprData),
                        nBins = 20, heterogeneity_threshold = 0.1) {

  bbreaks <- seq(from = from, to = to, length.out = nBins + 1)
  classes <- names(table(cut(exprData[1, ], breaks = bbreaks, include.lowest = TRUE)))

  ddata_cut <- apply(exprData, 1, cut, breaks = bbreaks, include.lowest = TRUE)
  rownames(ddata_cut) <- colnames(exprData)

  ttable <- matrix(0, nrow = ncol(exprData), ncol = nBins)
  colnames(ttable) <- classes
  rownames(ttable) <- colnames(exprData)

  for(k in seq_len(nrow(ddata_cut))) {
    tt <- table(ddata_cut[k,])
    ttable[k, names(tt)] <- tt
  }

  tmp <- apply(ttable, 1, heterogeneity, nBins = nBins)
  noisy_features <- which(tmp < heterogeneity_threshold)
  ans <- exprData[, -noisy_features]
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
