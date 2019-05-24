#' @export
gene_filter <- function(exprData, nBins = 20, shannon_cutoff = 0.1) {

  bbreaks <- seq(from = min(exprData), to = max(exprData), length.out = nBins + 1)
  halfDelta <- (bbreaks[2] - bbreaks[1])/2
  ddata_cut <- apply(exprData, 1, cut, breaks = bbreaks, include.lowest = TRUE)
  rownames(ddata_cut) <- colnames(exprData)

  tt <- unique(as.vector(ddata_cut))
  ttable <- matrix(0, nrow = ncol(exprData), ncol = length(tt))
  colnames(ttable) <- tt
  rownames(ttable) <- colnames(exprData)

  for(k in seq_len(nrow(ddata_cut))) {
    tt <- table(ddata_cut[k,])
    ttable[k, names(tt)] <- tt
  }

  tmp <- apply(ttable, 1, shannon.tt)
  which_genes <- which(tmp > shannon_cutoff)
  return(exprData[, which_genes])
}

shannon.tt <- function(tt, normalized = TRUE) {
  tt <- tt/sum(tt)
  tt[tt == 0] <- 1
  ans <- -sum(tt * log2(tt))
  if (normalized) {
    tmp <- rep(1/length(tt), length(tt))
    tmp <- -sum(tmp * log2(tmp))
    return(ans/tmp)
  }
  else return(ans)
}
