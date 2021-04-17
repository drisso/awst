#' Gene filtering based on heterogeneity
#'
#' This function filters out genes that show a low heterogeneity, as measured by
#' Shannon's entropy.
#'
#' @details Shannon's entropy is computed on the categorized data after AWST
#'    transformation. Those genes that show a lower entropy than the predefined
#'    threshold are deemed to carry too low information to be useful for the
#'    classification of the samples, and are hence removed.
#'
#' @references Risso and Pagnotta (2019). Within-sample standardization and
#'    asymmetric winsorization lead to accurate classification of RNA-seq
#'    expression profiles. Manuscript in preparation.
#'
#' @param x a matrix of transformed gene expression counts (typically the
#'    results of \code{\link{awst}}).
#' @param from the minimum value from which to start binning data.
#' @param to the maximum value for the binning of the data.
#' @param nBins the number of bins.
#' @param heterogeneity_threshold the trheshold used for the filtering.
#'
#' @return if `x` is a matrix, it returns a filtered matrix. If `x` is a
#'   `SummarizedExperiment`, it returns a filtered `SummarizedExperiment`
#'
#' @examples
#' set.seed(222)
#' x <- matrix(rpois(75, lambda=5), ncol=5, nrow=15)
#' a <- awst(x)
#' gene_filter(a)
#'
#' @export
#' @name gene_filter
NULL

#' @export
setGeneric("gene_filter", function(x, ...) standardGeneric("gene_filter"))

#' @describeIn gene_filter the input is a matrix of awst-transformed values.
#' @export
setMethod("gene_filter", "matrix",
          function(x, from = min(x, na.rm = TRUE),
                   to = max(x, na.rm = TRUE),
                   nBins = 20, heterogeneity_threshold = 0.1) {

              noisy_features <- get_noisy_features(x = x, from = from, to = to,
                                                   nBins = nBins,
                                                   heterogeneity_threshold = heterogeneity_threshold)
              ans <- x[-noisy_features, ]
              return(ans)
})

#' @export
#' @import SummarizedExperiment
#' @describeIn gene_filter the input is a SummarizedExperiment with
#'   awst-transformed values in one of its assays.
#' @param awst_values integer scalar or string indicating the assay that
#'   contains the awst-transformed values to use as input.
setMethod("gene_filter", "SummarizedExperiment",
          function(x, from = min(assay(x, awst_values), na.rm = TRUE),
                   to = max(assay(x, awst_values), na.rm = TRUE),
                   nBins = 20, heterogeneity_threshold = 0.1,
                   awst_values = "awst") {

              noisy_features <- get_noisy_features(x = assay(x, awst_values),
                                                   from = from, to = to,
                                                   nBins = nBins,
                                                   heterogeneity_threshold = heterogeneity_threshold)

              return(x[-noisy_features,])

})

get_noisy_features <- function(x, from, to, nBins, heterogeneity_threshold) {

    bbreaks <- seq(from = from, to = to, length.out = nBins + 1)
    classes <- levels(cut(x[, 1], breaks = bbreaks, include.lowest = TRUE))

    ddata_cut <- apply(x, 2, cut, breaks = bbreaks, include.lowest = TRUE)
    rownames(ddata_cut) <- rownames(x)

    ttable <- matrix(0, nrow = nrow(x), ncol = nBins)
    colnames(ttable) <- classes
    rownames(ttable) <- rownames(x)

    for(k in seq_len(nrow(ddata_cut))) {
        tt <- table(ddata_cut[k,])
        ttable[k, names(tt)] <- tt
    }

    tmp <- apply(ttable, 1, heterogeneity, nBins = nBins)
    return(which(tmp < heterogeneity_threshold))
}

heterogeneity <- function(empirical_probabilities, nBins = nBins) {
    # Shannon's entropy
    empirical_probabilities <- empirical_probabilities/sum(empirical_probabilities)
    empirical_probabilities[empirical_probabilities == 0] <- 1
    ans <- -sum(empirical_probabilities * log2(empirical_probabilities))
    ans <- ans/log2(nBins)
    return(ans)
}
