context("test-awst")

test_that("function gives consistent results", {
    old_retval <- matrix(data = c(-2, -2, 2.614586, -2, 2.541772,
                                  2.395290, 2.253305, 2.391610, -2, -1.824604,
                                  -2, -2, -1.999993, -1.867436, 2.596351,
                                  -2, 2.662462, -1.999993, -2, -1.824604,
                                  -2, -1.991209, -1.999993, -2, -1.824604),
                         nrow = 5, ncol = 5)

    set.seed(222)
    x <- t(matrix(rpois(25, lambda=5), ncol=5, nrow=5))
    retval <- t(awst(x, lambda = 5))
    expect_true(max(abs(retval - old_retval)) < 1e-5)
})

test_that("filtering gives consistent results", {
    mat <- matrix(data = c(-2, -2, 2.614586, -2, 2.541772,
                            2.395290, 2.253305, 2.391610, -2, -1.824604,
                            -2, -2, -1.999993, -1.867436, 2.596351,
                            -2, 2.662462, -1.999993, -2, -1.824604,
                            -2, -1.991209, -1.999993, -2, -1.824604),
                    nrow = 5, ncol = 5)

    old_retval <- mat[,1:4]

    filtered <- gene_filter(t(mat))
    expect_true(max(abs(t(filtered) - old_retval)) < 1e-5)
})

test_that("full-quantile version gives consistent results", {

    library(EDASeq)
    x <- t(matrix(rpois(25, lambda=5), ncol=5, nrow=5))
    fq <- betweenLaneNormalization(x, which="full")
    retval1 <- awst(fq)
    retval2 <- awst(fq, full_quantile=TRUE)
    expect_equal(retval1, retval2)

})

test_that("SummarizedExperiment and matrix methods give consistent results", {
    library(SummarizedExperiment)
    x <- t(matrix(rpois(25, lambda=5), ncol=5, nrow=5))
    se <- SummarizedExperiment(assays = list(counts = x))

    a <- awst(x)
    se <- awst(se)

    expect_equal(a, assay(se, "awst"))

    b <- gene_filter(a)
    c <- gene_filter(se)

    expect_equal(nrow(b), nrow(c))
})

test_that("SummarizedExperiment method works with non-default assays", {
    x <- t(matrix(rpois(25, lambda=5), ncol=5, nrow=5))
    se <- SummarizedExperiment(assays = list(counts = x, norm_counts = x))

    se <- awst(se)
    se <- awst(se, expr_values = "norm_counts", name = "awst2")

    expect_equal(assay(se, "awst"), assay(se, "awst2"))
})
