# awst

<!-- badges: start -->
[![R-CMD-check-bioc](https://github.com/drisso/awst/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/drisso/awst/actions)
<!-- badges: end -->

Asymmetric Winsorization per Sample Transformation

This package implements the AWST method presented in [Risso and Pagnotta (2020)](https://doi.org/10.1101/2020.06.04.134916).

## Installation

In virtually all cases, the package should be installed from Bioconductor, using the following command:

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("awst")
```

In the rare case that you need the development version from Github, you can install it with:

```{r}
BiocManager::install("drisso/awst")
```


## Usage

The AWST transformation can be applied either directly on raw counts, or (recommended) to normalized counts. The input should be either a matrix with genes in rows and samples in columns or a SummarizedExperiment object.

```{r}
library(awst)
xt <- awst(x)
```

For detailed examples and use cases see the vignette and https://github.com/drisso/awst_analysis.
