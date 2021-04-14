# awst

<!-- badges: start -->
[![R-CMD-check-bioc](https://github.com/drisso/awst/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/drisso/awst/actions)
<!-- badges: end -->

Asymmetric Winsorization per Sample Transformation

This package implements the AWST method presented in [Risso and Pagnotta (2020)](https://doi.org/10.1101/2020.06.04.134916).

## Installation

The easiest way to install the package is the following.

```{r}
install.packages("remotes")
remotes::install_github("drisso/awst")
```

## Usage

The AWST transformation can be applied either directly on raw counts, or (recommended) to normalized counts. The input should be a matrix with genes in rows and samples in columns.

```{r}
library(awst)
xt <- awst(x)
```

For detailed examples and use cases see https://github.com/drisso/awst_analysis.
