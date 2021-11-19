
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GREED : Bayesian greedy clustering

<!-- badges: start -->

[![R build
status](https://github.com/comeetie/greed/workflows/R-CMD-check/badge.svg)](https://github.com/comeetie/greed/actions)
[![CRAN
status](https://www.r-pkg.org/badges/version/greed)](https://CRAN.R-project.org/package=greed)
<!-- badges: end -->

Greed enables model based clustering of networks, matrices of count data
and much more with different types of generative models. Model selection
and clustering is performed in combination by optimizing the Integrated
Classification Likelihood. Details of the algorithms and methods
proposed by this package can be found in Côme, Jouvin, Latouche, and
Bouveyron (2021)
[10.1007/s11634-021-00440-z](https://doi.org/10.1007/s11634-021-00440-z).

<img src="man/figures/greed.png" width=200 align="right" />

The following generative models are available currently :

-   **Stochastic Block Models** (see \``?`Sbm-class\` ),
-   **Degree Corrected Stochastic Block Models** (see
    `` ?`DcSbm-class` ``),
-   **Multinomial Stochastic Block Models** (see
    `` ?`MultSbm-class` ``),
-   **Mixture of Multinomials** (see `` ?`MoR-class` ``),
-   **Latent Class Analysis** (see `` ?`Lca-class` ``),
-   **Gaussian Mixture Model** (see `` ?`Gmm-class` `` and
    `` ?`DiagGmm-class` ``),
-   **Multivariate Mixture of Gaussian Regression Model** (see
    `` ?`MoR-class` ``).

With the Integrated Classification Likelihood, the parameters of the
models are integrated out. This allows a natural regularization for
complex models. Since the Integrated Classification Likelihood penalizes
complex models it allows to automatically find a “natural” value for the
number of clusters *K*<sup>\*</sup>, the user only needs to provide an
initial guess as well as values for the prior parameters (sensible
default values are used if no prior information is available). The
optimization is performed by default thanks to a combination of a greedy
local search and a genetic algorithm. Several optimization algorithms
are available.

Eventually, the whole path of solutions from *K*<sup>\*</sup> to 1
cluster is extracted. This enables a partial ordering of the clusters,
and the evaluation of simpler clustering. The package also provides some
plotting functionality.

## Installation

You can install the development version of greed from
[GitHub](https://github.com/) with:

``` r
#GitHub
install.packages("devtools")
devtools::install_github("comeetie/greed")
```

Or use the CRAN version:

``` r
#CRAN
install.packages("greed")
```

## Usage

The main entry point for using the package is simply the greed function
(`?greed`). The generative model will be chosen automatically to fit
with the data provided, but you may specify another choice with the
model parameter. This is a basic example with the classical Jazz
network:

``` r
library(greed)
data(Jazz)
sol=greed(Jazz)
#> ------- guess DCSBM model fitting ------
#> ################# Generation  1: best solution with an ICL of -28614 and 18 clusters #################
#> ################# Generation  2: best solution with an ICL of -28601 and 17 clusters #################
#> ################# Generation  3: best solution with an ICL of -28588 and 13 clusters #################
#> ################# Generation  4: best solution with an ICL of -28575 and 16 clusters #################
#> ################# Generation  5: best solution with an ICL of -28575 and 16 clusters #################
#> ------- Final clustering -------
#> ICL clustering with a DCSBM model, 15 clusters and an icl of -28557.
```

Here Jazz is a square sparse matrix and a `` ?`DcSbm-class` `` model
will be used by default. Some plotting function enable the exploration
of the clustering results:

``` r
plot(sol)
```

<img src="man/figures/plot-1.png" width="60%" />

And the hierarchical structure between clusters:

``` r
plot(sol,type='tree')
```

<img src="man/figures/tree-1.png" width="60%" />

Eventually, one may explore some coarser clustering using the cut
function:

``` r
plot(cut(sol,5))
```

<img src="man/figures/cut-1.png" width="60%" />

For large datasets, it is possible to use parallelism to speed-up the
computation thanks to the
[future](https://github.com/HenrikBengtsson/future) package. You only
need to specify the type of backend you want to use.

``` r
library(future)
plan(multisession)
data("Blogs")
sol=greed(Blogs$X)
#> ------- guess DCSBM model fitting ------
#> ################# Generation  1: best solution with an ICL of -84503 and 16 clusters #################
#> ################# Generation  2: best solution with an ICL of -84413 and 17 clusters #################
#> ################# Generation  3: best solution with an ICL of -84341 and 19 clusters #################
#> ################# Generation  4: best solution with an ICL of -84293 and 17 clusters #################
#> ################# Generation  5: best solution with an ICL of -84275 and 18 clusters #################
#> ################# Generation  6: best solution with an ICL of -84275 and 18 clusters #################
#> ------- Final clustering -------
#> ICL clustering with a DCSBM model, 17 clusters and an icl of -84250.
plot(sol)
```

<img src="man/figures/future-1.png" width="60%" />
