
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GREED : Bayesian greedy clustering <img src="./inst/exdata/greed.png" width=200 align="right" />

<!-- badges: start -->

<!-- badges: end -->

Greed enable model based clustering of networks, counts data matrix and
much more with different type of generative models. Model selection and
clustering is performed in combination by optimizing the Integrated
Classification Likelihood.

Four generative models are availables currently :

  - sbm : Stochastick Block Models (for directed and undirected graphs),
  - dcsbm : Degree Corrected Stochastick Block Models (for directed and
    undirected graphs),
  - co\_dcsbm: Degree Corrected Latent Block Models,
  - mm: Mixture of Multinomials,
  - gmm : Gaussian Mixture Model.

With the Integrated Classification Likelihood the parameters of the
models are integrated out. Since the Integrated Classification
Likelihood introduces a natural regularisation for complex models such
strategie automaticaly find a “natural” value for the number of cluster,
the user needs only to provide an initial guess. The optimization is
performed by default thanks to a combination of greedy local search and
a genetic algorithm, several optimization algorithms are available.

Eventually, the whole path of solutions from K\* to 1 cluster is
extracted. This enable a partial ordering of the clusters, and the
evaluation of simpler clustering. The package also provides some ploting
functionality.

## Installation

You can install the released version of greed from
[GitHub](https://github.com/) with:

``` r
install.packages("devtools")
devtools::install_github("comeetie/greed")
```

## Example

This is a basic example whith the classical blogs network:

``` r
library(greed)
data(Jazz)
sol=greed(Jazz)
#> ------- Fitting a dcsbm model ------
#> ################# Generation  1: best solution with an ICL of -29391 and 17 clusters #################
#> ################# Generation  2: best solution with an ICL of -29308 and 19 clusters #################
#> ################# Generation  3: best solution with an ICL of -29308 and 19 clusters #################
```

The generative model will be chosen automatically to fit with the data
provided. Here Jazz is a square sparse matrix and a dcsbm model will be
used by default. Some plotting function enable the exploraiton of the
clustering results:

``` r
plot(sol)
```

<img src="man/figures/README-plot-1.png" width="60%" />

And the hierarhical structure between clusters:

``` r
plot(sol,type='tree')
```

<img src="man/figures/README-tree-1.png" width="60%" />

Eventually, one may explore some coarser clustering using the cut
function:

``` r
plot(cut(sol,5))
```

<img src="man/figures/README-cut-1.png" width="60%" />

For large datasets, it’s possible to use parallelism to speed-up the
computation thanks to the future package. You only need to specify the
type of parralelism.

``` r
library(future)
plan(multisession)
data("Blogs")
sol=greed(Blogs$X)
#> ------- Fitting a dcsbm model ------
#> ################# Generation  1: best solution with an ICL of -84490 and 16 clusters #################
#> ################# Generation  2: best solution with an ICL of -84209 and 17 clusters #################
#> ################# Generation  3: best solution with an ICL of -84152 and 17 clusters #################
#> ################# Generation  4: best solution with an ICL of -84146 and 17 clusters #################
#> ################# Generation  5: best solution with an ICL of -84146 and 17 clusters #################
plot(sol)
```

<img src="man/figures/README-future-1.png" width="60%" />
