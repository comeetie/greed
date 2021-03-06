% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gmm.R
\docType{class}
\name{gmm-class}
\alias{gmm-class}
\title{Gaussian mixture model description class}
\description{
An S4 class to represent a multivariate Gaussian mixture  model, extend \code{\link{icl_model-class}}.
The model corresponds to the following generative model:
\deqn{ \pi \sim Dirichlet(\alpha)}
\deqn{ Z_i  \sim \mathcal{M}(1,\pi)}
\deqn{ V_k \sim \mathcal{W}(\varepsilon^{-1},n_0)}
\deqn{ \mu_k \sim \mathcal{N}(\mu,(\tau V_k)^{-1})}
\deqn{ X_{i}|Z_{ik}=1 \sim \mathcal{N}(\mu_k,V_{k}^{-1})}
with \eqn{\mathcal{W}(\varepsilon^{-1},n_0)} the Whishart distribution.
}
\section{Slots}{

\describe{
\item{\code{name}}{name of the model}

\item{\code{alpha}}{Dirichlet over cluster proportions prior parameter (default to 1)}

\item{\code{tau}}{Prior parameter (inverse variance) default 0.01}

\item{\code{N0}}{Prior parameter (pseudo count) should be > number of features (default to NaN, in this case it will be estimated from data as the number of columns of X)}

\item{\code{epsilon}}{Prior parameter co-variance matrix prior (matrix of size D x D), (default to a matrix of NaN, in this case epsilon will be estimated from data and will corresponds to 0.1 times a diagonal matrix with the variances of the X columns)}

\item{\code{mu}}{Prior parameters for the means (vector of size D), (default to NaN, in this case mu will be estimated from the data and will be equal to the mean of X)}
}}

\examples{
new("gmm")
new("gmm",alpha=1,tau=0.1,N0=15)
}
\references{
Bertoletti, Marco & Friel, Nial & Rastelli, Riccardo. (2014). Choosing the number of clusters in a finite mixture model using an exact Integrated Completed Likelihood criterion. METRON. 73. 10.1007/s40300-015-0064-5. #'
}
