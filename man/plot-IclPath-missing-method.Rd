% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alg_classes.R
\name{plot,IclPath,missing-method}
\alias{plot,IclPath,missing-method}
\title{Plot an \code{\link{IclPath-class}} object}
\usage{
\S4method{plot}{IclPath,missing}(x, type = "tree")
}
\arguments{
\item{x}{a \code{\link{IclPath-class}}}

\item{type}{a string which specify plot type:
\itemize{
\item \code{'front'}: plot the extracted front ICL, log(alpha)
\item \code{'path'}: plot the evolution of ICL with respect to K
\item \code{'tree'}: plot the associated dendrogram
}}
}
\value{
a \code{\link[ggplot2]{ggplot}} graphic
}
\description{
Plot an \code{\link{IclPath-class}} object
}
