% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SL_stabs_fitfun.R
\name{SL_stabs_fitfun}
\alias{SL_stabs_fitfun}
\title{Wrapper for using Super Learner-based extrinsic selection within stability selection}
\usage{
SL_stabs_fitfun(x, y, q, ...)
}
\arguments{
\item{x}{the features.}

\item{y}{the outcome of interest.}

\item{q}{the number of features to select on average.}

\item{...}{other arguments to pass to \code{SuperLearner}.}
}
\value{
a named list, with elements: \code{selected} (a logical vector
indicating whether or not each variable was selected); and \code{path} (
a logical matrix indicating which variable was selected at each step).
}
\description{
A wrapper function for Super Learner-based extrinsic variable selection within
stability selection, using the \code{stabs} package.
}
\seealso{
\code{\link[stabs]{stabsel}} for general usage of stability selection.
}
