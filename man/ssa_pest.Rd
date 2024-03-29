% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ssa.R
\name{ssa_pest}
\alias{ssa_pest}
\title{Optimised point estimate (average final trajectory point) using SSA}
\usage{
ssa_pest(
  sys,
  trajectories = 1,
  k = NULL,
  force_compile = FALSE,
  rng_seed = NULL
)
}
\arguments{
\item{sys}{\code{\link{rsys}} instance}

\item{trajectories}{number of trajectories to generate}

\item{k}{\code{\link{vector}} of reaction rates corresponding to the reactions in \code{sys}, overrides those contained if \code{sys} if provided}

\item{force_compile}{if set to \code{TRUE}, forces the overwriting and recompilation of the network source file}

\item{rng_seed}{random number generator seed to use, if provided (default \code{NULL})}
}
\value{
\code{\link{numeric}} vector
}
\description{
Optimised point estimate (average final trajectory point) using SSA
}
