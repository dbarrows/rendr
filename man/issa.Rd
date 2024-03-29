% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/issa.R
\name{issa}
\alias{issa}
\title{Inhomogeneous Stochastic Simulation Algorithm (ISSA) solver}
\usage{
issa(
  sys,
  length.out = 100,
  all.out = FALSE,
  trajectories = 1,
  parallel = FALSE,
  cores = detectCores(),
  average = FALSE,
  verbose = TRUE,
  k = NULL,
  force_compile = FALSE
)
}
\arguments{
\item{sys}{an instance of the \code{\link{rdsys}} class}

\item{length.out}{length of solution output (table rows) (default 100)}

\item{all.out}{if \code{TRUE} (default \code{FALSE}), ignore \code{length.out} and return entire solution}

\item{trajectories}{number of trajectories to generate}

\item{parallel}{if \code{TRUE} (default \code{FALSE}) generate trajectories using multiple CPU cores}

\item{cores}{number of cores to use if \code{parallel} is \code{TRUE} (default is all system cores)}

\item{average}{if \code{TRUE} (default \code{FALSE}) and generating multiple trajectories, averages trajectories at sample times; incompatible with \code{all.out = TRUE}}

\item{verbose}{controls if output is generated during during run (default \code{TRUE})}

\item{k}{\code{\link{numeric}} vector of reaction rates corresponding to the reactions in \code{sys}, overrides those contained if \code{sys} if provided}

\item{force_compile}{if \code{TRUE} (default \code{FALSE}), force a recompile of the reaction network}
}
\value{
\code{\link{rdsol}} instance
}
\description{
Inhomogeneous Stochastic Simulation Algorithm (ISSA) solver
}
