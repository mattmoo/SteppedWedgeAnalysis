% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.helpers.R
\name{generate.cyclical.effect}
\alias{generate.cyclical.effect}
\title{Generates a cyclical effect based on a cosine function.}
\usage{
generate.cyclical.effect(
  time.vector,
  diff.amplitude = 1,
  max.diff.time = 0,
  period = 50
)
}
\arguments{
\item{time.vector}{A numeric vector of time values}

\item{diff.amplitude}{The numeric amplitude of the difference (Default: 1)}

\item{max.diff.time}{The period of the cycle (Default: 50)}
}
\value{
A numeric vector containing the effect at each provided time point.
}
\description{
Generates a cyclical effect based on a cosine function.
}
