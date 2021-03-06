% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.helpers.R
\name{generate.hypothetical.data.dt}
\alias{generate.hypothetical.data.dt}
\title{Precalculate the intervention and control groups for each site for each
cluster (i.e. sequence). This could fail for very large data sets.}
\usage{
generate.hypothetical.data.dt(
  data.dt,
  cluster.dt,
  sequence.dt,
  progress.bar = T
)
}
\arguments{
\item{data.dt}{data.table with columns participant, cluster, time, and
outcome. Outcome should be continuous.}

\item{cluster.dt}{data.table with the correspondence between cluster and
sequence, with columns cluster and sequence.}

\item{sequence.dt}{data.table with information about the sequences, with
columns sequence, transition.time, and intervention.time.}

\item{progress.bar}{Display a progress bar. A little bit of overhead.
Completion times will be echoed regardless.}

\item{...}{Passed on to the test}
}
\value{
A data.table with a row for each combination of cluster and sequence,
  and a data.table of the participants for that cluster assigned on that
  basis.
}
\description{
Precalculate the intervention and control groups for each site for each
cluster (i.e. sequence). This could fail for very large data sets.
}
