% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.helpers.R
\name{permute.clusters.to.sequence}
\alias{permute.clusters.to.sequence}
\title{Permutes clusters to a sequence, using the permutation in column perm.ind
from perm.dt, uses hypothetical.data.dt to assign conditions for speed.}
\usage{
permute.clusters.to.sequence(
  perm.ind,
  perm.dt,
  hypothetical.data.dt,
  exclude.transition = T
)
}
\arguments{
\item{perm.ind}{Number of permutation column in perm.dt that should be used.}

\item{perm.dt}{Precalculated table of how to permute clusters to sequences,
with first column being clusters, and then one column for each permutation
after that.}

\item{hypothetical.data.dt}{data.table with hypothetical data point condition
assignments if each cluster were in each sequence. Three columns, sequence,
cluster (containing all combinations), and data.dt, containing the
hypothetical data.dt from only that cluster.}

\item{exclude.transition}{boolean, should the result exclude data opints from
the transition period? (Default = T)}
}
\value{
A data.table with permuted data and group assignments.
}
\description{
Permutes clusters to a sequence, using the permutation in column perm.ind
from perm.dt, uses hypothetical.data.dt to assign conditions for speed.
}
