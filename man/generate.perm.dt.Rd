% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.helpers.R
\name{generate.perm.dt}
\alias{generate.perm.dt}
\title{Generate a table that holds a number of ways in which to permute sites to
different clusters.}
\usage{
generate.perm.dt(cluster.dt, max.r = 1000)
}
\arguments{
\item{cluster.dt}{A study info data.table, containing cluster and sequence
relations.}

\item{max.r}{Number of permutations.}
}
\value{
A data.table containing all the permutations of site to cluster.
}
\description{
Generate a table that holds a number of ways in which to permute sites to
different clusters.
}
