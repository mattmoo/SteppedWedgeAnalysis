% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.helpers.R
\name{generate.sim.parameters.dt}
\alias{generate.sim.parameters.dt}
\title{Generates parameters for simulating data for the study e.g. what is the
recruitment rate, mean, intervention effect etc in each cluster.}
\usage{
generate.sim.parameters.dt(
  cluster.dt,
  sequence.dt,
  ppt.per.unit.time.mean = 1,
  target.n.participants = NA,
  ppt.per.unit.time.sd = 0.1,
  ppt.per.unit.time.sd.proportion.of.mean = NA,
  base.outcome = 10,
  intervention.effect.mean = 1,
  intervention.effect.sd = 0.1,
  intervention.effect.sd.proportion.of.mean = NA,
  cluster.effect.mean = 0,
  cluster.effect.sd = 0.1,
  cluster.effect.sd.proportion.of.mean = NA,
  time.effect.per.unit.mean = 0.01,
  time.effect.per.unit.sd = 0.001,
  time.effect.per.unit.sd.proportion.of.mean = NA,
  individual.noise.mean = 0.2,
  individual.noise.sd = 0.04,
  individual.noise.sd.proportion.of.mean = NA,
  cycle.amplitude.mean = 0.8,
  cycle.amplitude.sd = 0.1,
  cycle.amplitude.sd.proportion.of.mean = NA,
  cycle.start.mean = 10,
  cycle.start.sd = 0,
  cycle.start.sd.proportion.of.mean = NA,
  cycle.period.mean = 20,
  cycle.period.sd = 0,
  cycle.period.sd.proportion.of.mean = NA,
  min.time = 0,
  max.time = NA,
  n.to.generate = 1
)
}
\arguments{
\item{cluster.dt}{A study info data.table, containing data about clusters.}

\item{sequence.dt}{A study info data.table, containing data about sequences.}

\item{ppt.per.unit.time.mean}{What is the mean recruitment rate across
clusters per one time.}

\item{target.n.participants}{Integer how many participants, will adjust
ppt.per.unit.time.mean in sim.parameters.dt by weighting them to get the
desired number (Default: NA)}

\item{ppt.per.unit.time.sd}{What is the variance in the recruitment rate
between clusters.}

\item{ppt.per.unit.time.sd.proportion.of.mean}{What is the variance in the
recruitment rate between clusters as a multiple of mean. Will override
ppt.per.unit.time.sd, disabled by default.}

\item{base.outcome}{The outcome at every cluster before adding any effects}

\item{intervention.effect.mean}{The mean effect of the intervention across
cluster.}

\item{intervention.effect.sd}{The variance in the  effect of the intervention
between clusters.}

\item{intervention.effect.sd.proportion.of.mean}{The variance in the  effect
of the intervention between clusters as a multiple of mean. Will override
intervention.effect.sd, disabled by default.}

\item{cluster.effect.mean}{The mean effect associated with each cluster.}

\item{cluster.effect.sd}{The variance in the  effect associated with cluster.}

\item{cluster.effect.sd.proportion.of.mean}{The variance in the  effect
associated with cluster as a multiple of mean. Will override
cluster.effect.sd, disabled by default.}

\item{time.effect.per.unit.mean}{The mean effect of time per one unit across
clusters.}

\item{time.effect.per.unit.sd}{The variance in the effect of time per one
unit between clusters.}

\item{time.effect.per.unit.sd.proportion.of.mean}{The variance in the effect
of time per one unit between clusters as a multiple of mean. Will override
time.effect.per.unit.sd, disabled by default.}

\item{individual.noise.mean}{The mean of the standard deviation used to
generate individual noise across clusters (i.e. generated by rnorm).}

\item{individual.noise.sd}{The variance in the standard deviation used to
generate individual noise across clusters.}

\item{individual.noise.sd.proportion.of.mean}{The variance in the standard
deviation used to generate individual noise across clusters as a multiple
of mean. Will override individual.noise.sd, disabled by default.}

\item{cycle.amplitude.mean}{The mean amplitude of the cyclical effect across
clusters.}

\item{cycle.amplitude.sd}{The variance in the amplitude of the cyclical
effect between clusters.}

\item{cycle.amplitude.sd.proportion.of.mean}{The variance in the cycle start
time between clusters as a multiple of mean. Will override cycle.start.sd,
disabled by default.}

\item{cycle.start.mean}{The mean cycle start time across clusters.}

\item{cycle.start.sd}{The variance in the cycle start time between clusters.}

\item{cycle.period.mean}{The mean cycle period across clusters.}

\item{cycle.period.sd}{The variance in the cycle period between clusters.}

\item{cycle.period.sd.proportion.of.mean}{The variance in the cycle period
between clusters as a multiple of mean. Will override cycle.period.sd,
disabled by default.}

\item{min.time}{Minimum time value (Default: 0)}

\item{max.time}{Maximum time value, if NA will be calculated from
sequence.dt, by assuming equal period lengths, which may not be correct
(Default: NA)}

\item{n.to.generate}{Integer how many different effects to generate. If more
than one, an integer column will be added (sim.number) to delineate.
(Default: 1)}
}
\value{
The provided cluster.dt parameter, with simulation parameters added
  as columns.
}
\description{
Generates parameters for simulating data for the study e.g. what is the
recruitment rate, mean, intervention effect etc in each cluster.
}
