#' Generate a difference that varies linearly over time. Transition times can be
#' equal, in which case there will be no transition (see generate.step.effect).
#' If both are zero, it is effectively a flat difference applied at every time
#' point (see generate.flat.effect).
#'
#' @param time.vector A numeric vector of time values
#' @param diff.amplitude The numeric amplitude of the difference (Default: 1)
#' @param transition.start.time When does the effect start, will not be applied
#'   if NA. (Default: 0)
#' @param transition.end.time When does the effect end, will not be applied
#'   if NA. (Default: 0)
#' @import data.table
#' @return A numeric vector containing the effect at each provided time point.
#' @export
generate.linear.effect = function(time.vector,
                                  diff.amplitude = 1,
                                  transition.start.time = 0,
                                  transition.end.time = 0) {

  if (is.na(transition.end.time)) {
    transition.end.time = 0
    diff.amplitude = 0
  }
  if (is.na(transition.start.time)) {
    transition.start.time = transition.end.time
  }


  if (transition.start.time > transition.end.time) {
    stop(paste0('transition.start.time (', transition.start.time, ') is greater than transition.end.time (', transition.end.time, ')'))
  }
  else if (transition.end.time == transition.start.time){
    transition.start.time = transition.end.time - 1
  }
  effects = (time.vector - transition.start.time + 1) * diff.amplitude/(transition.end.time - transition.start.time)
  effects = unlist(lapply(effects, function(x) min(x,diff.amplitude)))
  effects = unlist(lapply(effects, function(x) max(x,0)))
}

#' A difference that is a step, special case of generate.linear.effect where
#' transition.start and transition.end are equal.
#'
#' @param time.vector A numeric vector of time values
#' @param diff.amplitude The numeric amplitude of the difference (Default: 1)
#' @param step.time When does the step occur (Default: 0)
#' @return A numeric vector containing the effect at each provided time point.
#' @export
generate.step.effect = function(time.vector,
                                diff.amplitude = 1,
                                step.time = 0) {
  generate.linear.effect(time.vector = time.vector,
                             diff.amplitude = diff.amplitude,
                             transition.start.time = step.time,
                             transition.end.time = step.time)
}

#' A difference that is flat, special case of generate.step.effect where
#' step.effect is right at zero.
#'
#' @param time.vector A numeric vector of time values
#' @param diff.amplitude The numeric amplitude of the difference (Default: 1)
#' @return A numeric vector containing the effect at each provided time point.
#' @export
generate.flat.effect = function(time.vector,
                                diff.amplitude = 1) {
  generate.step.effect(time.vector = time.vector,
                       diff.amplitude = diff.amplitude,
                       step.time = 0)
}

#' Generates a cyclical effect based on a cosine function.
#'
#' @param time.vector A numeric vector of time values
#' @param diff.amplitude The numeric amplitude of the difference (Default: 1)
#' @param max.diff.time The time at which the amplitude of the function will be
#'   maximal (Default: 0)
#' @param max.diff.time The period of the cycle (Default: 50)
#' @return A numeric vector containing the effect at each provided time point.
#' @export
generate.cyclical.effect = function(time.vector,
                                    diff.amplitude = 1,
                                    max.diff.time = 0,
                                    period = 50) {
  diff.amplitude * cos(((2*pi)/period) * time.vector - max.diff.time)
}

#' Generates random normal noise. Doesn't really need a time vector, other than
#' to know how many data points there are.
#'
#' @param time.vector A numeric vector of time values
#' @param noise.sd Standard deviation of the noise (Default: 0.5)
#' @return A numeric vector containing the effect at each provided time point.
#' @export
generate.individual.noise = function(time.vector, noise.sd = 0.5) {
  rnorm(length(time.vector), mean = 0, sd = abs(noise.sd))
}



#' Generates parameters for simulating data for the study e.g. what is the
#' recruitment rate, mean, intervention effect etc in each cluster.
#'
#' @param cluster.dt A study info data.table, containing data about clusters.
#' @param sequence.dt A study info data.table, containing data about sequences.
#' @param ppt.per.unit.time.mean What is the mean recruitment rate across
#'   clusters per one time.
#' @param target.n.participants Integer how many participants, will adjust
#'   ppt.per.unit.time.mean in sim.parameters.dt by weighting them to get the
#'   desired number (Default: NA)
#' @param ppt.per.unit.time.sd What is the variance in the recruitment rate
#'   between clusters.
#' @param ppt.per.unit.time.sd.proportion.of.mean What is the variance in the
#'   recruitment rate between clusters as a multiple of mean. Will override
#'   ppt.per.unit.time.sd, disabled by default.
#' @param base.outcome The outcome at every cluster before adding any effects
#' @param intervention.effect.mean The mean effect of the intervention across
#'   cluster.
#' @param intervention.effect.sd The variance in the  effect of the intervention
#'   between clusters.
#' @param intervention.effect.sd.proportion.of.mean The variance in the  effect
#'   of the intervention between clusters as a multiple of mean. Will override
#'   intervention.effect.sd, disabled by default.
#' @param cluster.effect.mean The mean effect associated with each cluster.
#' @param cluster.effect.sd The variance in the  effect associated with cluster.
#' @param cluster.effect.sd.proportion.of.mean The variance in the  effect
#'   associated with cluster as a multiple of mean. Will override
#'   cluster.effect.sd, disabled by default.
#' @param time.effect.per.unit.mean The mean effect of time per one unit across
#'   clusters.
#' @param time.effect.per.unit.sd The variance in the effect of time per one
#'   unit between clusters.
#' @param time.effect.per.unit.sd.proportion.of.mean The variance in the effect
#'   of time per one unit between clusters as a multiple of mean. Will override
#'   time.effect.per.unit.sd, disabled by default.
#' @param individual.noise.mean The mean of the standard deviation used to
#'   generate individual noise across clusters (i.e. generated by rnorm).
#' @param individual.noise.sd The variance in the standard deviation used to
#'   generate individual noise across clusters.
#' @param individual.noise.sd.proportion.of.mean The variance in the standard
#'   deviation used to generate individual noise across clusters as a multiple
#'   of mean. Will override individual.noise.sd, disabled by default.
#' @param cycle.amplitude.mean The mean amplitude of the cyclical effect across
#'   clusters.
#' @param cycle.amplitude.sd The variance in the amplitude of the cyclical
#'   effect between clusters.
#' @param cycle.amplitude.sd.proportion.of.mean The variance in the amplitude of
#'   the cyclical effect between clusters as a multiple of mean. Will override
#'   cycle.amplitude.sd, disabled by default.
#' @param cycle.start.mean The mean cycle start time across clusters.
#' @param cycle.start.sd The variance in the cycle start time between clusters.
#' @param cycle.amplitude.sd.proportion.of.mean The variance in the cycle start
#'   time between clusters as a multiple of mean. Will override cycle.start.sd,
#'   disabled by default.
#' @param cycle.period.mean The mean cycle period across clusters.
#' @param cycle.period.sd The variance in the cycle period between clusters.
#' @param cycle.period.sd.proportion.of.mean The variance in the cycle period
#'   between clusters as a multiple of mean. Will override cycle.period.sd,
#'   disabled by default.
#' @param n.to.generate Integer how many different effects to generate. If more
#'   than one, an integer column will be added (sim.number) to delineate.
#'   (Default: 1)
#' @return The provided cluster.dt parameter, with simulation parameters added
#'   as columns.
#' @export
generate.sim.parameters.dt = function(cluster.dt,
                                      sequence.dt,
                                      #How many participants?
                                      ppt.per.unit.time.mean = 1,
                                      target.n.participants = NA,
                                      ppt.per.unit.time.sd = 0.1,
                                      ppt.per.unit.time.sd.proportion.of.mean = NA,
                                      #Base outcomes
                                      base.outcome = 10,
                                      intervention.effect.mean = 1,
                                      intervention.effect.sd = 0.1,
                                      intervention.effect.sd.proportion.of.mean = NA,
                                      #Site effect
                                      cluster.effect.mean = 0,
                                      cluster.effect.sd = .1,
                                      cluster.effect.sd.proportion.of.mean = NA,
                                      #Time effect per 1 time
                                      time.effect.per.unit.mean = 0.01,
                                      time.effect.per.unit.sd = 0.001,
                                      time.effect.per.unit.sd.proportion.of.mean = NA,
                                      #Individual noise
                                      individual.noise.mean = 0.2,
                                      individual.noise.sd = 0.04,
                                      individual.noise.sd.proportion.of.mean = NA,
                                      #Cyclical effects
                                      cycle.amplitude.mean = 0.8,
                                      cycle.amplitude.sd = 0.1,
                                      cycle.amplitude.sd.proportion.of.mean = NA,
                                      cycle.start.mean = 10,
                                      cycle.start.sd = 0,
                                      cycle.start.sd.proportion.of.mean = NA,
                                      cycle.period.mean = 20,
                                      cycle.period.sd = 0,
                                      cycle.period.sd.proportion.of.mean = NA,

                                      #How many sim parameters to generate?
                                      n.to.generate = 1) {

  #Function so that the generation of simulation parameters can be lapply.
  generate.single.sim.parameters.dt = function(sim.number = NULL) {

    sim.parameters.dt = copy(cluster.dt)

    #Calculate multiples if present.
    if (!is.na(ppt.per.unit.time.sd.proportion.of.mean)) {
      ppt.per.unit.time.sd = ppt.per.unit.time.mean * ppt.per.unit.time.sd.proportion.of.mean
    }
    if (!is.na(intervention.effect.sd.proportion.of.mean)) {
      intervention.effect.sd = intervention.effect.mean * intervention.effect.sd.proportion.of.mean
    }
    if (!is.na(cluster.effect.sd.proportion.of.mean)) {
      cluster.effect.sd = cluster.effect.mean * cluster.effect.sd.proportion.of.mean
    }
    if (!is.na(time.effect.per.unit.sd.proportion.of.mean)) {
      time.effect.per.unit.sd = time.effect.per.unit.mean * time.effect.per.unit.sd.proportion.of.mean
    }
    if (!is.na(individual.noise.sd.proportion.of.mean)) {
      individual.noise.sd = individual.noise.mean * individual.noise.sd.proportion.of.mean
    }
    if (!is.na(cycle.amplitude.sd.proportion.of.mean)) {
      cycle.amplitude.sd = cycle.amplitude.mean * cycle.amplitude.sd.proportion.of.mean
    }
    if (!is.na(cycle.start.sd.proportion.of.mean)) {
      cycle.start.sd = cycle.start.mean * cycle.start.sd.proportion.of.mean
    }
    if (!is.na(cycle.period.sd.proportion.of.mean)) {
      cycle.period.sd = cycle.period.mean * cycle.period.sd.proportion.of.mean
    }

    #Add simulation number of applicable
    if (!is.null(sim.number)) {
      sim.parameters.dt[,sim.number := sim.number]
    }

    #Recruitment rate
    sim.parameters.dt = merge(sim.parameters.dt, sequence.dt, by = 'sequence')
    sim.parameters.dt[, ppt.per.unit.time := rnorm(.N,
                                                   mean = ppt.per.unit.time.mean,
                                                   sd = abs(ppt.per.unit.time.sd))]
    #Set to target recruitment rate if required
    if (!is.na(target.n.participants)) {
      sim.parameters.dt[,ppt.per.unit.time := ppt.per.unit.time / (sum(ppt.per.unit.time * max.time)/2000)]
    }

    #Base outcome pre-intervention
    # print(base.outcome)
    sim.parameters.dt[, base.outcome := base.outcome]
    #Intervention effect
    sim.parameters.dt[, intervention.effect := rnorm(.N,
                                                     mean = intervention.effect.mean,
                                                     sd = abs(intervention.effect.sd))]
    #Cluster effect
    sim.parameters.dt[, cluster.effect := rnorm(.N,
                                                mean = cluster.effect.mean,
                                                sd = abs(cluster.effect.sd))]
    #Time effect
    sim.parameters.dt[, time.effect.per.unit := rnorm(.N,
                                                      mean = time.effect.per.unit.mean,
                                                      sd = abs(time.effect.per.unit.sd))]
    #Individual noise at each site
    sim.parameters.dt[, individual.noise := rnorm(.N,
                                                  mean = individual.noise.mean,
                                                  sd = abs(individual.noise.sd))]
    #Cycle amplitude
    sim.parameters.dt[, cycle.amplitude := rnorm(.N,
                                                 mean = cycle.amplitude.mean,
                                                 sd = abs(cycle.amplitude.sd))]
    #Cycle start
    sim.parameters.dt[, cycle.start := rnorm(.N,
                                             mean = cycle.start.mean,
                                             sd = abs(cycle.start.sd))]
    #Cycle period
    sim.parameters.dt[, cycle.period := rnorm(.N,
                                              mean = cycle.period.mean,
                                              sd = abs(cycle.period.sd))]

    return(sim.parameters.dt)
  }
  if (n.to.generate > 1) {
    sim.parameters.dt = rbindlist(l = mclapply(X = 1:n.to.generate,
                                               FUN = generate.single.sim.parameters.dt))
  } else {
    sim.parameters.dt = generate.single.sim.parameters.dt()
  }

  return(sim.parameters.dt)
}

#' Generates a table of confounding interventions (e.g. those from other
#' interventions or such). There can be a potentially random number of
#' confounding interventions.
#'
#' @param cluster.dt Earliest time of a random effect (Default: 0)
#' @param min.time Earliest time of a random effect (Default: 0)
#' @param max.time Latest time of a random effect (Default: 0)
#' @param n.confound.intervention.mean How many random effects total, mean parameter for
#'   rnorm (Default: 5)
#' @param n.confound.intervention.sd How many random effects total, sd parameter for
#'   rnorm (Default: 1)
#' @param confound.intervention.transition.duration.mean Mean time taken for transition
#'   of effect across all random effects (Default: 5)
#' @param confound.intervention.transition.duration.mean Variance of time taken for
#'   transition of effect between random effects (Default: 5)
#' @param confound.intervention.mean Mean random effect (Default: 0.7)
#' @param confound.intervention.sd Variance of the random effects (Default: 0.1)
#' @param n.to.generate Integer how many different effects to generate. If more
#'   than one, an integer column will be added (sim.number) to delineate.
#'   (Default: 1)
#' @return A numeric vector containing the effect at each provided time point.
#' @export
generate.confound.interventions.parameters.dt = function(cluster.dt,
                                                         min.time = 0,
                                                         max.time = 0,
                                                         n.confound.intervention.mean = 5,
                                                         #The number of random effects
                                                         n.confound.intervention.sd = 1,
                                                         confound.intervention.transition.duration.mean = 5,
                                                         confound.intervention.transition.duration.sd = 0.6,
                                                         confound.intervention.mean = 0.7,
                                                         confound.intervention.sd = 0.1,

                                                         n.to.generate = 1) {


  #Function so that the generation of confounding effect parameters can be lapply.
  generate.single.confound.interventions.parameters.dt = function(sim.number = NULL) {
    #MAke a data.table with zero or more confounding effects, chuck it in the simulation
    #table as a list.
    n.confound.interventions = round(rnorm(n = 1,
                                   mean = n.confound.intervention.mean,
                                   sd = n.confound.intervention.sd))
    n.confound.interventions = max(n.confound.interventions,0)
    confound.interventions.parameters.dt = data.table(
      confound.intervention.transition.duration = rnorm(
        n = n.confound.interventions,
        mean = confound.intervention.transition.duration.mean,
        sd = abs(confound.intervention.transition.duration.sd)
      ),
      confound.intervention.start = sample(
        x = min.time:max.time,
        size = n.confound.interventions,
        replace = T
      ),
      confound.intervention = rnorm(
        n = n.confound.interventions,
        mean = confound.intervention.mean,
        sd = abs(confound.intervention.sd)
      ),
      cluster = sample(x = cluster.dt[, cluster],
                       size = n.confound.interventions,
                       replace = T)
    )


    #Add simulation number of applicable
    if (!is.null(sim.number)) {
      confound.interventions.parameters.dt[,sim.number := sim.number]
    }

    return(confound.interventions.parameters.dt)
  }
  if (n.to.generate > 1) {
    confound.interventions.parameters.dt = rbindlist(l = mclapply(X = 1:n.to.generate,
                                               FUN =generate.single.confound.interventions.parameters.dt))
  } else {
    confound.interventions.parameters.dt = generate.single.confound.interventions.parameters.dt()
  }

  return(confound.interventions.parameters.dt)
}


#' Creates a multiple regression model predicting ICC so that data can be
#' simulated fulfilling criteria.
#'
#' @param target.icc Numeric target ICC ([0-1], Default: 0.02)
#' @param iterations Number of ICC that will be calculated to fit the regression
#'   model (Default: 500)
#' @param sim.parameters.args.list A list of the arguments for
#'   generate.sim.parameters.dt, with names being the argument names, and values
#'   being the values, with vectors for the values that should be optimised.
#' @param confound.interventions.args.list A list of the arguments for
#'   generate.confound.interventions.parameters.dt, with names being the argument names,
#'   and values being the values, with vectors for the values that should be
#'   optimised.
#' @param interactions.list A list of lists of interactions to assess (Default:
#'   none)
#' @return A regression model with the provided predictors and ICC as the
#'   outcome.
#' @export
fit.icc.model = function(target.icc = 0.02,
                         sim.parameters.args.list = NULL,
                         confound.interventions.args.list = NULL,
                         interactions.list = NULL) {

}


#' Generate a data table populated with participants from simulated recruitment
#' rates,
#'
#' @param sim.parameters.dt The data.table generated by
#'   generate.sim.parameters.dt
#' @return A data.table including participant IDs, cluster, a base outcome and
#'   time.
#' @export
generate.base.data.dt = function(sim.parameters.dt, min.time = 0, max.time = NULL) {
  setorder(sim.parameters.dt[, .(outcome = rep(.SD[, base.outcome], #Create the correct number of participants for each cluster
                                               each = round(.SD[, ppt.per.unit.time] * max.time))),
                             by = cluster,
                             .SDcols = c('cluster', 'ppt.per.unit.time', 'base.outcome')][, time := sample(x = min.time:max.time,
                                                                                                           #Give participants uniformly random time.
                                                                                                           size = .N,
                                                                                                           replace = T)], cluster, time)[, participant := 1:.N
                                                                                                                                         ]
}

#' Generates a data.table that describes the sequence of interventions in
#' clusters.
#'
#' @param n.sequence Integer number of sequences.
#' @param duration.of.period Integer length of period.
#' @param transition.period Integer transition period. This is the number of
#'   days before intervention that are considered transitional (Default: 0, no
#'   transition)
#' @param intervention.sequence Boolean add a sequence that is entirely
#'   intervention (Default: F)
#' @param control.sequence Boolean add a sequence that is entirely
#'   control (Default: F)
#' @return A data.table with sequence number, transition period start, and
#'   intervention start.
#' @export
generate.sequence.dt = function(n.sequence,
                                duration.of.period,
                                transition.period = 0,
                                intervention.sequence = F,
                                control.sequence = F) {


  #Generate table
  sequence.dt = data.table(
    sequence = factor(paste('Sequence', 1:(n.sequence)),
                      levels = paste('Sequence', 1:(n.sequence)))
  )
  # n.sequence = n.sequence - intervention.sequence - control.sequence

  #Create intervention times, with an intervention period if requested.
  sequence.dt[,intervention.time := duration.of.period * (1:(n.sequence) - intervention.sequence)]
  sequence.dt[,transition.time := intervention.time - transition.period]
  if (intervention.sequence) {
    sequence.dt[1,transition.time := 0]
    sequence.dt[1,intervention.time := 0]
  }

  # Create control period if applicable.
  if (control.sequence) {
    sequence.dt[.N,intervention.time := NA]
    sequence.dt[.N,transition.time := NA]
  }
  #
  #Set keys
  setkey(sequence.dt, sequence)

  return(sequence.dt)
}

#' Generates a data.table that describes how clusters are assigned to sequences.
#'
#' @param sequence.dt data.table with information about the sequences, with
#'   columns sequence, transition.time, and intervention.time.
#' @param n.cluster.per.sequence Integer number of clusters assigned to each
#'   sequence.
#' @return A data.table with cluster and sequence ID columns.
#' @export
generate.cluster.dt = function(sequence.dt, n.cluster.per.sequence = 1) {
  #Generate table
  cluster.dt = data.table(cluster = factor(paste(
    'Cluster', 1:(nrow(sequence.dt) * n.cluster.per.sequence)
  ), levels = paste(
    'Cluster', 1:(nrow(sequence.dt) * n.cluster.per.sequence)
  )),
  sequence = factor(paste(
    rep(x = sequence.dt[, sequence], each = n.cluster.per.sequence)
  )))

  #Set keys
  setkey(cluster.dt, sequence)

  return(cluster.dt)
}

#' Creates a data.table with participants and simulates data for them.
#'
#' @param cluster.dt data.table with the correspondence between cluster and
#'   sequence, with columns cluster and sequence.
#' @param sequence.dt data.table with information about the sequences, with
#'   columns sequence, transition.time, and intervention.time.
#' @param sim.parameters.dt The data.table generated by
#'   generate.sim.parameters.dt
#' @param confound.interventions.parameters.dt The data.table generated by
#'   generate.confound.interventions.parameters.dt
#' @param return.intermediaries If true, intermediaries will be returned (return
#'   type will be list), useful for plotting.
#' @return A data.table with simulated data, or a list with simulated data
#'   (data.dt), all separate effects (effects.dt), and summed effects
#'   (summed.effects.dt)
#' @export
generate.data.dt = function(cluster.dt,
                            sequence.dt,
                            sim.parameters.dt,
                            confound.interventions.parameters.dt,
                            return.intermediaries = F,
                            max.time = NULL) {

  #Generate base data
  data.dt = generate.base.data.dt(sim.parameters.dt, max.time = max.time)

  #Stick periods on each observation.
  cut.data.dt.into.periods(data.dt, sequence.dt)

  #Reorder
  setcolorder(data.dt, c("participant", "time", "cluster", "period", "outcome"))

  #Set key
  setkey(data.dt, participant)

  #Calculate all effects incident by participant.
  effects.dt = generate.effects.dt(data.dt, sim.parameters.dt, confound.interventions.parameters.dt)

  #Sum effects for participant.
  sum.effects.dt = calculate.sum.effects.dt(effects.dt)

  #Sum simulated effects for each participant
  data.dt = apply.sum.effects(data.dt, sum.effects.dt)


  if (return.intermediaries) {
    return(list(data.dt = data.dt,
                effects.dt = effects.dt,
                sum.effects.dt = sum.effects.dt))
  } else {
    return(data.dt)
  }
}

#' Creates a data.table of all effects that will be applied.
#'
#' @param data.dt data.table with minimally participant ID and cluster.
#' @param sim.parameters.dt The data.table generated by
#'   generate.sim.parameters.dt
#' @param confound.interventions.parameters.dt The data.table generated by
#'   generate.confound.interventions.parameters.dt
#' @return A data.table including participant IDs, effect names, effect
#'   sizes, and the order in which they will be applied.
#' @export
generate.effects.dt = function(data.dt, sim.parameters.dt, confound.interventions.parameters.dt) {


  #Add intervention, cluster, time, and seasonal effects
  effects.dt = rbindlist(
    list(
      data.table(
        participant = data.dt[,participant],
        name = 'Intervention effect',
        order = 1,
        effect = merge(data.dt,sim.parameters.dt, by = 'cluster')[
          , .(effects = generate.linear.effect(time.vector = time,
                                               diff.amplitude = intervention.effect,
                                               transition.start.time = transition.time,
                                               transition.end.time = intervention.time)),
          by = c('cluster', 'transition.time', 'intervention.time')
          ][,effects]
      ),
      data.table(
        participant = data.dt[,participant],
        name = 'Cluster effect',
        order = 2,
        effect = merge(data.dt,sim.parameters.dt, by = 'cluster')[
          , .(effects = generate.flat.effect(time.vector = time,
                                                 diff.amplitude = cluster.effect)),
          by = c('cluster')
          ][,effects]
      ),
      data.table(
        participant = data.dt[,participant],
        name = 'Time effect',
        order = 3,
        effect = merge(data.dt,sim.parameters.dt, by = 'cluster')[
          , .(effects = generate.linear.effect(time.vector = time,
                                                   transition.start.time = min(time),
                                                   transition.end.time = max(time))),
          by = c('cluster')
          ][,effects]
      ),
      data.table(
        participant = data.dt[,participant],
        name = 'Seasonal effect',
        order = 4,
        effect = merge(data.dt,sim.parameters.dt, by = 'cluster')[
          , .(effects = generate.cyclical.effect(time.vector = time,
                                                     diff.amplitude = cycle.amplitude,
                                                     max.diff.time = cycle.start,
                                                     period = cycle.period)),
          by = c('cluster', 'cycle.amplitude', 'cycle.start', 'cycle.period')
          ][,effects]
      )
    )
  )
  #Add random effects to relevant clusters
  if (nrow(confound.interventions.parameters.dt) > 0) {
    for (confound.intervention.ind in 1:nrow(confound.interventions.parameters.dt)) {
      effects.dt = rbindlist(
        list(
          effects.dt,
          data.table(
            participant = data.dt[cluster==confound.interventions.parameters.dt[confound.intervention.ind,cluster], participant],
            name = paste('Random effect', confound.intervention.ind),
            order = effects.dt[,max(order)+1],
            effect = merge(data.dt, confound.interventions.parameters.dt[confound.intervention.ind], by = 'cluster')[
              , .(effects = generate.linear.effect(time.vector = time,
                                                   diff.amplitude = confound.intervention,
                                                   transition.start.time = confound.intervention.start,
                                                   transition.end.time = confound.intervention.start + confound.intervention.transition.duration
              )
              ),
              by = c('cluster', 'confound.intervention.transition.duration', 'confound.intervention.start', 'confound.intervention')][, effects]
          )
        )
      )
    }
  }

  #Add individual noise
  effects.dt = rbindlist(
    list(
      effects.dt,
      data.table(
        participant = data.dt[, participant],
        name = 'Individual noise',
        order = effects.dt[,max(order)+1],
        effect = merge(data.dt, sim.parameters.dt, by = 'cluster')[
          , .(effects = generate.individual.noise(time.vector = time,
                                                  noise.sd = individual.noise
          )
          ),
          by = c('cluster', 'individual.noise')][, effects]
      )
    )
  )
}


#' Sums the effects for each participant.
#'
#' @param effects.dt The data.table generated by generate.effects.dt
#' @return A data.table including participant IDs, and summed effects.
#' @export
calculate.sum.effects.dt = function(effects.dt) {
  sum.effects.dt = effects.dt[,.(sum.effect = sum(effect)),by = participant]
  setkey(sum.effects.dt, participant)
}


#' Apply the sum of effects to each participant.
#'
#' @param data.dt data.table with minimally participant ID and outcome.
#' @param sum.effects.dt The data.table generated by calculate.sum.effects.dt,
#'   including an effect to be applied to each participant.
#' @return A data.table including participant IDs, and outcomes after summed
#'   effects are applied.
#' @export
apply.sum.effects = function(data.dt, sum.effects.dt) {
  result.dt = merge(data.dt, sum.effects.dt)
  result.dt[, outcome := outcome + sum.effect]
  result.dt[, sum.effect := NULL]
  return(result.dt)
}

#' Converts a continuous variable to a binary one. Zero or negative values of
#' the continuous variable will be considered to have the first level of the
#' variable.
#'
#' @param data.dt data.table with minimally participant ID and outcome.
#' @param max.chance The chance that the highest value of the continuous
#'   variable will be assigned the second level of the binary variable.
#'   (Default: 0.5)
#' @param var.levels Vector containing the values that the binary outcome will
#'   assume (Default: c(0,1))
#' @param delete.outcome.orig Delete original continuous outcome (Default: TRUE)
#' @param delete.prob Delete calculated probability of var.levels[2] (Default:
#'   TRUE)
#' @return A data.table with the continuous outcome column replaced by the
#'   specified binary variable. Not necessary to use return value, as data.dt
#'   will be modified by reference,
#' @export
convert.continuous.to.binary = function(data.dt, max.chance = 0.5,
                                        var.levels = c(0,1),
                                        delete.outcome.orig = T,
                                        delete.prob = T) {
  #No chance less than zero
  data.dt[outcome<0, outcome := 0]
  #Calculate chance of var.levels[2] for each participant
  data.dt[,outcome.chance := outcome * max.chance/max(outcome)]
  #Remove old outcome
  if (!delete.outcome.orig) {
    data.dt[,outcome.orig := outcome]
  }

  data.dt[,outcome := NULL]

  #Calculate new binary outcome with calculated chance
  # print(data.dt[,mean(outcome.orig),by = cluster])
  data.dt[,outcome := sample(var.levels,
                                       size=1,
                                       prob = c(1-outcome.chance, outcome.chance)),
          by = participant]

  #Delete the chance that var.levels[2] would be chosen.
  if (delete.prob) {
    data.dt[,outcome.chance := NULL]
  }
}

#' Generates a number of simulated datasets, and assesses the ICC of each.
#' Parameter lists are as provided to generate.sim.parameters.dt and
#' generate.confound.interventions.parameters.dt, but if there are vectors for
#' any of the values, these will be randomly sampled, generating multiple
#' scenarios.
#'
#' @param sequence.dt data.table with information about the sequences, with
#'   columns sequence, transition.time, and intervention.time.
#' @param cluster.dt data.table with the correspondence between cluster and
#'   sequence, with columns cluster and sequence.
#' @param sim.parameters.args.list List of parameters, as provided to
#'   generate.sim.parameters.dt
#' @param confound.interventions.args.list List of parameters, as provided to
#'   confound.interventions.args.list
#' @param n.iterations Integer number of *total* iterations to assess (i.e. not
#'   the number per scenario)
#' @param binary.data Generate binary data (Default: F)
#' @return A data.table with a column for each simulation parameter, and a
#'   column named ICC added.
#' @export
calculate.icc.dt.of.simulations = function(cluster.dt,
                                           sequence.dt,
                                           sim.parameters.args.list,
                                           confound.interventions.args.list,
                                           n.iterations = 20,
                                           progress.bar = T,
                                           binarise.data = F,
                                           max.time = NULL) {


  #Find what parameters are modifiable.
  sim.modifiable.parameters = names(sim.parameters.args.list[lapply(sim.parameters.args.list,
                                                                    length) > 1])
  confound.modifiable.parameters = names(confound.interventions.args.list[lapply(confound.interventions.args.list,
                                                                                 length) > 1])
  if ((length(sim.modifiable.parameters) + length(confound.modifiable.parameters)) == 0) {
    stop('There are no modifiable parameters in sim.parameters.args.list or confound.interventions.args.list')
  }

  #Generates a data.table with as many rows as iterations of data that will be
  #used for the regression.
  generate.parameters.args.dt = function(parameters.args.list,
                                         n.iterations = n.iterations,
                                         progress.bar = progress.bar) {

    generate.list.of.regression.parameters = function(iteration,
                                                      parameters.args.list,
                                                      progress.bar = NULL) {
      #Function to sample from parameters list, using sample if parameter is a vector.
      sampling.func = function(x) ifelse(length(x) > 1,
                                         yes = sample(x, size = 1, replace = T),
                                         no = x)
      #Get parameters for this iteration.
      parameters.args.list = lapply(parameters.args.list, sampling.func)

      #Update progress bar.
      if (!is.null(progress.bar)) {
        update.progress.bar(progress.bar = progress.bar,
                            index = iteration,
                            max.r = n.iterations-1)
      }

      return(parameters.args.list.dt = as.data.table(parameters.args.list)[,iteration := iteration])
    }

    parameters.args.dt = rbindlist(
      mclapply(
        X = 1:n.iterations,
        FUN = generate.list.of.regression.parameters,
        parameters.args.list = parameters.args.list,
        progress.bar = progress.bar
      )
    )

  }


  message(paste0("Generating simulation parameters for ", n.iterations, " different scenarios..."))
  if (progress.bar == T) {
    pb = txtProgressBar(style = 3, char = '|')
  } else {
    pb = NULL
  }
  t = Sys.time()
  sim.parameters.args.dt = generate.parameters.args.dt(sim.parameters.args.list,
                                                       n.iterations,
                                                       progress.bar = pb)
  #Close progress bar, get time elapsed
  elapsed.time.message(start.time = t,
                       n.iterations = n.iterations,
                       progress.bar = pb)

  message(paste0("Calculating confounding intervention effects for ", n.iterations, " different scenarios..."))
  if (progress.bar == T) {
    pb = txtProgressBar(style = 3, char = '|')
  }
  t = Sys.time()
  confound.interventions.args.dt = generate.parameters.args.dt(confound.interventions.args.list,
                                                               n.iterations,
                                                               progress.bar = pb)
    #Close progress bar, get time elapsed
  elapsed.time.message(start.time = t,
                       n.iterations = n.iterations,
                       progress.bar = pb)

  #Generates a list of data.tables for effects and confounding interventions to be
  #fed to generate.effects.dt
  calculate.icc.from.iteration = function(iteration.num = 1,
                                          sequence.dt = sequence.dt,
                                          cluster.dt = cluster.dt,
                                          progress.bar = NULL) {


    sim.parameters.dt = do.call(generate.sim.parameters.dt,
                                args = c(as.list(sim.parameters.args.dt[iteration.num, !"iteration"]),
                                         cluster.dt = list(cluster.dt),
                                         sequence.dt = list(sequence.dt)))

    confound.interventions.parameters.dt = do.call(generate.confound.interventions.parameters.dt,
                                                   args = c(
                                                     as.list(confound.interventions.args.dt[iteration.num,!"iteration"]),
                                                     cluster.dt = list(cluster.dt)
                                                   ))
                                       # as.list(confound.interventions.args.dt[iteration.num,!"iteration"]))

    data.dt = generate.data.dt(cluster.dt = cluster.dt,
                               sequence.dt = sequence.dt,
                               sim.parameters.dt = sim.parameters.dt,
                               confound.interventions.parameters.dt = confound.interventions.parameters.dt,
                               max.time = max.time)


    #Binarise data if requested and not already done.
    if (binarise.data) {
      bin.data.dt = convert.continuous.to.binary(
        copy(data.dt),
        delete.outcome.orig = T,
        delete.prob = T
      )
      icc = calculate.icc(bin.data.dt)
    } else {
      icc = calculate.icc(data.dt)
    }

    #Update progress bar.
    if (!is.null(progress.bar)) {
      update.progress.bar(progress.bar = progress.bar,
                          index = iteration.num,
                          max.r = n.iterations-1)
    }

    return(icc)
    # return(list(sim.parameters.dt = sim.parameters.dt,
    #             confound.parameters.dt = confound.parameters.dt))
  }

  #Make a wee progress bar.
  message(paste0("Generating and calculating ICC for ", nrow(sim.parameters.args.dt), " different scenarios..."))
  if (progress.bar == T) {
    pb = txtProgressBar(style = 3, char = '|')
  }
  t = Sys.time()

  #Create a data.table with columns being the simulation parameters and calculated ICC.
  icc.dt = cbind(
    merge(sim.parameters.args.dt,
          confound.interventions.args.dt,
          by = 'iteration'),
    icc = unlist(
      mclapply(
        X = 1:n.iterations,
        FUN = calculate.icc.from.iteration,
        sequence.dt = sequence.dt,
        cluster.dt = cluster.dt,
        progress.bar = pb
      )
    )
  )

  #Close progress bar, get time elapsed
  elapsed.time.message(start.time = t,
                       n.iterations = n.iterations,
                       progress.bar = pb)

  return(icc.dt)

}

#' Detects what columns in a data.table have more than a certain number of
#' columns.
#'
#' @param data.dt data.table with minimally participant ID and outcome.
#' @param threshold The number of levels over which the column name will be
#'   returned (Default: 1)
#' @param cols.to.ignore Columns to not count (Default: c('iteration', 'icc'))
#' @return A vector of strings which are the names of columns with numers of
#'   levels exceeding the threshold.
#' @export
cols.exceed.level.threshold = function(data.dt,
                                       threshold = 1,
                                       cols.to.ignore = c('iteration', 'icc')) {
  #Get unique levels of all variables
  n.levels.per.col = apply(data.dt, 2, function(x) length(unique(x)))
  #Get those values that vary.
  vars.to.plot = names(n.levels.per.col)[n.levels.per.col > 1]
  # print(vars.to.plot)
  return(vars.to.plot[!vars.to.plot %in% cols.to.ignore])
}

#' Fits a regression model predicting an outcome variable from provided variables.
#'
#' @param icc.dt A data.table with the parameters of the various simulations and
#'   ICC.
#' @param outcome.var Character column name to predict (Default: ICC)
#' @param predictor.vars String of the variable(s) that will be used as
#'   predictors of ICC, if NULL all variables with more than one level will be
#'   used. (Default: NULL)
#' @param predictor.vars Calculate interation effects, or just simple effects
#'   (Default: F)
#' @param model.type Type of model to fit, currently supports linear ('lm') and
#'   quadratic ('loess') (Default: 'lm')
#' @return A linear model object.
#' @export
fit.regression.model = function(icc.dt,
                                predictor.vars = NULL,
                                outcome.var = 'icc',
                                include.interactions = F,
                                model.type = c('lm', 'loess')[1]) {


  #Auto detect levels to plot against if requested
  if (is.null(predictor.vars)) {
    predictor.vars = cols.exceed.level.threshold(icc.dt,
                                                 threshold = 1,
                                                 cols.to.ignore = c('iteration', 'icc'))
  }

  include.interactions = F
  sep = ifelse(include.interactions, '*', '+')
  model.formula = as.formula(paste(outcome.var,'~', paste(predictor.vars, collapse = sep)))
  if (model.type == 'lm') {
    model = lm(formula = model.formula, data = icc.dt)
  } else if (model.type == 'loess') {
    model = loess(formula = model.formula, data = icc.dt)
  } else {
    stop(paste0('Model type ', model.type, ' not supported.'))
  }

  return(model)
}

#' Fits a regression model predicting ICC from provided variables.
#'
#' @param icc.dt A data.table with the parameters of the various simulations and
#'   ICC.
#' @param predictor.vars String of the variable(s) that will be used as
#'   predictors of ICC, if NULL all variables with more than one level will be
#'   used. (Default: NULL)
#' @param predictor.vars Calculate interation effects, or just simple effects
#'   (Default: F)
#' @param model.type Type of model to fit, currently supports linear ('lm') and
#'   quadratic ('loess') (Default: 'lm')
#' @return A linear model object.
#' @export
fit.icc.regression.model = function(icc.dt,
                                    predictor.vars = NULL,
                                    include.interactions = F,
                                    model.type = c('lm','loess')[1]) {

  fit.regression.model(
    icc.dt = icc.dt,
    predictor.vars = predictor.vars,
    outcome.var = 'icc',
    include.interactions = include.interactions,
    model.type = model.type
  )
}
