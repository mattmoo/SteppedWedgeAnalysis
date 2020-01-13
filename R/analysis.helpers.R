#' Generate a table that holds a number of ways in which to permute sites to
#' different clusters.
#'
#' @param cluster.dt A study info data.table, containing cluster and sequence
#'   relations.
#' @param max.r Number of permutations.
#' @return A data.table containing all the permutations of site to cluster.
#'
#' @export
generate.perm.dt = function(cluster.dt, max.r = 1000) {
  #Find out roughly how many permutations are possible.
  max.poss.r = factorial(nrow(cluster.dt))

  #Adjust max iterations if necessary.
  if (max.poss.r < max.r) {
    max.r = max.poss.r
  }

  #If there aren't heaps more possible iterations than the number requested,
  if (max.r*10 > max.poss.r) {
    #Get all permutations
    perm.dt = transpose(do.call(data.table,combinat::permn(cluster.dt$sequence)))
    #Get a random sample of the number that you really want.
    perm.dt = perm.dt[sample(x = 1:nrow(perm.dt), size = max.r)]
  } else {
    #If there are heaps of possible permutations, get a random sample without ensuring that they're unique, it's probably good enough.
    perm.dt = transpose(do.call(data.table,lapply(rep(1,max.r), function(n) sample(cluster.dt$sequence))))
  }

  #Now clusters are rows, with each permutation column containing sequences.
  perm.dt = cbind(data.table(cluster = cluster.dt$cluster), transpose(perm.dt))
  setkey(perm.dt, cluster)

  return(perm.dt)
}

#' Do a wilcox test on a data.table. More efficient that coin for large numbers,
#' sort the input by outcome to make it faster.
#'
#' @param data.dt data.table with two-level grouping factor 'condition' and
#'   outcome 'outcome''
#' @param two.sided Do a two sided test?
#' @param tie.correction Apply tie correction (slightly slower, but sometimes
#'   necessary). If NULL, tie correction will be applied if there are.
#' @param sort.input Sorts the input by reference, this enables faster ranking
#'   (i.e. speeds permutation test)
#' @return A data.table containing statistics, including z score and theoretical
#'   p value (e.g. z.dt$z).
#'
#' @export
test.wilcox.dt = function(data.dt,
                          two.sided = T,
                          tie.correction = NULL,
                          sort.input = F) {

  #What proportion of values should be unique to apply tie correction (if not
  #explicitly enabled/disabled)
  tie.correction.threshold = .9

  t = Sys.time()

  if (sort.input) {
    setorder(data.dt, outcome)
  }
  if (is.null(tie.correction)) {
    tie.correction = length(data.dt[,unique(outcome)]) < nrow(data.dt)*tie.correction.threshold
  }

  # message(paste('t0 = ', Sys.time()-t ))

  #Do WMW
  #Assign ranks to outcomes
  data.dt[, rank := frank(outcome, ties.method = 'average')]


  #Make a data.table with rank sums and n per group.
  wmw.dt = data.dt[,.(R = sum(rank), n = .N), by = condition]

  #Calculate U for each group
  data.table::set(x = wmw.dt,
                  j = 'u',
                  value = wmw.dt[,R - (n*(n+1))/2])

  #Return if no valid comparison.
  if (wmw.dt[,.N] == 1) {
    return(NULL)
  }

  #Get minimum U, n for each group, and also theoretical mean.
  z.dt = cbind(wmw.dt[u == min(u), .(n.a=as.numeric(n), u = u)],
               wmw.dt[u == max(u), .(n.b=as.numeric(n))],
               wmw.dt[,.(m.u = prod(n)/2)])
  # print(z.dt)

  #Calculate theoretical SD, correcting for errors if requested.
  if (tie.correction == F) {
    z.dt[,sd.u := wmw.dt[,sqrt(prod(n)*(sum(n)+1)/12)]]
  } else {
    z.dt[,sd.u := sqrt((n.a*n.b/((n.a+n.b)*(n.a+n.b-1))) * ((((n.a+n.b)^3-(n.a+n.b))/12)-sum(data.dt[,(.N^3-.N)/12, by = outcome])))]
  }


  #Calculate z score for U
  data.table::set(x = z.dt,
                  j = 'z',
                  value = z.dt[,(u - m.u)/sd.u])


  data.table::set(x = z.dt,
                  j = 'p',
                  value = z.dt[,pnorm(z) * (1+two.sided)])

  return(z.dt)
}


#' Do a difference in means test on a data.table, as outlined in Thompson,
#' Davey, et al 2018
#'
#' @param data.dt data.table with two-level grouping factor 'condition' and
#'   outcome 'outcome'
#' @return A data.table containing mean difference by condition for each period
#'   and a variance weighting for calculating the mean, along with means,
#'   variances, and numbers of periods for each level of condition.
#'
#' @export
test.f.effect.dt = function(data.dt) {

  #Get summary of variable for each cluster
  cluster.summ.dt = data.dt[, .(cp.mean = mean(outcome, na.rm = T)),
                            by = .(condition, cluster)]

  #Summarise the summary variable, getting mean, number of summaries, and
  #variance.
  slice.long.dt = cluster.summ.dt[, .(
    mean = mean(cp.mean, na.rm = T),
    n = .N,
    var = var(cp.mean, na.rm = T)
  ),
  by = .(condition)]


  if (slice.long.dt[,.N] == 1) {
    return(NULL)
  }

  #Make condition numeric
  slice.long.dt[,condition.num := factor(as.numeric(condition))]

  #Cast to wide to compare two conditions
  summary.dt = dcast(data = slice.long.dt,
                     formula = . ~ condition.num,
                     value.var = c('mean','n','var'))
  summary.dt[,`.` := NULL]

  #Remove periods with no comparison
  summary.dt = summary.dt[!is.na(mean_1) & !is.na(mean_2)]

  #Get difference in means and variance weight
  summary.dt[, diff.mean := mean_2 - mean_1]
  summary.dt[, wgt :=
               ((((n_1 - 1) * var_1 + (n_2 - 1) * var_2) /
                   (n_1 + n_2 - 2)) * (1 / n_1 + 1 / n_2)) ^ -1]

}


#' Make a period column in a data.table, with each period starting at an
#' intervention time.
#'
#' @param data.dt data.table with columns participant, cluster, time, and
#'   outcome. Outcome should be continuous.
#' @param sequence.dt data.table with information about the sequences, with
#'   columns sequence, transition.time, and intervention.time.
#' @return A data.table with period column assigned. This will be done by
#'   reference for speed, so you don't need to use the return value.
#'
#' @export
cut.data.dt.into.periods = function(data.dt, sequence.dt) {
  #Stick periods on each observation.
  sequence.dt = sequence.dt[!is.na(intervention.time)]
  data.dt[, period := cut(
    x = data.dt[, time],
    breaks = c(-Inf, sequence.dt[,intervention.time], Inf),
    right = F,
    labels = c(1:(nrow(sequence.dt) + 1))
  )]
}

#' Assigns control, transition, and intervention condition column on the basis
#' of experiment info. Can optionally include transition period data points.
#'
#' @param data.dt data.table with columns participant, cluster, time, and
#'   outcome. Outcome should be continuous.
#' @param sequence.dt data.table with information about the sequences, with
#'   columns sequence, transition.time, and intervention.time.
#' @param cluster.dt data.table with the correspondence between cluster and
#'   sequence, with columns cluster and sequence.
#' @param exclude.transition boolean, should the result exclude data points from
#'   the transition period? (Default = T)
#' @param cut.labels Character vector of labels for control, transition, and
#'   intervention periods. Probably don't change this (Default:
#'   c('control','transition','intervention'))
#' @return The input data.dt with conditions attached.
#' @export
cut.data.dt.into.conditions = function(data.dt,
                                       cluster.dt,
                                       sequence.dt,
                                       exclude.transition = T,
                                       cut.labels = c('control',
                                                      'transition',
                                                      'intervention')) {

  #Remove sequence and cluster columns to be replaced if present.
  data.dt = data.dt[,c(names(data.dt)[!names(data.dt) %in% c(names(sequence.dt),names(cluster.dt))],
                       'cluster'), with = F]

  #Hack so that there is no malfunctioning when a cluster has no intervention.
  sequence.dt = copy(sequence.dt)[is.na(transition.time), transition.time := max(data.dt[,time]*1000) ]
  sequence.dt = copy(sequence.dt)[is.na(intervention.time), intervention.time := max(data.dt[,time]*1000) ]

  #Do a bit of merging with the cluster and sequence info.
  cut.data.dt = merge(merge(data.dt,
                            cluster.dt,
                            by = 'cluster',
                            all.x = T),
                      sequence.dt,
                      by = 'sequence',
                      all.x = T)
  #Attach condition labels
  cut.data.dt[, condition := cut.labels[1]]
  cut.data.dt[time >= transition.time, condition := cut.labels[2]]
  cut.data.dt[time >= intervention.time, condition := cut.labels[3]]
  cut.data.dt[, condition := as.factor(condition)]
  # [, condition := cut(
  #                       time,
  #                       breaks = c(-Inf,
  #                                  unique(transition.time),
  #                                  unique(intervention.time),
  #                                  Inf),
  #                       right = F,
  #                       labels = cut.labels
  #                     ),
  #                     by = c('cluster')]

  if (exclude.transition) {
    cut.data.dt = cut.data.dt[condition != "transition"]
    cut.data.dt[,condition := droplevels(condition)]
  }

  return(cut.data.dt)
}


#' Cut timepoints from a single cluster into separate conditions, return list of
#' conditions.
#'
#' @param time numeric vector of times.
#' @param transition.time numeric indicating when the transition period should
#'   start.
#' @param intervention.time numeric indicating when the intervention period
#'   should start.
#' @param time.start numeric indicating when labelling the transition time
#'   should start. (Default = -Inf)
#' @param time.end numeric indicating when labelling the intervention time
#'   should end. (Default = Inf)
#' @param cut.labels What should the labels be? (Default:
#'   c('control','transition','intervention'))
#' @return List of conditions in the same order as presented.
#' @export
cut.time.into.conditions = function(time,
                                    transition.time,
                                    intervention.time,
                                    time.start = -Inf,
                                    time.end = Inf,
                                    cut.labels = c('control',
                                                   'transition',
                                                   'intervention')) {

  #Cut on the time column to assign groups.
  cut.breaks = c(time.start,
                 transition.time,
                 intervention.time,
                 time.end)

  result.dt = data.table(time = time,
                         condition = cut.labels[1])
  result.dt[time >= transition.time, condition := cut.labels[2]]
  result.dt[time >= intervention.time, condition := cut.labels[3]]
  result.dt[, condition := as.factor(condition)]


  return(result.dt$condition)
}


#' A function to perform a permutation step, so that the process can be deployed
#' in parallel.
#' @param perm.ind Index of permutation column in perm.dt that should be used.
#' @param perm.dt Precalculated table of how to permute clusters to sequences,
#'   with first column being clusters, and then one column for each permutation
#'   after that.
#' @param data.dt data.table with columns participant, cluster, time, and
#'   outcome. Outcome should be continuous.
#' @param cluster.dt data.table with the correspondence between cluster and
#'   sequence, with columns cluster and sequence.
#' @param sequence.dt data.table with information about the sequences, with
#'   columns sequence, transition.time, and intervention.time.
#' @param hypothetical.data.dt data.table with hypothetical data point condition
#'   assignments if each cluster were in each sequence. Three columns, sequence,
#'   cluster (containing all combinations), and data.dt, containing the
#'   hypothetical data.dt from only that cluster.
#' @param comparison.within Will comparisons be within cluster (i.e. between
#'   period), or within period (i.e. within cluster).
#' @param stat.func What statistic will be used, there is a wilcox and f.effect
#'   as I'm writing this. (Default: Wilcox)
#' @param progress.bar Should be a text progress bar if you want one.
#' @return A data.table with the stat separated by period and permutation
#'   number.
#'
#' @export
perform.permutation.step = function(perm.ind,
                                    perm.dt,
                                    data.dt,
                                    cluster.dt,
                                    sequence.dt,
                                    hypothetical.data.dt,
                                    comparison.within = c('cluster','period')[1],
                                    stat.func = c(test.wilcox.dt,test.f.effect.dt)[[1]],
                                    progress.bar = NULL) {

  #Permute the data, permutation 0 is the actual data.
  if (perm.ind==0) {
    perm.data.dt = copy(data.dt)
  } else {
    perm.data.dt = permute.clusters.to.sequence(perm.ind, perm.dt, hypothetical.data.dt)
  }
  perm.data.dt = cut.data.dt.into.conditions(perm.data.dt, cluster.dt, sequence.dt)

  #Remove rows with no comparison.
  # perm.data.dt = remove.rows.with.no.comparison(perm.data.dt, cluster.dt, sequence.dt, remove.type = comparison.within)

  #Put the statistics in a table for later binding.
  if (comparison.within == 'cluster') {
    #########################################
    #Only data.table optimised wmw.u
    perm.stat.dt = stat.func(perm.data.dt)[,perm.num := perm.ind]
  } else if (comparison.within == 'period') {

    #Get a data.table out with one statistic per period.
    perm.stat.dt = perm.data.dt[,stat.func(copy(.SD)),
                                .SDcols=c('condition','outcome','period','cluster'),
                                by=period][
                                  ,perm.num := perm.ind]

  } else {
    stop(paste0('Do not recognise within ', comparison.within, ' comparison.'))
  }

  if (!is.null(progress.bar)) {
    update.progress.bar(progress.bar = progress.bar,
                        index = perm.ind,
                        max.r = ncol(perm.dt)-1)
  }


  return(perm.stat.dt)
}

#' Precalculate the intervention and control groups for each site for each
#' cluster (i.e. sequence). This could fail for very large data sets.
#'
#' @param data.dt data.table with columns participant, cluster, time, and
#'   outcome. Outcome should be continuous.
#' @param cluster.dt data.table with the correspondence between cluster and
#'   sequence, with columns cluster and sequence.
#' @param sequence.dt data.table with information about the sequences, with
#'   columns sequence, transition.time, and intervention.time.
#' @param progress.bar Display a progress bar. A little bit of overhead.
#'   Completion times will be echoed regardless.
#' @param ... Passed on to the test
#' @return A data.table with a row for each combination of cluster and sequence,
#'   and a data.table of the participants for that cluster assigned on that
#'   basis.
#' @export

generate.hypothetical.data.dt = function(data.dt,
                                         cluster.dt,
                                         sequence.dt,
                                         progress.bar = T) {


  hypothetical.data.dt = data.table::data.table(expand.grid(sequence = cluster.dt[,levels(sequence)],
                                                            cluster = cluster.dt[,levels(cluster)]))
  #This could be vectorised, but won't be the most compute-heavy part anyway.

  #Iterate through the different assignments of site to clusters, making tables of
  #what the groups would look like.

  #Make a wee progress bar.
  message(paste0("Calculating ", nrow(hypothetical.data.dt), " hypothetical site data tables..."))
  if (progress.bar == T) {
    pb = txtProgressBar(style = 3, char = '|')
  }
  t = Sys.time()

  for (row.ind in 1L:nrow(hypothetical.data.dt)) {
    #Get data from a single site, for which to calculate who would be in what
    #condition depending on what sequence it were assigned to.
    sequence.cluster.dt = data.dt[cluster == hypothetical.data.dt[row.ind, cluster],]

    condition.vector = cut.time.into.conditions(time = sequence.cluster.dt[,time],
                                                transition.time = sequence.dt[sequence == hypothetical.data.dt[row.ind, sequence]][,c(transition.time)],
                                                intervention.time = sequence.dt[sequence == hypothetical.data.dt[row.ind, sequence]][,c(intervention.time)])

    sequence.cluster.dt[,condition := condition.vector]

    #Add that data.table to the hypothetical data.
    data.table::set(x = hypothetical.data.dt,
                    i = row.ind,
                    j = "data.dt",
                    value = list(list(sequence.cluster.dt)))

    #Update progress
    if (progress.bar == T) {
      update.progress.bar(progress.bar = pb,
                          index = row.ind,
                          max.r = nrow(hypothetical.data.dt))
      # setTxtProgressBar(pb, row.ind/nrow(hypothetical.data.dt))
    }
  }
  #Close progress bar, get time elapsed
  elapsed.time.message(start.time = t,
                       n.iterations = nrow(hypothetical.data.dt),
                       progress.bar = pb)

  return(hypothetical.data.dt)
}


#' Construct a statistic distribution generated by permuting sites to different
#' clusters.
#'
#' Statistic is Wilcoxon-Mann-Whitney at the moment.
#'
#' @param data.dt data.table with columns participant, cluster, time, and
#'   outcome. Outcome should be continuous.
#' @param cluster.dt data.table with the correspondence between cluster and
#'   sequence, with columns cluster and sequence.
#' @param sequence.dt data.table with information about the sequences, with
#'   columns sequence, transition.time, and intervention.time.
#' @param hypothetical.data.dt Precalculated table of how outcomes would be
#'   assigned to groups if clusters were in each different sequence, will be
#'   generated if NULL (Default: NULL)
#' @param perm.dt Precalculated table of how to permute clusters to sequences,
#'   with first column being clusters, and then one column for each permutation
#'   after that.
#' @param comparison.within Will comparisons be within cluster (i.e. between
#'   period), or within period (i.e. within cluster).
#' @param max.r How many permutations?
#' @param sort.input Will sort the input by outcome (by reference), which
#'   slightly speeds ranking
#' @param exclude.transition boolean, should the result exclude data points from
#'   the transition period? (Default = T)
#' @param stat.func What statistic will be used, there is a wilcox and f.effect
#'   as I'm writing this. (Default: Wilcox)
#' @param progress.bar Display a progress bar. A little bit of overhead.
#'   Completion times will be echoed regardless.
#' @param ... Passed on to the test
#' @return A data.table with the statistic value at each permutation (with zero
#'   as the unpermuted comparison).
#' @export
generate.stat.dt = function(data.dt,
                            cluster.dt,
                            sequence.dt,
                            comparison.within = c('cluster','period')[1],
                            hypothetical.data.dt = NULL,
                            perm.dt = NULL,
                            max.r = 1000,
                            sort.input = T,
                            exclude.transition = T,
                            stat.func = c(test.wilcox.dt,test.f.effect.dt)[[1]],
                            progress.bar = T) {

  if (sort.input == T) {
    setorder(data.dt, outcome)
  }

  #It will probably streamline the calculations if I precalculate the intervention
  #and control groups for each site for each cluster (i.e. sequence).
  #This could fail for very large data sets.
  #This could be vectorised, but won't be the most compute-heavy part anyway.
  if (is.null(hypothetical.data.dt)) {
    hypothetical.data.dt = generate.hypothetical.data.dt(data.dt,
                                                         cluster.dt,
                                                         sequence.dt,
                                                         progress.bar = T)
  }

  #Now let's figure out how to permute clusters to sequences. You might provide
  #one if there is already one generated.
  perm.dt = NULL
  if (is.null(perm.dt)) {
    perm.dt = generate.perm.dt(cluster.dt, max.r = max.r)
  }
  if ((ncol(perm.dt)-1) != max.r) {
    message(paste0("Number of permutations in perm.dt: ", (ncol(perm.dt)-1)))
    message(paste0("Number of requested permutations: ", max.r))
    if ((ncol(perm.dt)-1) > max.r) {
      message(paste0("First ", max.r, " permutations will be used."))
    } else {
      message(paste0("max.r changed to "), (ncol(perm.dt)-1))
      max.r = (ncol(perm.dt)-1)
    }
  }

  #Make a wee progress bar.
  message(paste0("Calculating ", max.r, " permutations..."))

  if (progress.bar == T) {
    pb = txtProgressBar(style = 3, char = '|')
  } else {
    pb = NULL
  }
  t = Sys.time()

  #Apply the permutation step, parallel if available.
  stat.dt = rbindlist(mclapply(0L:max.r,
                               function(perm.ind)
                                 perform.permutation.step(
                                   perm.ind,
                                   perm.dt = perm.dt,
                                   data.dt = data.dt,
                                   cluster.dt = cluster.dt,
                                   sequence.dt = sequence.dt,
                                   hypothetical.data.dt = hypothetical.data.dt,
                                   comparison.within = comparison.within,
                                   stat.func = stat.func,
                                   progress.bar = pb
    )))

  #Close progress bar, get time elapsed
  elapsed.time.message(start.time = t,
                       n.iterations = ncol(perm.dt),
                       progress.bar = pb)

  return(stat.dt)
}

#' Calculate a p-value from a list of statistics, try and automatically detect
#' the appropriate one. If there is a column named 'wgt', use that to weight the
#' average per permutation number.
#'
#' @param stat.dt data.table with a column called perm.num which has the
#'   comparison from the actual data with 0
#' @param stat.col character name of the column which is the stat
#' @return Two-sided p-value
#' @export
calculate.p.from.stat.dt = function(stat.dt) {
  #Find the statistic.
  if ('z' %in% colnames(stat.dt)) {
    stat.col = 'z'
  } else if ('diff.mean' %in% colnames(stat.dt)) {
    stat.col = 'diff.mean'
  }
  #Maybe there's a column for weighting summaries.
  if ('wgt' %in% colnames(stat.dt)) {
    wgt.col = 'wgt'
  } else {
    wgt.col = NULL
  }

  if (!is.null(wgt.col)) {
    stat.vector = stat.dt[,weighted.mean(x = get(stat.col),
                                         w = get(wgt.col)),
                          by = perm.num][,V1]
  } else {
    stat.vector = stat.dt[,mean(get(stat.col)),
                          by = perm.num][,V1]
  }
  sum(abs(stat.vector)>=abs(stat.vector[1]))/length(stat.vector)
  # stat.dt[abs(get(stat.col))>=stat.dt[perm.num==0,abs(get(stat.col))],.N/nrow(stat.dt)]
}


#' Permutes clusters to a sequence, using the permutation in column perm.ind
#' from perm.dt, uses hypothetical.data.dt to assign conditions for speed.
#'
#' @param perm.ind Number of permutation column in perm.dt that should be used.
#' @param perm.dt Precalculated table of how to permute clusters to sequences,
#'   with first column being clusters, and then one column for each permutation
#'   after that.
#' @param hypothetical.data.dt data.table with hypothetical data point condition
#'   assignments if each cluster were in each sequence. Three columns, sequence,
#'   cluster (containing all combinations), and data.dt, containing the
#'   hypothetical data.dt from only that cluster.
#' @param exclude.transition boolean, should the result exclude data opints from
#'   the transition period? (Default = T)
#' @return A data.table with permuted data and group assignments.
#' @export
permute.clusters.to.sequence = function(perm.ind,
                                        perm.dt,
                                        hypothetical.data.dt,
                                        exclude.transition = T) {

  #Rename the relevant permutation's column to pick it out.
  old.col.name = names(perm.dt)[perm.ind+1]
  data.table::setnames(x = perm.dt,
                       old = perm.ind+1,
                       new = "sequence.perm")
  #Shuffle the clusters between sites as dictated in that permutation, and
  #pick out the hypothetical condition groups calculated earlier.
  perm.data.dt = rbindlist(
    merge(x = perm.dt[, .(cluster, sequence.perm)],
          y = hypothetical.data.dt,
          by.x = c("cluster", "sequence.perm"),
          by.y = c("cluster", "sequence")
    )$data.dt)

  if (exclude.transition) {
    perm.data.dt = perm.data.dt[condition != "transition"]
  }

  #Reset name to what it was before
  data.table::setnames(x = perm.dt,
                       old = perm.ind+1,
                       new = old.col.name)


  return(perm.data.dt)
}


#' Removes rows from periods/clusters that do not have participants in both the
#' control and intervention conditions. Needs the conditions attached, so will
#' attempt to do so if not already done.
#'
#' @param data.dt data.table with columns participant, cluster, time, and
#'   outcome. Outcome should be continuous.
#' @param sequence.dt data.table with information about the sequences, with
#'   columns sequence, transition.time, and intervention.time.
#' @param cluster.dt data.table with the correspondence between cluster and
#'   sequence, with columns cluster and sequence.
#' @param remove.type character indicating what to compare on.
#' @return data.dt with rows removed that did not have a corresponding control
#'   or intervention.
#' @export
remove.rows.with.no.comparison = function(data.dt,
                                          cluster.dt = NULL,
                                          sequence.dt = NULL,
                                          remove.type = c('cluster', 'period')[1]) {

  if (!'condition' %in% names(data.dt) & !is.null(cluster.dt) & !is.null(sequence.dt)) {

    data.dt = cut.data.dt.into.conditions(data.dt = data.dt,
                                          cluster.dt = cluster.dt,
                                          sequence.dt = sequence.dt)

  } else {
    if (!'condition' %in% names(data.dt)) {
      stop('Need to have conditions attached to data.dt to see what conditions participants are in.')
    }
  }


  data.dt = data.dt[get(remove.type) %in% data.dt[,.SD[,.N],by = .(remove.type = get(remove.type), condition)][,.SD[,.N, by = .(remove.type)]][N>1,remove.type]]
  return(data.dt)
}


#FROM THOMPSON, DAVEY ET AL 2018

f.effect.TD = function(data.dt, cluster.dt){

  #Summarise for each condition, period, and cluster
  cpsummary.dt = data.dt[,.(summary = mean(outcome)), by = c('condition','period','cluster')]

  #Aggregate data by time and intervention condition
  slices.long = aggregate(summary ~ period + condition, data = cpsummary.dt,
                           FUN = function(X) c(mean = mean(X, na.rm = TRUE),
                                               n = sum(!is.na(X)),
                                               var = var(X, na.rm = TRUE)))

  slices.long$mean = slices.long$summary[,1]
  slices.long$n = slices.long$summary[,2]
  slices.long$var = slices.long$summary[,3]

  #reshape so one row per time slice
  slices = reshape(direction = "wide",
                    data = slices.long[,c("period", "condition", "mean", "n", "var")],
                    v.names = c("mean", "n", "var"),
                    idvar = c("period"),
                    timevar = "condition")

  #calculate difference
  slices$diff = slices$mean.intervention - slices$mean.control

  #Calculate a weight assuming the same variance in both arms
  slices$wgt =
    ((((slices$n.control - 1) * slices$var.control +
         (slices$n.intervention - 1) *  slices$var.intervention) /
        (slices$n.control + slices$n.intervention - 2)) * (1 / slices$n.control + 1 /slices$n.intervention)) ^ -1

  weighted.mean(slices$diff, slices$wgt, na.rm = TRUE)
}

#' Calculates intracluster correlation coefficient for observations
#'

#' @param data.dt data.table with minimally columns for cluster and outcome.
#' @param mode Type of variable for which ICC will be calculated. Can be
#'   automatically detected on the basis of number of levels. (Default: auto)
#' @param return.type Can return either the single value rho, or a list of all
#'   values calculated (Default: 'value')
#' @return Numeric ICC
#' @export
calculate.icc = function(data.dt,
                         mode = c('binary','continuous','auto')[3],
                         return.type = c('value', 'variables')[1]) {

  #Automatically detect variable type.
  if (mode == 'auto') {
    if (length(unique(data.dt[,outcome])) <= 2) {
      mode = 'binary'
    } else {
      mode = 'continuous'
    }
  }


  if (mode == 'continuous') {
    summary_aov = summary(aov(outcome ~ cluster,data=data.dt))
    ssw = summary_aov[[1]][1,2]
    ssb = summary_aov[[1]][2,2]
    rho = ssw/(ssw + ssb)
    if (return.type == 'variables') {
      return(list(rho = rho,
                  ssw = ssw,
                  ssb = ssb))
    } else {
      return(rho)
    }
    return(summary_aov[[1]][1,2]/sum(summary_aov[[1]][,2]))
  } else if (mode == 'binary') {

    #From iccbin package

    # Number off clusters
    k = length(unique(data.dt[,cluster]))
    # Number of observations in each cluster
    ni = as.vector(table(data.dt[,cluster]))
    # Total number of observations
    N = sum(ni)

    n0 = (1/(k - 1))*(N - sum((ni^2)/N))
    yi = aggregate(data.dt[,outcome], by = list(data.dt[,cluster]), sum)[ , 2]
    yisq = yi^2
    msb = (1/(k - 1))*(sum(yisq/ni) - (1/N)*(sum(yi))^2)
    msw = (1/(N - k))*(sum(yi) - sum(yisq/ni))
    rho = (msb - msw)/(msb + (n0 - 1)*msw)

    if (return.type == 'variables') {
      return(list(rho = rho,
                  k = k,
                  ni = list(ni),
                  N = N,
                  n0 = n0,
                  yi = yi,
                  yisq = yisq,
                  msb = msb,
                  msw = msw))
    } else {
      return(rho)
    }

  }
}

#' Takes a linear regression model and provided values of predictors and
#' outcome, and calculates a provided (or missing) predictor (or outcome).
#' Probably doesn't work for models with interaction terms.
#'
#' @param model A lm generated model.
#' @param constant.values.list List of provided values for outcome and
#'   predictors.
#' @param predict.var Character name of the variable to predict. If NULL the
#'   function will look for a missing item in outcome.values in
#'   constant.values.list. If predict.var exists in constant.values.list, its
#'   value there will be ignored (Default : NULL)
#' @return A scatter ggplot
#'
#' @export
calculate.predictor.value.from.regression.model = function(model,
                                                           constant.values.list,
                                                           predict.var = NULL) {

  #Create data.table of model coefficients, get outcome name, rearrange the
  #equation to add to zero
  coef.dt = data.table(var.name = c(names(model$coefficients), colnames(model$model[1])),
                       coef = c(model$coefficients, -1))


  #Intercept is always there (multiple by 1)
  constant.values.list = unlist(c(constant.values.list, '(Intercept)' = 1))

  #Auto detect predictor variables if requested.
  if (is.null(predict.var)) {
    predict.var = coef.dt[!var.name %in% names(constant.values.list), var.name]
  }
  #If none are detected, stop.
  if (length(predict.var) == 0) {
    stop('No variables detected for predicting.')
  }


  #Check that predictor variables exist in the model.
  if (sum(predict.var %in% coef.dt[,var.name]) == 0) {
    stop(paste('No ',
               paste(predict.var, collapse = ', '),
               'in model (only ',
               paste(coef.dt$var.name, collapse = ', '),
               ').'))
  }

  #Cannot predict more than one variable
  if (length(predict.var) > 1) {
    stop(paste0('Cannot predict more than one variable value (requested: ',
                paste0(predict.var, collapse = ', '),
                ')'))
  }

  #Make a similar data.table for coefficients (excluding those to be predicted.)
  constant.dt = data.table(var.name = names(constant.values.list),
                           var.value = constant.values.list)


  #Merge the two tables, and get product of coefficients and variable values.
  coef.dt = merge(coef.dt, constant.dt, by = 'var.name', all = T)
  coef.dt[,value := coef * var.value]

  #Calculate predicted value of target predictor term, and divide by predictor coefficient.
  predicted.value = -sum(coef.dt[!var.name %in% predict.var]$value, na.rm = T) / coef.dt[var.name == predict.var, coef]
  names(predicted.value) = predict.var

  return(predicted.value)
}

