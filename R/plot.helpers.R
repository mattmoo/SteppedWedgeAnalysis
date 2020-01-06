#' Plots a line and scatter plot faceted by sequence and cluster with a few
#' visual aids. If it's a permutation, maybe provide the actual clusters so that
#' the theming can work properly. sort the input by outcome to make it faster.
#'
#' @param data.dt data.table with columns participant, cluster, time, and
#'   outcome. Outcome should be continuous.
#' @param sequence.dt data.table with information about the sequences, with
#'   columns sequence, transition.time, and intervention.time.
#' @param cluster.dt data.table with the correspondence between cluster and
#'   sequence, with columns cluster and sequence.
#' @param actual.cluster.dt This function is partly designed to visualise
#'   permutation statistics, if cluster.dt is permuted, then you might provide a
#'   similar data.table with the actual correpsondence between cluster and
#'   sequence so that the graph can be themed as such (i.e. data colours
#'   matching original sequence).
#' @param min.time Numeric minimum time on x axis (Default: 0)
#' @param max.time Numeric maximum time on x axis, will try to automatically
#'   generate from sequence.dt if NULL (Default: NULL)
#' @return A ggplot with data faceted by sequence and cluster.
#'
#' @export
faceted.line.plot = function(data.dt,
                             sequence.dt,
                             cluster.dt,
                             actual.cluster.dt = NULL,
                             min.time = 0,
                             max.time = NULL) {
  
  #You might need the actual clusters for visualising when they have been permuted.
  if (!is.null(actual.cluster.dt)) {
    actual.cluster.dt = copy(actual.cluster.dt)
  } else {
    actual.cluster.dt = copy(cluster.dt)
  }
  setnames(actual.cluster.dt, old = 'cluster', new = 'actual.cluster')
  
  #Estimate max.time if it is not provided.
  if(is.null(max.time)) {
    max.time = 2 * sequence.dt[,max(intervention.time, na.rm = T)] - 
      sequence.dt[intervention.time != max(intervention.time, na.rm = T)][
        ,max(intervention.time, na.rm = T)] 
  }
  
  #Someone might have already attached sequence and cluster to data.dt, if so,
  #remove them.
  scatter.dt = data.dt[,c(names(data.dt)[!names(data.dt) %in% c(names(sequence.dt),names(cluster.dt))],
                       'cluster'), with = F]
  
  #Hack so that there is no malfunctioning when a cluster has no intervention.
  # sequence.dt = copy(sequence.dt)[is.na(transition.time), transition.time := max(data.dt[,time]*1000) ]
  # sequence.dt = copy(sequence.dt)[is.na(intervention.time), intervention.time := max(data.dt[,time]*1000) ]
  
  #Data for scatter plots
  #Associate participant and cluster
  scatter.dt = merge(scatter.dt, cluster.dt, by = 'cluster')
  #Associate cluster and sequence
  scatter.dt = merge(scatter.dt, sequence.dt, by = 'sequence')
  #Attach original sequence data
  scatter.dt = merge(scatter.dt,
                     merge(
                       merge(
                         scatter.dt[, .(cluster, participant, time)],
                         actual.cluster.dt,
                         by.x = 'cluster',
                         by.y = 'actual.cluster'
                       ),
                       sequence.dt[, .(sequence, intervention.time)],
                       by = 'sequence'
                     )[, .(participant = participant,
                           actual.intervention = time >= intervention.time,
                           actual.sequence = sequence)],
                     by = 'participant')
  scatter.dt[is.na(actual.intervention), actual.intervention := F]
  
  
  #Sequence-specific intervention time data
  vline.dt = merge(cluster.dt,sequence.dt, by = 'sequence')

  #Remove clusters without intervention from rect/line plots.
  vline.dt = vline.dt[!is.na(intervention.time)]
  #If there is not transition.time, assume there was no transition.
  vline.dt[is.na(transition.time), transition.time := intervention.time]
  
  #Generate all period time data
  period.dt = data.table(expand.grid(sequence = sequence.dt[,sequence], intervention.time = sequence.dt[,intervention.time]))
  period.dt = merge(period.dt, vline.dt[,.(sequence, cluster)], by = 'sequence', allow.cartesian = T)
  period.dt = period.dt[!is.na(intervention.time)]
  
  
  output.plot = ggplot2::ggplot() +
    #Background colours marking period
    ggplot2::geom_rect(data = vline.dt,
                       aes(fill = sequence,
                           xmin = -Inf, 
                           xmax = transition.time, 
                           ymin = -Inf, 
                           ymax = Inf),
                       alpha = .05) +
    ggplot2::geom_rect(data = vline.dt,
                       aes(fill = sequence,
                           xmin = transition.time, 
                           xmax = intervention.time, 
                           ymin = -Inf, 
                           ymax = Inf),
                       alpha = .2) +
    ggplot2::geom_rect(data = vline.dt,
                       aes(fill = sequence,
                           xmin = intervention.time, 
                           xmax = Inf, 
                           ymin = -Inf, 
                           ymax = Inf),
                       alpha = .4) +
    #Lines marking intervention for that sequence
    ggplot2::geom_vline(data = vline.dt,
                        aes(xintercept = intervention.time,
                            colour = sequence),
                        alpha = .9,
                        size = 1)  +
    #Lines between periods
    ggplot2::geom_vline(data = period.dt,
                        aes(xintercept = intervention.time,
                            colour = sequence),
                        alpha = .5,
                        size = .5)  +
    #New theme so that colours for clusters and sequences can be separated
    scale_fill_brewer(palette="Dark2") +
    new_scale_fill() +
    #Best fit lines
    ggplot2::geom_line(data=scatter.dt,
                       aes(x=time, 
                           y=outcome), 
                       stat = "smooth", 
                       method = 'loess', 
                       formula = 'y ~ x', 
                       alpha = 0.4, 
                       span = 0.2) +
    #Confidence intervals of best fit
    ggplot2::geom_ribbon(data=scatter.dt,
                         aes(x=time, 
                             y=outcome, 
                             fill=actual.sequence), 
                         stat='smooth', 
                         method = 'loess', 
                         se=TRUE, 
                         alpha=0.3, 
                         span = 0.2) +
    #Scatter plot
    ggplot2::geom_point(data=scatter.dt,
                        aes(shape = actual.intervention,
                            fill=actual.sequence, 
                            x=time, 
                            y=outcome), 
                        alpha=.6, 
                        size = 2) +
    #Themes
    scale_fill_brewer(palette="Dark2") +
    scale_colour_brewer(palette="Dark2") +
    scale_shape_manual(values = c(25,24)) +
    ggplot2::theme(legend.position="none",
                   panel.grid.minor = element_blank()) +
    ggplot2::scale_x_continuous(breaks = c(min.time,
                                           sequence.dt[,intervention.time],
                                           max.time),
                                limits = c(min.time, max.time),
                                minor_breaks = NULL)+
    ggplot2::facet_grid(sequence + cluster ~ .)
  
  return(output.plot)
}


#' Creates a scatter plot of ICC vs one or two variables, ICC on y-axis. Chooses
#' automatically between a 3D plotly and 2D ggplot. If only one or two variables
#' in the parameters has multiple levels across scenarios, these will be
#' automatically detected.
#'
#' @param icc.dt A data.table with the parameters of the various simulations and
#'   ICC.
#' @param vars.to.plot String name of the variable that will be on the x-axis.
#'   Will autodetect if NULL, but that will fail if there are more than three
#'   potential variables to plot ICC against. (Default: NULL)
#' @param plot.mean If true, plot the mean value at each point of var.to.plot,
#'   otherwise plot all values individually.
#' @param cols.to.ignore Columns that will be ignored by detection of
#'   vars.to.plot (Default: c('iteration', 'icc'))
#' @return A scatter ggplot
#'
#' @export
icc.simulation.scatter.plot = function(icc.dt, 
                                       vars.to.plot = NULL, 
                                       plot.mean = T, 
                                       cols.to.ignore = c('iteration', 'icc')) {
  
  #Auto detect levels to plot against if requested
  if (is.null(vars.to.plot)) {
    vars.to.plot = cols.exceed.level.threshold(icc.dt,
                                               threshold = 1,
                                               cols.to.ignore = c('iteration', 'icc'))
  }
  
  if (length(vars.to.plot) >= 3) {
    stop(paste0('Autodetection of vars.to.plot failed, more than two potential variables to plot against (',
                vars.to.plot, ').'))
  } else if (length(vars.to.plot) <= 0) {
    stop('Autodetection of vars.to.plot failed, no variables with more than one level.')
  } else if (length(vars.to.plot) == 1) {
    output.plot = icc.simulation.2D.scatter.plot(icc.dt = icc.dt,
                                                 var.to.plot = vars.to.plot,
                                                 plot.mean = plot.mean)
  } else if (length(vars.to.plot) == 2) {
    output.plot = icc.simulation.3D.scatter.plot(icc.dt = icc.dt,
                                                 vars.to.plot = vars.to.plot,
                                                 plot.mean = plot.mean)
  }
  
  
}

#' Creates a scatter plot of ICC vs one variable, ICC on y-axis.
#'
#' @param icc.dt A data.table with the parameters of the various simulations and
#'   ICC.
#' @param var.to.plot String name of the variable that will be on the x-axis.
#' @param plot.mean If true, plot the mean value at each point of var.to.plot,
#'   otherwise plot all values individually.
#' @return A scatter ggplot
#'
#' @export
icc.simulation.2D.scatter.plot = function(icc.dt, var.to.plot, plot.mean = T) {
  
  if (plot.mean) {
    plot.dt = icc.dt[,.(icc = mean(icc)), by = var.to.plot]
  } else {
    plot.dt = icc.dt
  }
  
  plot.dt = icc.bin.dt[get(var.to.plot)==1]
  
  
  output.plot = ggplot(plot.dt, aes_string(x = var.to.plot,
                                           y = icc), 
                       alpha = 0.15) +
    geom_point()
  
  
  return(output.plot)
}

#' Creates a scatter plot of ICC vs two other variables, ICC on y-axis.
#'
#' @param icc.dt A data.table with the parameters of the various simulations and
#'   ICC.
#' @param vars.to.plot Vector of two string names of the variable that will be on
#'   the x- and y-axes.
#' @param plot.mean If true, plot the mean value at each point of var.to.plot,
#'   otherwise plot all values individually.
#' @return A plotly 3D scatter plot
#'
#' @export
icc.simulation.3D.scatter.plot = function(icc.dt, vars.to.plot, plot.mean = T) {
  
  if (plot.mean) {
    plot.dt = icc.dt[,.(icc = mean(icc)), by = vars.to.plot]
  } else {
    plot.dt = icc.dt
  }

  output.plot = plot_ly(
    plot.dt,
    x = ~ get(vars.to.plot[1]),
    y = ~ get(vars.to.plot[2]),
    z = ~ icc,
    alpha = 0.8,
    color = ~ icc
  ) %>%
    layout(
      scene = list(
        xaxis = list(title = vars.to.plot[1]),
        yaxis = list(title = vars.to.plot[2]),
        zaxis = list(title = "ICC"),
        marker = list(size = 1)
      )) %>%
    add_markers()
  
  return(output.plot)
}

#' Creates a 3D scatter plot of ICC vs two other variables, ICC on y-axis.
#'
#' @param icc.model A model for predicting ICC from given variables.
#' @param vars.to.plot Vector of two string names of the variable that will be
#'   on the x- and y-axes. Will autodetect if NULL, but that will fail if there
#'   are more than three potential variables to plot ICC against. (Default:
#'   NULL)
#' @param x.vals x values to generate values for (Default: 1:10 by 0.1)
#' @param y.vals y values to generate values for (Default: 1:10 by 0.1)
#' @return A plotly 3D scatter plot of predicted values from a model.
#'
#' @export
icc.model.3D.scatter.plot = function(icc.model,
                                     vars.to.plot = NULL,
                                     x.vals = seq(0, 10, by = 0.1),
                                     y.vals = seq(0, 10, by = 0.1)) {
  

  plot.dt = data.table(expand.grid(x = x.vals,
                                   y = y.vals))
  setnames(plot.dt, old = c('x','y'), new = vars.to.plot)
  
  plot.dt[, icc := predict(icc.model, plot.dt)]

  output.plot = icc.simulation.3D.scatter.plot(icc.dt = plot.dt,
                                               vars.to.plot = vars.to.plot,
                                               plot.mean = F)

  return(output.plot)
}

