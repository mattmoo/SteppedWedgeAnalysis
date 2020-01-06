
#' Update a progress bar.
#'
#' @param progress.bar A text progress bar item
#' @param index Numeric index of the permutation you are in (out of max.r)
#' @param max.r Numeric how many you are expecting total
#' @param max.n.updates Maximum number of updates, there is a bit of overhead,
#'   so don't want it updating every iteration. (Default: 50)
update.progress.bar = function(progress.bar, 
                               index, 
                               max.r, 
                               max.n.updates = 50) {
  
  progress.proportion = index/max.r
  if ((progress.proportion %% (1/max.n.updates)) <= 1/max.r) {
    setTxtProgressBar(progress.bar, progress.proportion)
  }
}

#' Calculates and outputs a message of elapsed time.
#'
#' @param start.time Numeric start time from Sys.time()
#' @param n.iterations Integer how many iterations (time per iteration won't be
#'   calculated if NULL, Default: NULL)
#' @param progress.bar Progress bar to close (ignored if NULL, Default: NULL)
elapsed.time.message = function(start.time, n.iterations = NULL, progress.bar = NULL) {
  
  #Close progress bar if requested.
  if (!is.null(progress.bar)) {
    close(progress.bar)
    # message("")
  }
  elapsed.time = Sys.time()-start.time
  message(paste("Time elapsed: "),hms::as.hms(elapsed.time))
  #Calculate time per iteration if requested.
  if (!is.null(n.iterations)) {
    message(paste("Time per iteration: "),hms::as.hms(elapsed.time/n.iterations))
  }
  message("")
}