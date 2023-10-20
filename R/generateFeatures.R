#' @export
getAggregate <- function(CF, aggSize, aggSlot = "auto"){
  if (aggSlot == "auto"){
    if ("ref_data" %in% names(CF)){
      agg <- as.matrix((dplyr::bind_rows(CF$ref_data)))
      # Perform extra sampling in case extra data has been added
      set.seed(42)
      agg <- agg[sample(c(1:nrow(agg), aggSize)),]
    }
    else if ("ref_paths" %in% names(CF)){
      CF <- addReferencedata(CF, CF$ref_paths, read = TRUE, reload = TRUE, 
                             aggSize = aggSize)
      agg <- as.matrix((dplyr::bind_rows(CF$ref_data)))
    }
    else if ("test_data" %in% names(CF)){
      agg <- as.matrix((dplyr::bind_rows(CF$test_data)))
      set.seed(42)
      agg <- agg[sample(c(1:nrow(agg), aggSize)),]
    }
    else if ("test_paths" %in% names(CF)){
      CF <- addTestdata(CF, CF$test_paths, read = TRUE, reload = TRUE, 
                        aggSize = aggSize)
      agg <- as.matrix((dplyr::bind_rows(CF$test_data)))
    }
  }
  else {
    agg <- as.matrix((dplyr::bind_rows(CF[[paste0(aggSlot, "_data")]])))
    set.seed(42)
    agg <- agg[sample(c(1:nrow(agg), aggSize)),]
  }
  return(agg)
}


#' Generate features and attach to CytoFlag object
#'
#' @param CF CytoFlag object
#' @param channels Channels to calculate features for
#' @param featMethod Feature generation method to use
#' @param n How many cells to use for feature generation
#' @param aggSize How many cells to use in total for aggregate sample (default = 10000)
#' @param aggSlot Whether to use reference or test data as aggregate (default = "auto")
#' @param cores How many cores to use for parallelization (default = 50%)
#' @param recalculate Whether to recalculate features for existing data
#' @param nRecursions Number of recursions to use for fingerprinting (default = 4)
#'
#' @return Dataframe with features (columns) for samples (rows)
#' 
#' @export
generateFeatures <- function(CF, channels, featMethod = "summary", n = 1000,
                              aggSize = 10000, aggSlot = "auto", cores = "auto", 
                              recalculate = FALSE, nRecursions = 4){
  # Determine the number of cores to use
  if (cores == "auto"){
      cores = parallel::detectCores() / 2 
      message(paste("Using 50% of cores:", cores))
  }
  
  # Summary statistics and landmarks can be generated for individual files
  if (featMethod %in% c("summary", "landmarks", "quantiles")){
    if (featMethod == "summary"){
      func <- summaryStats
    }
    if (featMethod == "landmarks"){
      func <- LandmarkStats
    }
    if (featMethod == "quantiles"){
      func <- Quantiles
    }
    # Generate features for the ref and test paths slots (if in CF object)
    for (slot in c("ref", "test")){
      if (paste0(slot, "_paths") %in% names(CF)){
        if (!featMethod %in% names(CF$features[[slot]]) || recalculate == TRUE){
          # Generate features for all paths
          CF$features[[slot]][[featMethod]] <- func(CF, CF[[paste0(slot, "_paths")]],
                                                    channels, n, cores)
        }
        else {
          # Generate features for the new file paths
          new_paths <- c()
          for (path in CF[[paste0(slot, "_paths")]]){
            if (!path %in% rownames(CF$features[[slot]][[featMethod]])){
              message("Generating additional features")
              new_paths <- c(new_paths, path)
            }
          }
          if (length(new_paths) == 0){
            message(paste("Did not detect any new files for", slot, "slot"))
            message("Force re-calculation of statistics in this slot using recalculate = TRUE.")
            next
          } 
          else {
            new_stats <- func(CF, new_paths, channels, n, cores)
          }
          CF$features[[slot]][[featMethod]] <- rbind(CF$features[[slot]][[featMethod]], 
                                                     new_stats)
        }
      }
    }
    return(CF)
  }
  
  # Everything from here depends on aggregated data
  agg <- getAggregate(CF, aggSize = aggSize, aggSlot = aggSlot)

  for (slot in c("ref", "test")){
    if (paste0(slot, "_paths") %in% names(CF)){
      if (!featMethod %in% names(CF$features[[slot]]) || recalculate == TRUE){
        paths <- CF[[paste0(slot, "_paths")]]
        if (featMethod == "binning"){
          CF$features[[slot]][[featMethod]] <- Bin(CF, paths, agg, channels, n, 
                                                   cores)
        }
        else if (featMethod == "EMD"){
          CF$features[[slot]][[featMethod]] <- EMD(CF, paths, agg, channels, n, 
                                                   cores)
        }
        else if (featMethod == "fingerprint"){
          CF$features[[slot]][[featMethod]] <- Fingerprint(CF, paths, agg, 
                                                           channels, n, cores, 
                                                           nRecursions)
        }
      }
    }
  }
  return(CF)
}

#' @export
calculateQuantiles <- function(CF, path, channels, n){
  ff <- readInput(CF, path, n)
  df <- data.frame(ff@exprs[,channels], check.names = FALSE)
  stats <- list()
  # Calculate quantiles for every variable
  percentiles <- lapply(df, function(col) quantile(col,
                                                   probs = seq(0.1, 0.9, 0.1)))
  stats <- list()
  for (i in names(percentiles)){
    for (j in names(percentiles[[i]])){
      stats[[paste0(i, "_", j)]] <- percentiles[[i]][[j]]
    }
  }
  stats <- data.frame(stats, check.names=FALSE)
  return(stats)
}

#' @export
Quantiles <- function(CF, input, channels, n, cores){
  if (cores > 1){
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    parallel::clusterExport(cl, c(CF[["parallel"]][["parallelVars"]], "calculateQuantiles"))
    `%dopar%` <- foreach::`%dopar%`
    all_stats <- foreach::foreach(path = input, .combine = "c", 
                                  .packages = CF[["parallel"]][["parallelPackages"]]) %dopar% {
                                    stats <- list(calculateQuantiles(CF, path, channels, n))
                                    names(stats) <- path
                                    return(stats)
                                  }
    parallel::stopCluster(cl)
  }
  else {
    all_stats <- list()
    for (path in input){
      stats <- calculateQuantiles(CF, path, channels, n)
      all_stats[[path]] <- stats
    }
  }
  stats <- data.frame(dplyr::bind_rows(all_stats), check.names = FALSE)
  rownames(stats) <- names(all_stats)
  return(stats)
}

#' @export
calculateSummary <- function(CF, path, channels, n){
  ff <- readInput(CF, path, n)
  stats <- list()
  for (channel in channels){
    stats[paste0(channel,"_mean")] <- base::mean(ff@exprs[,channel])
    stats[paste0(channel,"_sd")] <- stats::sd(ff@exprs[,channel])
    stats[paste0(channel,"_median")] <- stats::median(ff@exprs[,channel])
    stats[paste0(channel,"_IQR")] <- stats::IQR(ff@exprs[,channel])
  }
  return(stats)
}

#' @export
summaryStats <- function(CF, input, channels, n, cores){
  if (cores > 1){
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    parallel::clusterExport(cl, c(CF[["parallel"]][["parallelVars"]], "calculateSummary"))
    `%dopar%` <- foreach::`%dopar%`
    all_stats <- foreach::foreach(path = input, .combine = "c", 
                           .packages = CF[["parallel"]][["parallelPackages"]]) %dopar% {
                           stats <- list(calculateSummary(CF, path, channels, n))
                           names(stats) <- path
                           return(stats)
                         }
    parallel::stopCluster(cl)
  }
  else {
    all_stats <- list()
    for (path in input){
      stats <- calculateSummary(CF, path, channels, n)
      all_stats[[path]] <- stats
    }
  }
  stats <- data.frame(dplyr::bind_rows(all_stats), check.names = FALSE)
  rownames(stats) <- names(all_stats)
  return(stats)
}


#' @export
calculateEMD <- function(CF, path, agg, channels, n){
  ff <- readInput(CF, path, n = n)
  stats <- list()
  for (channel in channels){
    stats[paste0(channel,'_', 'EMD')] <- transport::wasserstein1d(ff@exprs[, channel], 
                                                                  agg[, channel])
  }
  return(stats)
}


#' @export
EMD <- function(CF, input, agg, channels, n, cores){
  agg <- agg
  if (cores > 1){
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    parallel::clusterExport(cl, c(CF[["parallel"]][["parallelVars"]], "agg", "calculateEMD"),
                            envir=environment())
    `%dopar%` <- foreach::`%dopar%`
    all_stats <- foreach::foreach(path = input, .combine = "c", 
                                  .packages = c(CF[["parallel"]][["parallelPackages"]], "transport")) %dopar% {
                                    stats <- list(calculateEMD(CF, path, agg, channels, n))
                                    names(stats) <- path
                                    return(stats)
                                  }
    parallel::stopCluster(cl)
  }
  else {
    all_stats <- list()
    for (path in input){
      all_stats[[path]] <- calculateEMD(CF, path, agg, channels, n)
    }
  }
  stats <- data.frame(dplyr::bind_rows(all_stats), check.names = FALSE)
  rownames(stats) <- names(all_stats)
  return(stats)
}


#' @export
calculateFingerprint <- function(CF, path, model, channels, n){
  ff <- readInput(CF, path, n = n)
  call = flowFP::flowFP(ff[, channels], model)
  bin_counts = flowFP::counts(call)
  # Convert counts to bin frequencies
  stats = data.frame(t(apply(bin_counts, 1, function(x) x/sum(x))))
  stats = list(stats)
  return(stats)
}


#' @export
Fingerprint <- function(CF, input, agg, channels, n, cores, nRecursions = 4){
  # Convert aggregated matrix to flowframe
  agg <- flowCore::flowFrame(agg[,channels])
  model <- flowFP::flowFPModel(agg, parameters = channels, 
                               nRecursions = nRecursions)
  if (cores > 1){
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    parallel::clusterExport(cl, c(CF[["parallel"]][["parallelVars"]], "model", 
                                  "calculateFingerprint"), envir=environment())
    `%dopar%` <- foreach::`%dopar%`
    all_stats <- foreach::foreach(path = input, .combine = "c", 
                                  .packages = c(CF[["parallel"]][["parallelPackages"]], "flowFP")) %dopar% {
                                  stats <- calculateFingerprint(CF, path, model, channels, n)
                                  names(stats) <- path
                                  return(stats)
                                  }
    parallel::stopCluster(cl)
  } else {
    all_stats <- list()
    for (path in input){
      all_stats[[path]] <- calculateFingerprint(CF, path, model, channels, n)
    }
  }
  stats <- data.frame(dplyr::bind_rows(all_stats), check.names = FALSE)
  rownames(stats) <- names(all_stats)
  return(stats)
}


#' @export
calculateBins <- function(CF, path, bin_boundaries, channels, n){
  ff <- readInput(CF, path, n = n)
  df <- data.frame(ff@exprs[,channels], check.names=FALSE)
  stats <- c()
  for (i in seq_along(channels)){
    # Calculate the frequencies per bin
    counts <- as.numeric(t(data.frame(table(cut(df[, channels[i]], 
                                                breaks = bin_boundaries[, i]))))[2,])
    freqs <- counts / sum(counts)
    names(freqs) <- paste(channels[i], "_bin", 1:10, sep = "")
    stats <- c(stats, freqs)
  }
  stats <- list(stats)
  return(stats)
}


#' @export
Bin <- function(CF, input, agg, channels, n, cores){
  # Determine the bins on the aggregated data
  bin_boundaries <- apply(agg[,channels], 2, 
                          function(x) stats::quantile(x, probs = seq(0, 1, by = 0.1)))
  bin_boundaries <- data.frame(bin_boundaries)
  if (cores > 1){
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    parallel::clusterExport(cl, c(CF[["parallel"]][["parallelVars"]], "bin_boundaries", 
                                  "calculateBins"), envir=environment())
    `%dopar%` <- foreach::`%dopar%`
    all_stats <- foreach::foreach(path = input, .combine = "c", 
                                  .packages = CF[["parallel"]][["parallelPackages"]]) %dopar% {
                                    stats <- calculateBins(CF, path, bin_boundaries, 
                                                           channels, n)
                                    names(stats) <- path
                                    return(stats)
                                  }
    parallel::stopCluster(cl)
  }
  else {
    all_stats <- list()
    for (path in input){
      all_stats[[path]] <- calculateBins(CF, path, bin_boundaries, channels, n)
    }
  }
  stats <- data.frame(dplyr::bind_rows(all_stats), check.names = FALSE)
  rownames(stats) <- names(all_stats)
  return(stats)
}
