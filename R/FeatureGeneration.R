#' Generate features and attach to CytoFlag object
#'
#' @param CF CytoFlag object
#' @param channels Channels to calculate features for
#' @param featMethod Feature generation method to use.
#' Either one of: "summary", "EMD", "binning" or "fingerprint".
#' @param nRecursions Number of recursions to use for fingerprinting (default = 4)
#'
#' @return Dataframe with features (columns) for samples (rows)
#' 
#' @export
FeatureGeneration <- function(CF, channels, featMethod = "summary", 
                              nRecursions = 4, cores = "auto"){
  # Determine the number of cores to use
  if (cores == "auto"){
      cores = detectCores() / 2 
      message(paste("Using 50% of cores:", cores))
  }
  
  if (!"test_data" %in% names(CF)){
    if (!"test_paths" %in% names(CF)){
      message("No test data or paths found in CytoFlag object")
    }
  }
  
  if (featMethod == "summary"){
    # TO-DO: CHECK IF NOT EMPTY
    # For summary statistics, use the paths directly
    if ("ref_paths" %in% names(CF)){
      CF$features$ref$summary <- SummaryStats(CF, CF$ref_paths, channels, cores)
    }
    if ("test_paths" %in% names(CF)){
      CF$features$test$summary <- SummaryStats(CF, CF$test_paths, channels, cores)
    }
    return(CF)
  }
  
  if (featMethod == "landmarks"){
    # TO-DO:
    # CHECK IF LANDMARKS HAVE BEEN GENERATED FOR ALL CHANNELS
    if ("ref_data" %in% names(CF)){
      CF$features$ref$landmarks <- LandmarkStats(CF, CF$ref_paths, channels,
                                                 cores)
      CF$features$test$landmarks <- LandmarkStats(CF, CF$test_paths, channels,
                                                  cores)
    }
    else if ("test_data" %in% names(CF)){
      CF$features$test$landmarks <- LandmarkStats(CF, CF$test_paths, channels,
                                                  cores)
    }
  }
  
  # Everything from here on depends on aggregated data for features
  if ("ref_data" %in% names(CF)){
    message("Using aggregated reference data to calculate features")
    agg <- as.matrix((dplyr::bind_rows(CF$ref_data)))
  }
  else if ("ref_paths" %in% names(CF)){
    message("Aggregating reference data to calculate features")
    CF <- AddReferenceData(CF, CF$ref_paths, read = TRUE)
    message("Using aggregated reference data to calculate features")
    agg <- as.matrix((dplyr::bind_rows(CF$ref_data)))
  }
  else if ("test_data" %in% names(CF)){
    message("Using aggregated test data to calculate features")
    agg <- as.matrix((dplyr::bind_rows(CF$test_data)))
  }
  else if ("test_paths" %in% names(CF)){
    message("Aggregating test data to calculate features")
    CF <- AddTestData(CF, CF$test_paths, read = TRUE)
    message("Using aggregated test data to calculate features")
    agg <- as.matrix((dplyr::bind_rows(CF$test_data)))
  }
  
  if (featMethod == "binning"){
    if ("ref_data" %in% names(CF)){
      CF$features$ref$binning <- Bin(CF, CF$ref_paths, agg, channels)
      CF$features$test$binning <- Bin(CF, CF$test_paths, agg, channels)
    }
    else if ("test_data" %in% names(CF)){
      CF$features$test$binning <- Bin(CF, CF$test_paths, agg, channels)
    }
  }
  
  if (featMethod == "EMD"){
    if ("ref_data" %in% names(CF)){
      CF$features$ref$EMD <- EMD(CF, CF$ref_paths, agg, channels, cores)
      CF$features$test$EMD <- EMD(CF, CF$test_paths, agg, channels, cores)
    }
    else if ("test_data" %in% names(CF)){
      CF$features$test$EMD <- EMD(CF, CF$test_paths, agg, channels, cores)
    }
  }
  
  if (featMethod == "fingerprint"){
    if ("ref_data" %in% names(CF)){
      CF$features$ref$fingerprint <- Fingerprint(CF, CF$ref_paths, agg, channels, 
                                                 cores, nRecursions)
      CF$features$test$fingerprint <- Fingerprint(CF, CF$test_paths, agg, channels, 
                                                  cores, nRecursions)
    }
    else if ("test_data" %in% names(CF)){
      CF$features$test$fingerprint <- Fingerprint(CF, CF$test_paths, agg, channels, 
                                                  cores, nRecursions)
    }
  }
  return(CF)
}


calculateSummary <- function(CF, path, channels, n = 1000){
  ff <- ReadInput(CF, path, n = n)
  stats <- list()
  for (channel in channels){
    stats[paste0(channel,"_mean")] <- base::mean(ff@exprs[,channel])
    stats[paste0(channel,"_sd")] <- stats::sd(ff@exprs[,channel])
    stats[paste0(channel,"_median")] <- stats::median(ff@exprs[,channel])
    stats[paste0(channel,"_IQR")] <- stats::IQR(ff@exprs[,channel])
  }
  return(stats)
}


SummaryStats <- function(CF, input, channels, cores){
  if (cores > 1){
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    parallel::clusterExport(cl, c("channels", "CF", "ReadInput", "calculateSummary"))
    all_stats <- foreach::foreach(path = input, .combine = "c", 
                           .packages=c("flowCore","PeacoQC")) %dopar% {
                           stats <- list(calculateSummary(CF, path, channels))
                           names(stats) <- path
                           return(stats)
                         }
    parallel::stopCluster(cl)
  }
  else {
    all_stats <- list()
    for (path in input){
      stats <- calculateSummary(CF, path, channels)
      all_stats[[path]] <- stats
    }
  }
  stats <- data.frame(dplyr::bind_rows(all_stats), check.names = FALSE)
  rownames(stats) <- names(all_stats)
  return(stats)
}


calculateLandmarks <- function(CF, path, channels, n = 1000){
  ff <- ReadInput(CF, path, n = n)
  df <- data.frame(ff@exprs[,channels], check.names = FALSE)
  stats <- list()
  for (channel in channels){
    agg_peaks <- CF$landmarks$ref[[channel]]$landmarks
    for (i in seq(1, nrow(agg_peaks))){
      # Sample all the values in the expression range
      min <- agg_peaks[i, "exprs_start"]
      max <- agg_peaks[i, "exprs_end"]
      peak_exprs <- df[df[,channel] > min & df[,channel] < max, ]
      if (is.null(peak_exprs)){
        print("too little cells")
        stats[paste0(channel, "_peak", i, "_median")] <- NA
      }
      else if (nrow(peak_exprs) < 20){
        stats[paste0(channel, "_peak", i, "_median")] <- NA
      }
      else {
        stats[paste0(channel, "_peak", i, "_median")] <- stats::median(peak_exprs[,channel])
      }
    }
  }
  return(stats)
}


LandmarkStats <- function(CF, input, channels, cores){
  if (cores > 1){
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    parallel::clusterExport(cl, c("channels", "CF", "ReadInput", "calculateLandmarks"))
    all_stats <- foreach::foreach(path = input, .combine = "c", 
                                  .packages=c("flowCore","PeacoQC")) %dopar% {
                                  stats <- list(calculateLandmarks(CF, path, channels))
                                  names(stats) <- path
                                  return(stats)
                                  }
    parallel::stopCluster(cl)
  } else {
    all_stats <- list()
    for (path in input){
      all_stats[[path]] <- calculateLandmarks(CF, path, channels)
    }
  }
  stats <- data.frame(dplyr::bind_rows(all_stats), check.names = FALSE)
  # Impute missing values
  stats <- VIM::kNN(stats, k = 1, imp_var	= FALSE)
  rownames(stats) <- names(all_stats)
  return(stats)
}


calculateEMD <- function(CF, path, agg, channels, n = 1000){
  ff <- ReadInput(CF, path, n = n)
  stats <- list()
  for (channel in channels){
    stats[paste0(channel,'_', 'EMD')] <- transport::wasserstein1d(ff@exprs[, channel], 
                                                                  agg[, channel])
  }
  return(stats)
}


EMD <- function(CF, input, agg, channels, cores){
  agg <- agg
  if (cores > 1){
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    parallel::clusterExport(cl, c("channels", "CF", "agg", "ReadInput", "calculateEMD"),
                            envir=environment())
    all_stats <- foreach::foreach(path = input, .combine = "c", 
                                  .packages=c("flowCore","PeacoQC", "transport")) %dopar% {
                                    stats <- list(calculateEMD(CF, path, agg, channels))
                                    names(stats) <- path
                                    return(stats)
                                  }
    parallel::stopCluster(cl)
  }
  else {
    all_stats <- list()
    for (path in input){
      all_stats[[path]] <- calculateEMD(CF, path, agg, channels)
    }
  }
  stats <- data.frame(dplyr::bind_rows(all_stats), check.names = FALSE)
  rownames(stats) <- names(all_stats)
  return(stats)
}


calculateFingerprint <- function(CF, path, model, channels, n = 1000){
  ff <- ReadInput(CF, path, n = n)
  call = flowFP::flowFP(ff[, channels], model)
  bin_counts = flowFP::counts(call)
  # Convert counts to bin frequencies
  stats = data.frame(t(apply(bin_counts, 1, function(x) x/sum(x))))
  stats = list(stats)
  return(stats)
}


Fingerprint <- function(CF, input, agg, channels, cores, nRecursions = 4){
  # Convert aggregated matrix to flowframe
  agg <- flowCore::flowFrame(agg[,channels])
  message("Fitting flowFP fingerprinting model")
  model <- flowFP::flowFPModel(agg, parameters = channels, 
                               nRecursions = nRecursions)
  if (cores > 1){
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    parallel::clusterExport(cl, c("channels", "CF", "model", "ReadInput", 
                                  "calculateFingerprint"),
                            envir=environment())
    all_stats <- foreach::foreach(path = input, .combine = "c", 
                                  .packages=c("flowCore","PeacoQC", "flowFP")) %dopar% {
                                  stats <- calculateFingerprint(CF, path, model, channels)
                                  names(stats) <- path
                                  return(stats)
                                  }
    parallel::stopCluster(cl)
  } else {
    all_stats <- list()
    for (path in input){
      all_stats[[path]] <- calculateFingerprint(CF, path, model, channels)
    }
  }
  stats <- data.frame(dplyr::bind_rows(all_stats), check.names = FALSE)
  rownames(stats) <- names(all_stats)
  return(stats)
}


Bin <- function(CF, input, agg, channels){
  # Determine the bins on the aggregated data
  message("Determining bin boundaries on aggregated data")
  bin_boundaries <- apply(agg[,channels], 2, 
                          function(x) stats::quantile(x, probs = seq(0, 1, by = 0.1)))
  bin_boundaries <- data.frame(bin_boundaries)
  
  all_stats <- list()
  for (path in input){
    ff <- ReadInput(CF, path, n = NULL)
    df <- data.frame(ff@exprs[,channels], check.names=FALSE)
    stats <- c()
    for (i in seq_along(channels)){
      # Calculate the frequencies per bin
      counts <- as.numeric(t(data.frame(table(cut(df[, channels[i]], breaks = bin_boundaries[, i]))))[2,])
      freqs <- counts / sum(counts)
      names(freqs) <- paste(channels[i], "_bin", 1:10, sep = "")
      stats <- c(stats, freqs)
    }
    all_stats[[path]] <- stats
  }
  stats <- data.frame(dplyr::bind_rows(all_stats), check.names = FALSE)
  return(stats)
}


#' Reduce features to first two principal components
#'
#' @param testFeatures Dataframe of features from test samples
#' @param flagStrat Which anomaly detection strategy to use ("outlier" or "novelty").
#' @param refFeatures Dataframe of features from reference samples
#'
#' @return List containing PC loadings, variance and reduced features
reduceDim <- function(testFeatures, flagStrat, refFeatures = NULL){
  if (flagStrat == "outlier"){
    pca <- stats::prcomp(testFeatures, center = TRUE, scale = TRUE)
    testFeatures$PC1 <- pca$x[, 1]
    testFeatures$PC2 <- pca$x[, 2]
    testFeatures <- testFeatures[c("PC1", "PC2")]
    return(list("loadings" = as.data.frame(pca$rotation[, 1:2]),
                "PC1_var" = pca$sdev[1]^2,
                "PC2_var" = pca$sdev[2]^2,
                "testFeatures" = testFeatures))
  }
  if (flagStrat == "novelty"){
    # Fit the PCA on the reference features before reducing the test features
    pca <- stats::prcomp(refFeatures, center = TRUE, scale = TRUE)
    test_PCA <- stats::predict(pca, newdata = testFeatures)
    testFeatures$PC1 <- test_PCA[,"PC1"]
    testFeatures$PC2 <- test_PCA[,"PC2"]
    testFeatures <- testFeatures[c("PC1", "PC2")]
    ref_PCA <- stats::predict(pca, newdata = refFeatures)
    refFeatures$PC1 <- ref_PCA[,"PC1"]
    refFeatures$PC2 <- ref_PCA[,"PC2"]
    refFeatures <- refFeatures[c("PC1", "PC2")]
    return(list("loadings" = as.data.frame(pca$rotation[, 1:2]),
                "PC1_var" = pca$sdev[1]^2,
                "PC2_var" = pca$sdev[2]^2,
                "testFeatures" = testFeatures,
                "refFeatures" = refFeatures))
  }
}