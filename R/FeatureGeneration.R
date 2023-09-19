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
                              nRecursions = 4){
  if (featMethod == "summary"){
    # TO-DO: CHECK IF NOT EMPTY
    # For summary statistics, use the paths directly
    if ("ref_paths" %in% names(CF)){
      CF$features$ref$summary <- SummaryStats(CF$ref_paths, channels)
    }
    if ("test_paths" %in% names(CF)){
      CF$features$test$summary <- SummaryStats(CF$test_paths, channels)
    }
    return(CF)
  }
  
  if (featMethod == "peaks"){
    if ("ref_paths" %in% names(CF)){
      CF$features$ref$peaks <- PeakExtraction(CF$ref_paths, channels)
    }
    if ("test_paths" %in% names(CF)){
      CF$features$test$peaks <- PeakExtraction(CF$test_paths, channels)
    }
    return(CF)
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
      CF$features$ref$binning <- Bin(CF$ref_paths, agg, channels)
      CF$features$test$binning <- Bin(CF$test_paths, agg, channels)
    }
    else if ("test_data" %in% names(CF)){
      CF$features$test$binning <- Bin(CF$test_paths, agg, channels)
    }
  }
  
  if (featMethod == "EMD"){
    if ("ref_data" %in% names(CF)){
      CF$features$ref$EMD <- EMD(CF$ref_paths, agg, channels)
      CF$features$test$EMD <- EMD(CF$test_paths, agg, channels)
    }
    else if ("test_data" %in% names(CF)){
      CF$features$test$EMD <- EMD(CF$test_paths, agg, channels)
    }
  }
  
  if (featMethod == "fingerprint"){
    if ("ref_data" %in% names(CF)){
      CF$features$ref$fingerprint <- Fingerprint(CF$ref_paths, agg, channels, nRecursions)
      CF$features$test$fingerprint <- Fingerprint(CF$test_paths, agg, channels, nRecursions)
    }
    else if ("test_data" %in% names(CF)){
      CF$features$test$fingerprint <- Fingerprint(CF$test_paths, agg, channels, nRecursions)
    }
  }
  return(CF)
}


#' Calculate summary statistics (mean, SD, median, IQR)
#'
#' @param input List of FCS file paths
#' @param channels Channels to calculate features for
#'
#' @return Dataframe with features (columns) for samples (rows)
SummaryStats <- function(input, channels){
  all_stats <- list()
  for (path in input){
    ff <- ReadInput(path, n = NULL)
    stats <- list()
    for (channel in channels){
      stats[paste0(channel,"_mean")] <- base::mean(ff@exprs[,channel])
      stats[paste0(channel,"_sd")] <- stats::sd(ff@exprs[,channel])
      stats[paste0(channel,"_median")] <- stats::median(ff@exprs[,channel])
      stats[paste0(channel,"_IQR")] <- stats::IQR(ff@exprs[,channel])
    }
    all_stats[[path]] <- stats
  }
  stats <- data.frame(dplyr::bind_rows(all_stats), check.names = FALSE)
  return(stats)
}


#' Calculate earth mover's distance (EMD) between samples and aggregated dataset
#'
#' @param input List of FCS file paths
#' @param agg Dataframe of aggregated data
#' @param channels Channels to calculate features for
#'
#' @return Dataframe with features (columns) for samples (rows)
EMD <- function(input, agg, channels){
  all_stats <- list()
  for (path in input){
    ff <- ReadInput(path, n = NULL)
    stats <- list()
    for (channel in channels){
      stats[paste0(channel,'_', 'EMD')] <- transport::wasserstein1d(ff@exprs[, channel], 
                                                                    agg[, channel])
    }
    all_stats[[path]] <- stats
  }
  stats <- data.frame(dplyr::bind_rows(all_stats), check.names = FALSE)
  return(stats)
}


#' Calculate flowFP bin counts for samples based on aggregated data
#'
#' @param input List of FCS file paths
#' @param agg Dataframe of aggregated data
#' @param channels Channels to calculate features for
#' @param nRecursions Amount of flowFP recursions to use (default = 4)
#'
#' @return Dataframe with features (columns) for samples (rows)
Fingerprint <- function(input, agg, channels, nRecursions = 4){
  # Convert aggregated matrix to flowframe
  agg <- flowCore::flowFrame(agg[,channels])
  
  message("Fitting flowFP fingerprinting model")
  model <- flowFP::flowFPModel(agg, parameters = channels, 
                               nRecursions = nRecursions)
  all_stats <- list()
  for (path in input){
    ff <- ReadInput(path, n = NULL)
    call = flowFP::flowFP(ff[, channels], model)
    bin_counts = flowFP::counts(call)
    # Convert counts to bin frequencies
    stats = data.frame(t(apply(bin_counts, 1, function(x) x/sum(x))))
    all_stats[[path]] <- stats
  }
  stats <- data.frame(dplyr::bind_rows(all_stats), check.names = FALSE)
  return(stats)
}


#' Calculate flowFP bin counts for samples based on aggregated data
#'
#' @param input List of FCS file paths
#' @param agg Dataframe of aggregated data
#' @param channels Channels to calculate features for
#'
#' @return Dataframe with features (columns) for samples (rows)
Bin <- function(input, agg, channels){
  # Determine the bins on the aggregated data
  message("Determining bin boundaries on aggregated data")
  bin_boundaries <- apply(agg[,channels], 2, function(x) stats::quantile(x, probs = seq(0, 1, by = 0.1)))
  bin_boundaries <- data.frame(bin_boundaries)
  
  all_stats <- list()
  for (path in input){
    ff <- ReadInput(path, n = NULL)
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


#' Calculate peak features
#'
#' @param input List of FCS file paths
#' @param channels Channels to calculate features for
#'
#' @return Dataframe with features (columns) for samples (rows)
PeakExtraction <- function(input, channels){
  all_features <- list()
  for (path in input){
    ff <- ReadInput(path, n = NULL)
    features <- list()
    for (channel in channels){
      # Identify the 0.01 and 0.999 quantiles to cut the KDE density
      # This prevents the findpeaks function from identifying small peaks at extremes
      min_x <- quantile(ff@exprs[,channel], 0.001)
      max_x <- quantile(ff@exprs[,channel], 0.999)
      # Calculate y-axis density
      dens <- stats::density(ff@exprs[,channel], from = min_x, to = max_x)
      # Identify peaks
      peaks <- pracma::findpeaks(dens$y, minpeakheight = 0.01)
      # Account for cases where only one peak is detected
      if (nrow(peaks) == 1){
        min_x_peak <- peaks[1,]
        max_x_peak <- peaks[1,]
        min_y_peak <- peaks[1,]
        max_y_peak <- peaks[1,]
      }
      else{
        # Identify the smallest and largest peak based on the x-axis
        min_x_peak <- peaks[1,]
        max_x_peak <- peaks[nrow(peaks),]
        # Identify the smallest and largest peak based on the y-axis
        sorted_idx <- order(peaks[, 1])
        sorted_peaks <- peaks[sorted_idx, ]
        min_y_peak <- sorted_peaks[1,]
        max_y_peak <- sorted_peaks[nrow(peaks),]
      }
      # Get the expression values of all the peaks
      features[paste0(channel,"_min_x_peak_center")] <- dens$x[min_x_peak[2]]
      features[paste0(channel,"_min_x_peak_start")] <- dens$x[min_x_peak[3]]
      features[paste0(channel,"_min_x_peak_end")] <- dens$x[min_x_peak[4]]
      features[paste0(channel,"_max_x_peak_center")] <- dens$x[max_x_peak[2]]
      features[paste0(channel,"_max_x_peak_start")] <- dens$x[max_x_peak[3]]
      features[paste0(channel,"_max_x_peak_end")] <- dens$x[max_x_peak[4]]
      features[paste0(channel,"_min_y_peak_center")] <- dens$x[min_y_peak[2]]
      features[paste0(channel,"_min_y_peak_start")] <- dens$x[min_y_peak[3]]
      features[paste0(channel,"_min_y_peak_end")] <- dens$x[min_y_peak[4]]
      features[paste0(channel,"_max_y_peak_center")] <- dens$x[max_y_peak[2]]
      features[paste0(channel,"_max_y_peak_start")] <- dens$x[max_y_peak[3]]
      features[paste0(channel,"_max_y_peak_end")] <- dens$x[max_y_peak[4]]
    }
    all_features[[path]] <- features
  }
  features <- data.frame(dplyr::bind_rows(all_features), check.names = FALSE)
  return(features)
}