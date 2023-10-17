registerLandmarks <- function(CF, channel, aggSize = 10000, aggSlot = "auto",
                              adjust = 1, min_peak_height = 0.05, 
                              min_peak_dist = 5, plot = TRUE, save = TRUE){
  agg <- getAggregate(CF, aggSize = aggSize, aggSlot = aggSlot)
  data <- as.matrix(agg)
  
  # Kernel density estimation
  kde <- stats::density(data[,channel], adjust = adjust)
  bw <- kde$bw
  dens <- kde$y
  
  # Calculate the median percentile distance in the data
  # Depending on the distribution of the data, percentile distances can vary
  # We assume that the median is a robust estimate
  percentiles <- sapply(1:100, function(p) {
    stats::quantile(data[,channel], probs = p/100)})
  percentile_dist <- diff(percentiles)
  median_dist <- stats::median(percentile_dist)
  min_peak_dist <- median_dist * min_peak_dist
  # Calculate how many indices in the density profile this equates to
  min_peak_indices <- round(min_peak_dist / diff(kde$x)[1], 0)

  # Find peaks
  peaks <- pracma::findpeaks(dens, 
                             minpeakheight = min_peak_height,
                             minpeakdistance = min_peak_indices)
  peaks <- data.frame(peaks)
  colnames(peaks) <- c("peak_height", "index_max", "index_start", "index_end")
  
  # Add the expression values at the start, center, and end of the peaks
  peaks[,"exprs_max"] <- kde$x[peaks[,"index_max"]]
  peaks[,"exprs_start"] <- kde$x[peaks[,"index_start"]]
  peaks[,"exprs_end"] <- kde$x[peaks[,"index_end"]]
  
  if (save){
    CF$landmarks$ref[[channel]]$landmarks <- peaks
    CF$landmarks$ref[[channel]]$bw <- bw
  }
  
  if (plot){
    # Plot the reference distribution with peak information
    kde_ref <- data.frame(cbind(kde$x, kde$y))
    kde_ref$index <- 1
    colnames(kde_ref) <- c("x", "y", "index")
    # Mark all the peaks based on expression values
    mark_agg_peaks <- function(value) {
      for (i in 1:nrow(peaks)) {
        if (value >= peaks$exprs_start[i] & value <= peaks$exprs_end[i]) {
          return(i)
        }
      }
      return(NA)
    }
    kde_ref$peak <- as.factor(sapply(kde_ref$x, mark_agg_peaks))
    p1 <- ggplot2::ggplot(kde_ref, aes(x, index, height = y, group = index, fill = peak)) +
      ggridges::geom_ridgeline_gradient(scale=5) +
      ggplot2::theme_minimal() + 
      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank(),
            legend.position = "none") +
      ggplot2::labs(x = channel) +
      ggplot2::ggtitle("Aggregated data")
    
    # Run KDE for each individual file and mark peaks from aggregate
    kde_list <- list()
    cnt <- 0
    set.seed(42)
    for (i in sample(seq(1, length(CF$test_paths)), 20)){
      cnt <- cnt + 1
      # Read input
      ff <- ReadInput(CF, CF$test_paths[[i]], n=1000)
      df <- data.frame(ff@exprs, check.names=FALSE)
      df <- df[df[,channel] > stats::quantile(df[,channel], 0.001) & df[,channel] < stats::quantile(df[,channel], 0.999), ]
      
      # WARNING
      # Make sure to use the same bandwidth (???)
      kde_sample <- stats::density(df[,channel])
      kde_sample <- data.frame(cbind(kde_sample$x, kde_sample$y))
      kde_sample$index <- cnt
      colnames(kde_sample) <- c("x", "y", "index")
      
      # Mark all the peaks based on expression values
      mark_agg_peaks <- function(value) {
        for (i in 1:nrow(peaks)) {
          if (value >= peaks$exprs_start[i] & value <= peaks$exprs_end[i]) {
            return(i)
          }
        }
        return(NA)
      }
      kde_sample$peak <- sapply(kde_sample$x, mark_agg_peaks)
      kde_list[[cnt]] <- kde_sample
    }
    data <- do.call("rbind", kde_list)
    data$peak <- as.factor(data$peak)
    
    p2 <- ggplot2::ggplot(data, aes(x, index, height = y, group = index, fill = peak)) +
      ggridges::geom_ridgeline_gradient(scale=5) +
      ggplot2::theme_minimal() + 
      ggplot2::theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank()) +
      ggplot2::labs(x = channel) +
      ggplot2::ggtitle("Individual samples")
    
    gridExtra::grid.arrange(p1, p2, ncol=2)
  }
  return(CF)
}
