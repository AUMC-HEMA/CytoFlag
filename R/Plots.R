PrepareFeaturePlot <- function(CF, featMethod, plotRef, flagSlot){
  testFeatures <- CF$features$test[[featMethod]]
  refFeatures <- NULL
  testAnomaly <- NULL
  testLabels <- NULL
  refLabels <- NULL

  # Check for the presence of labels
  if ("test" %in% names(CF$labels)){
    testLabels <- CF$labels$test
  }
  
  if (plotRef){
    refFeatures <- CF$features$ref[[featMethod]]
    if ("ref" %in% names(CF$labels)){
      refLabels <- CF$labels$ref
    }
  }
  
  if (!is.null(flagSlot)){
    if (featMethod %in% names(CF[[flagSlot]])){
      testAnomaly <- CF[[flagSlot]][[featMethod]]
    }
  }
  return(list("refFeatures" = refFeatures,
              "testFeatures" = testFeatures,
              "testAnomaly" = testAnomaly,
              "testLabels" = testLabels,
              "refLabels" = refLabels))
}


#' Plot Heatmap
#'
#' @param CF CytoFlag object
#' @param featMethod Which generated features to plot
#' @param plotRef Whether to visualize reference, if available (default = TRUE)
#' @param flagSlot Flags to plot (NULL, "outliers", "novelties") (default = NULL)
#' @param plotLabels Whether to plot supplied labels
#'
#' @return Heatmap plot
PlotHeatmap <- function(CF, featMethod, plotRef = FALSE, 
                        flagSlot = NULL, plotLabels = FALSE){
  features <- PrepareFeaturePlot(CF, featMethod, plotRef, flagSlot)
  testFeatures <- features$testFeatures
  featNames <- colnames(testFeatures)
  refFeatures <- features$refFeatures
  testAnomaly <- features$testAnomaly
  testLabels <- features$testLabels
  refLabels <- features$refLabels
  
  if (plotRef){
    refFeatures$category <- "reference"
    testFeatures$category <- "test"
    testFeatures$labels <- testLabels
    refFeatures$labels <- refLabels
    
    if (!is.null(flagSlot)){
      testFeatures$anomaly <- testAnomaly
      refFeatures$anomaly <- FALSE
    }
    else {
      testFeatures$anomaly <- "NA"
      refFeatures$anomaly <- "NA"
    }
    features <- rbind(testFeatures, refFeatures)
  }
  else{
    testFeatures$category <- "test"
    testFeatures$labels <- testLabels
    if (!is.null(testAnomaly)){
      testFeatures$anomaly <- testAnomaly
    }
    else{
      testFeatures$anomaly <- "NA"
    }
    features <- testFeatures
  }
  
  # Annotation data
  annot_colors <- list(category = c("reference" = "green",
                                    "test" = "blue"),
                       anomaly = c("NA" = "grey",
                                   "TRUE" = "red",
                                   "FALSE" = "blue"))
  cols <- c("category", "anomaly")
  if (plotLabels){
    cols <- c(cols, "labels")
  }
  annot_data <- features[, cols]
  colnames(annot_data) <- cols
  
  # Prepare pheatmap
  mat <- as.matrix(features[,featNames])
  rownames(mat) <- seq(1, nrow(mat))
  
  # Plot
  p <- ComplexHeatmap::pheatmap(mat[,featNames], cluster_cols = FALSE,
                                cluster_rows = TRUE,
                                # color = grDevices::colorRampPalette(rev(c("red", "white", "blue")))(100), 
                                show_rownames = TRUE, show_colnames = TRUE, 
                                annotation_row = annot_data,
                                annotation_colors = annot_colors)
  return(p)
}


#' Plot PCA
#'
#' @param CF CytoFlag object
#' @param featMethod Which generated features to plot
#' @param plotRef Whether to visualize reference, if available (default = TRUE)
#' @param flagSlot Flags to plot (NULL, "outliers", "novelties") (default = NULL)
#' @param ids Whether to plot indices of anomalies (default = FALSE)
#' @param plotLabels Whether to plot supplied labels, overrides flagSlot
#'
#' @return PCA plot
PlotPCA <- function(CF, featMethod, plotRef = FALSE,
                    flagSlot = NULL, ids = FALSE, plotLabels = FALSE){
  features <- PrepareFeaturePlot(CF, featMethod, plotRef, flagSlot)
  testFeatures <- features$testFeatures
  featNames <- colnames(testFeatures)
  refFeatures <- features$refFeatures
  testAnomaly <- features$testAnomaly
  testLabels <- features$testLabels
  refLabels <- features$refLabels
  
  if (plotRef){
    # Reduce to first two principal components
    PCAOutput <- reduceDim(testFeatures = testFeatures[, featNames],
                           flagStrat = "novelty",
                           refFeatures[, featNames])
    testFeatures <- PCAOutput[["testFeatures"]]
    refFeatures <- PCAOutput[["refFeatures"]]
    testFeatures$index <- seq(1, nrow(testFeatures))
    refFeatures$index <- seq(1, nrow(refFeatures))
    refFeatures$category <- "reference"
    refFeatures$labels <- refLabels

    if (!is.null(flagSlot)){
      testFeatures$anomaly <- testAnomaly
      refFeatures$anomaly <- FALSE
    }
    testFeatures$category <- "test"
    testFeatures$labels <- testLabels
    features <- rbind(testFeatures, refFeatures)
  }
  else{
    # Reduce to first two principal components
    PCAOutput <- reduceDim(testFeatures = testFeatures[,featNames],
                           flagStrat = "outlier")
    testFeatures <- PCAOutput[["testFeatures"]]
    testFeatures$index <- seq(1, nrow(testFeatures))
    testFeatures$category <- "test"
    testFeatures$labels <- testLabels
    if (!is.null(flagSlot)){
      testFeatures$anomaly <- testAnomaly
    }
    features <- testFeatures
  }
  loadings <- PCAOutput[["loadings"]]
  loadings$variable <- rownames(loadings)
  PC1_var <- PCAOutput[["PC1_var"]]
  PC2_var <- PCAOutput[["PC2_var"]]
  
  loadScale <- 3
  if (is.null(flagSlot)){
    # Don't color the dots if there are no labels and outliers to flag
    p <- ggplot2::ggplot(features, ggplot2::aes(x = .data[["PC1"]],
                                                y = .data[["PC2"]]))
  }
  # Color the dots by either labels or flagged anomalies
  if (plotLabels){
    color <- "labels"
  }
  else {
    color <- "anomaly"
  }
  p <- ggplot2::ggplot(features, ggplot2::aes(x = .data[["PC1"]],
                                              y = .data[["PC2"]],
                                              colour = .data[[color]])) 

  p <- p + ggplot2::xlab(paste0("PC1 (", round(PC1_var, 1), "%)")) +
           ggplot2::ylab(paste0("PC2 (", round(PC2_var, 1), "%)")) +
           ggplot2::theme(panel.background = ggplot2::element_blank(),
                          plot.background = ggplot2::element_blank(),
                          panel.border = ggplot2::element_rect(colour = "black", fill = NA,
                                                              linewidth = 0.5))
  if (plotRef) {
    p <- p + ggplot2::geom_point(size = 3, ggplot2::aes(shape = category))
  }
  else {
    p <- p + ggplot2::geom_point(size = 3)
  }
  
  if (ids){
    p <- p + ggplot2::geom_label(data = features,
                                 ggplot2::aes(label = index),
                                 nudge_x = 0.1, nudge_y = 0.1,
                                 label.padding = ggplot2::unit(0.2, "lines"),
                                 label.size = 0.25, # Add border of size 0.25
                                 label.r = ggplot2::unit(0.1, "lines"), # Round corners slightly
                                 fill = "white",
                                 color = "black")
  }
  return(p)
}


PlotAnomaly <- function(CF, file, channel, n = 500, includeRef = "auto"){
  # Read the anomalous data
  ff_anom <- ReadInput(CF, file, n = n)
  df_anom <- data.frame(ff_anom@exprs, check.names=FALSE)
  df_anom$index <- file
  df_anom$category <- "Outlier"
  
  if (includeRef == "auto"){
    includeRef <- "ref" %in% names(CF$features)
  }
  # Sample n other samples from test data
  set.seed(42)
  if (includeRef){
    test_files <- CF$ref_paths
    background_label <- "Reference"
  } else {
    test_files <- CF$test_paths[CF$test_paths != file]
    background_label <- "Other test samples"
  }
  if (length(test_files) > 20){
    test_files <- sample(test_files, 20, replace = FALSE)
  }
  
  test_data <- list()
  for (test_file in test_files){
    ff <- ReadInput(CF, test_file, n  = n)
    df <- data.frame(ff@exprs, check.names=FALSE)
    df$index <- as.character(test_file)
    df$category <- background_label
    test_data[[test_file]] <- df
  }
  test_data <- data.frame(dplyr::bind_rows(test_data), check.names = FALSE)
  all_data <- rbind(df_anom, test_data)
  
  ggplot2::ggplot(all_data, ggplot2::aes(x = .data[[channel]], y = .data[["index"]],
                                         fill = .data[["category"]])) +
    ggridges::geom_density_ridges(scale = 5, alpha = 0.7) +
    ggplot2::scale_fill_manual(values = c("Outlier" = "red", "Reference" = "blue",
                                          "Other test samples" = "blue")) +
    ggplot2::labs(y = "Index") +
    ggridges::theme_ridges() # to make it pretty
}
