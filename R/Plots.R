PrepareFeaturePlot <- function(CF, featMethod, includeRef = "auto", 
                               flagMethod = NULL, flagSlot = "auto"){
  testFeatures <- CF$features$test[[featMethod]]
  refFeatures <- NULL
  testAnomaly <- NULL
  featNames <- colnames(testFeatures)
  # Check if reference is included in object
  if (includeRef == "auto") {
    includeRef <- "ref" %in% names(CF$features)
  }
  if (includeRef){
    refFeatures <- CF$features$ref[[featMethod]]
  }
  
  if (!is.null(flagMethod)){
    # Detect if novelties and/or outliers have already been flagged
    if (flagSlot == "auto"){
      if ("novelties" %in% names(CF)){
        flagSlot <- "novelties"
      }
      else if ("outliers" %in% names(CF)){
        flagSlot <- "outliers"
      }
    }
    testAnomaly <- CF[[flagSlot]][[flagMethod]]
  }
  return(list("refFeatures" = refFeatures,
              "testFeatures" = testFeatures,
              "testAnomaly" = testAnomaly))
}


#' Plot Heatmap
#'
#' @param CF CytoFlag object
#' @param featMethod Which generated features to plot
#' @param includeRef Whether to visualize reference data (default = "auto")
#' @param flagMethod Which flagging method to plot (default = NULL)
#' @param flagSlot Whether to use outlier or novelty flags (default = "auto")
#'
#' @return Heatmap plot
PlotHeatmap <- function(CF, featMethod, includeRef = "auto", flagMethod = NULL, 
                    flagSlot = "auto"){
  features <- PrepareFeaturePlot(CF, featMethod, includeRef, flagMethod, flagSlot)
  testFeatures <- features$testFeatures
  featNames <- colnames(testFeatures)
  refFeatures <- features$refFeatures
  testAnomaly <- features$testAnomaly
  
  if (!is.null(refFeatures)){
    refFeatures$category <- "reference"
    testFeatures$category <- "test"
    if (!is.null(testAnomaly)){
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
    if (!is.null(testAnomaly)){
      testFeatures$anomaly <- testAnomaly
    }
    else{
      testFeatures$anomaly <- "NA"
    }
    features <- testFeatures
  }
  
  # Annotation data
  annot_data <- features[,c("category", "anomaly")]
  colnames(annot_data) <- c("category", "anomaly")
  annot_colors <- list(category = c("reference" = "green",
                                    "test" = "blue"),
                       anomaly = c("NA" = "grey",
                                   "TRUE" = "red",
                                   "FALSE" = "blue"))

  # Prepare pheatmap
  mat <- as.matrix(features[,featNames])
  rownames(mat) <- seq(1, nrow(mat))

  # Plot
  p <- ComplexHeatmap::pheatmap(mat[,featNames], cluster_cols = FALSE,
                           color = colorRampPalette(rev(c("red", "white", "blue")))(100), 
                           show_rownames = TRUE, show_colnames = TRUE, 
                           scale = "column",
                           annotation_row = annot_data,
                           annotation_colors = annot_colors)
  return(p)
}
  

#' Plot PCA
#'
#' @param CF CytoFlag object
#' @param featMethod Which generated features to plot
#' @param includeRef Whether to visualize reference data (default = "auto")
#' @param flagMethod Which flagging method to plot (default = NULL)
#' @param flagSlot Whether to use outlier or novelty flags (default = "auto")
#' @param arrows Whether to plot loadings (default = FALSE)
#' @param labels Whether to plot labels of anomalies (default = TRUE)
#'
#' @return PCA plot
PlotPCA <- function(CF, featMethod, includeRef = "auto", flagMethod = NULL, 
                    flagSlot = "auto", arrows = FALSE, labels = TRUE){
  features <- PrepareFeaturePlot(CF, featMethod, includeRef, flagMethod, flagSlot)
  testFeatures <- features$testFeatures
  featNames <- colnames(testFeatures)
  refFeatures <- features$refFeatures
  testAnomaly <- features$testAnomaly

  if (!is.null(refFeatures)){
    # Reduce to first two principal components
    PCAOutput <- reduceDim(testFeatures = testFeatures[,featNames], 
                           flagStrat = "novelty",
                           refFeatures[,featNames])
    testFeatures <- PCAOutput[["testFeatures"]]
    refFeatures <- PCAOutput[["refFeatures"]]
    testFeatures$index <- rownames(testFeatures)
    refFeatures$index <- rownames(refFeatures)
    refFeatures$category <- "reference"
    if (!is.null(testAnomaly)){
      testFeatures$anomaly <- testAnomaly
      refFeatures$anomaly <- FALSE
    }
    testFeatures$category <- "test"
    features <- rbind(testFeatures, refFeatures)
  }
  else{
    # Reduce to first two principal components
    PCAOutput <- reduceDim(testFeatures = testFeatures[,featNames], 
                           flagStrat = "outlier")
    testFeatures <- PCAOutput[["testFeatures"]]
    testFeatures$index <- rownames(testFeatures)
    refFeatures <- PCAOutput[["refFeatures"]]
    testFeatures$category <- "test"
    if (!is.null(testAnomaly)){
      testFeatures$anomaly <- testAnomaly
    }
    features <- testFeatures
  }
  loadings <- PCAOutput[["loadings"]]
  loadings$variable <- rownames(loadings)
  PC1_var <- PCAOutput[["PC1_var"]]
  PC2_var <- PCAOutput[["PC2_var"]]

  # Plot
  loadScale <- 3
  p <- ggplot2::ggplot(features, ggplot2::aes(x = PC1, y = PC2, colour = category)) +
    ggplot2::xlab(paste0("PC1 (", round(PC1_var, 1), "%)")) +
    ggplot2::ylab(paste0("PC2 (", round(PC2_var, 1), "%)")) +
    ggplot2::theme(panel.background = ggplot2::element_blank(), 
          plot.background = ggplot2::element_blank(), 
          panel.border = ggplot2::element_rect(colour = "black", fill = NA, 
                                      linewidth = 0.5))
  if (is.null(testAnomaly)){
    p <- p + ggplot2::geom_point(size = 3)
  }
  else {
    p <- p + ggplot2::geom_point(size = 3, ggplot2::aes(shape = anomaly))
    
    if (labels){
      p <- p + ggplot2::geom_label(data = subset(features, anomaly == TRUE),
                                   ggplot2::aes(label = index), 
                                   nudge_x = 0.3, nudge_y = 0.3, 
                                   label.padding = unit(0.2, "lines"), 
                                   label.size = 0.25, # Add border of size 0.25
                                   label.r = unit(0.1, "lines"), # Round corners slightly
                                   fill = "white", 
                                   color = "black")
    }
  }
  
  if (arrows){
    p <- p + ggplot2::geom_segment(data = loadings, ggplot2::aes(x = 0, y = 0, xend = (PC1*loadScale), 
                                                                 yend = (PC2*loadScale)),
                                   arrow = ggplot2::arrow(length = ggplot2::unit(1/2, "picas")), 
                                   color = "black") +
      ggplot2::annotate("text", x = (loadings$PC1*loadScale), y = (loadings$PC2*loadScale),
                        label = loadings$variable)
  }
  return(p)
}


PlotAnomaly <- function(CF, idx, channel, n = 500, includeRef = "auto"){
  # Read the anomalous data
  ff_anom <- ReadInput(CF$test_paths[idx], n = n)
  df_anom <- data.frame(ff_anom@exprs[,channels], check.names=FALSE)
  df_anom$index <- as.character(idx)
  df_anom$category <- "Outlier"
  
  if (includeRef == "auto"){
    includeRef <- "ref" %in% names(CF$features)
  }
  # Sample n other samples from test data
  set.seed(42)
  if (includeRef){
    test_indices <- 1:length(CF$ref_paths)
    background_label <- "Reference"
  } else {
    test_indices <- setdiff(1:length(CF$test_paths), idx)
    background_label <- "Other test samples"
  }
  if (length(test_indices) > 20){
    test_indices <- sample(test_indices, 20, replace = FALSE)
  }
  
  test_data <- list()
  for (test_idx in test_indices){
    print(test_idx)
    ff <- ReadInput(CF$test_paths[[test_idx]], n  = n)
    df <- data.frame(ff@exprs[,channels], check.names=FALSE)
    df$index <- as.character(test_idx)
    df$category <- background_label
    test_data[[test_idx]] <- df
  }
  test_data <- do.call("rbind", test_data)
  all_data <- rbind(df_anom, test_data)
  
  ggplot2::ggplot(all_data, ggplot2::aes(x = .data[[channel]], y = .data[["index"]],
                                         fill = .data[["category"]])) +
    ggridges::geom_density_ridges(scale = 5) +
    ggplot2::scale_fill_manual(values = c("Outlier" = "red", "Reference" = "blue",
                                          "Other test samples" = "blue")) +
    ggplot2::labs(y = "Index") +
    ggridges::theme_ridges() # to make it pretty
}
