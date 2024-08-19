#' @export
prepareFeatureplot <- function(CF, featMethod, plotRef, flagSlot){
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
#' @export
plotHeatmap <- function(CF, featMethod, plotRef = FALSE, 
                        flagSlot = NULL, plotLabels = FALSE){
  features <- prepareFeatureplot(CF, featMethod, plotRef, flagSlot)
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
#' @param featMethod Which features to use
#' @param fitData Which data to fit PCA on ("test", "reference" or "all")
#' @param plotData Which data to plot based on PCA fit ("test", "reference" or "all")
#' @param PCx Which PC to plot along x-axis
#' @param PCy Which PC to plot along y-axis
#' @param color Category used for point color ("outlier", "novelty", "labels")
#' @param shape Category used for point shape ("outlier", "novelty", "labels")
#' @param plotLoadings Whether to plot loadings along x and y-axis
#'
#' @return PCA plot
#' @export
plotPCA <- function(CF, featMethod, fitData, plotData, PCx = 1, PCy = 2,
                    color = NULL, shape = NULL, plotLoadings = TRUE){
  # Fit PCA on input data
  if (fitData == "test"){
    inputData <- CF$features$test[[featMethod]]
  } else if (fitData == "reference"){
    inputData <- CF$features$reference[[featMethod]]
  } else if (fitData == "all"){
    inputData <- rbind(CF$features$test[[featMethod]], CF$features$ref[[featMethod]])
  }
  pca <- stats::prcomp(inputData, center = TRUE, scale = TRUE)
  x_var <- pca$sdev[PCx]^2
  y_var <- pca$sdev[PCy]^2
  
  features <- colnames(inputData)
  
  if (plotData == "test" | plotData == "all"){
    testData <- CF$features$test[[featMethod]]
    testData$slot <- "test"
    if (is.null(color)){
      testData$color <- "test"
    } else if (color == "outlier"){
      testData$color <- CF$outliers[[featMethod]]
    } else if (color == "novelty"){
      testData$color <- CF$novelties[[featMethod]]
    } else if (color == "labels"){
      testData$color <- CF$labels$test
    }
    if (is.null(shape)){
      testData$shape <- "test"
    } else if (shape == "outlier"){
      testData$shape <- CF$outliers[[featMethod]]
    } else if (shape == "novelty"){
      testData$shape <- CF$novelties[[featMethod]]
    } else if (shape == "labels"){
      testData$shape <- CF$labels$test
    }
  }
  if (plotData == "reference" | plotData == "all"){
    referenceData <- CF$features$reference[[featMethod]]
    referenceData$slot <- "reference"
    if (is.null(color)){
      referenceData$color <- "reference"
    } else if (color == "outlier"){
      referenceData$color <- "reference"
    } else if (color == "novelty"){
      referenceData$color <- "reference"
    } else if (color == "labels"){
      referenceData$color <- CF$labels$reference
    }
    if (is.null(shape)){
      referenceData$shape <- "reference"
    } else if (shape == "outlier"){
      referenceData$shape <- "reference"
    } else if (shape == "novelty"){
      referenceData$shape <- "reference"
    } else if (shape == "labels"){
      referenceData$shape <- CF$labels$reference
    }
  }
  if (plotData == "test"){
    plotData <- testData
  } else if (plotData == "reference"){
    plotData <- referenceData
  } else if (plotData == "all"){
    plotData <- rbind(testData, referenceData)
  }
  output <- stats::predict(pca, newdata = plotData[, features])
  plotData$x <- output[, paste0("PC", as.character(PCx))]
  plotData$y <- output[, paste0("PC", as.character(PCy))]
  plotData <- data.frame(plotData, check.names = FALSE)
  
  PCAPlot <- ggplot2::ggplot(plotData, ggplot2::aes(x = x, y = y, color = color, shape = shape)) +
    ggplot2::geom_point() +
    ggplot2::xlab(paste0("PC", PCx, " (", round(x_var, 1), "%)")) +
    ggplot2::ylab(paste0("PC", PCy, " (", round(y_var, 1), "%)")) +
    ggplot2::theme(panel.background = ggplot2::element_blank(),
                   plot.background = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(colour = "black", fill = NA))
  
  if (plotLoadings){
    # Extract the loadings
    loadings <- as.data.frame(pca$rotation)
    sorted_PCx <- loadings[order(abs(loadings[, paste0("PC", as.character(PCx))]),
                                 decreasing = TRUE), ]
    sorted_PCy <- loadings[order(abs(loadings[, paste0("PC", as.character(PCy))]),
                                 decreasing = FALSE), ]
    loadings_x <- data.frame(Variable = rownames(sorted_PCx), 
                             Loading = sorted_PCx[, paste0("PC", as.character(PCx))])
    loadings_x$Variable <- factor(loadings_x$Variable, levels = loadings_x$Variable)
    loadings_y <- data.frame(Variable = rownames(sorted_PCy), 
                             Loading = sorted_PCy[, paste0("PC", as.character(PCy))])
    loadings_y$Variable <- factor(loadings_y$Variable, levels = loadings_y$Variable)
    
    xBar <- ggplot2::ggplot(loadings_x, ggplot2::aes(x = Variable, y = Loading)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
                     axis.title.x = ggplot2::element_blank(), 
                     axis.title.y = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(),
                     plot.background = ggplot2::element_blank(),
                     panel.border = ggplot2::element_rect(colour = "black", fill = NA))
    
    yBar <- ggplot2::ggplot(loadings_y, ggplot2::aes(x = Variable, y = Loading)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::coord_flip() +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(), 
                     axis.title.y = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(),
                     plot.background = ggplot2::element_blank(),
                     panel.border = ggplot2::element_rect(colour = "black", fill = NA))
    
    pGrid <- gridExtra::grid.arrange(yBar, PCAPlot, NULL, xBar,     
                                     ncol = 2,
                                     heights = c(1, 0.35),
                                     widths = c(0.35, 1)
    )
    return(pGrid)
  } else {
    return(PCAPlot)
  }
}


#' @export
plotAnomaly <- function(CF, file, channel, n = 500, includeRef = "auto"){
  # Read the anomalous data
  ff_anom <- readInput(CF, file, n = n)
  df_anom <- data.frame(ff_anom@exprs, check.names=FALSE)
  df_anom$index <- file
  df_anom$category <- "Outlier"
  
  if (includeRef == "auto"){
    includeRef <- "ref" %in% names(CF$features)
  }
  # Sample n other samples from test data
  set.seed(42)
  if (includeRef){
    test_files <- CF$paths$reference
    background_label <- "Reference"
  } else {
    test_files <- CF$paths$test[CF$paths$test != file]
    background_label <- "Other test samples"
  }
  if (length(test_files) > 20){
    test_files <- sample(test_files, 20, replace = FALSE)
  }
  
  test_data <- list()
  for (test_file in test_files){
    ff <- readInput(CF, test_file, n  = n)
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
