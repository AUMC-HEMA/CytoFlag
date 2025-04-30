#' Plot heatmap
#'
#' @param CF CytoFlag object
#' @param featMethod Which features to use
#' @param plotData Which data to plot ("test", "reference" or "all")
#'
#' @return Heatmap
#' @export
plotHeatmap <- function(CF, featMethod, plotData){
  if (plotData == "test"){
    inputData <- CF$features$test[[featMethod]]
  } else if (fiplotDatatData == "reference"){
    inputData <- CF$features$reference[[featMethod]]
  } else if (plotData == "all"){
    inputData <- rbind(CF$features$test[[featMethod]], CF$features$ref[[featMethod]])
  }
  
  # Construct dendogram for raw generated features
  dend <- as.dendrogram(hclust(dist(inputData), method = "ward.D"))
  
  # Standardize the features within each channel
  channels <- CF$metadata[[featMethod]]$channels
  for (channel in channels){
    cols <- grep(channel, colnames(inputData), value = TRUE)
    mat <- as.matrix(inputData[,cols])
    inputData[,cols] <- (mat - mean(mat)) / sd(mat)
  }
  
  # Get # features per channnel
  split <- rep(1:length(channels), 
               each = ncol(inputData) / length(channels))
  print(split)
  
  # Create heatmap
  g <- ComplexHeatmap::Heatmap(
    as.matrix(inputData),
    cluster_rows = dend,    
    cluster_columns = FALSE, 
    column_split = split,    
    show_column_names = FALSE, 
    show_row_names = FALSE,
    column_title = channels,
    heatmap_legend_param = list(title = "Z-score"),
    row_dend_width = grid::unit(3, "cm"))
  return(g)
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
#' @param plotBars Whether to plot loading coefficients as bars along
#' @param plotArrows Whether to plot loadings as arrows
#' @param nLoadings Number of loadings to plot. Selects most influential
#'
#' @return PCA plot
#' @export
plotPCA <- function(CF, featMethod, fitData, plotData, PCx = 1, PCy = 2,
                    color = NULL, shape = NULL, plotBars = TRUE, 
                    plotArrows = TRUE, nLoadings = NULL){
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
  
  if (plotArrows) {
    loadings <- as.data.frame(pca$rotation[, c(PCx, PCy)])
    colnames(loadings) <- c("PCx", "PCy")
    loadings$feature <- rownames(loadings)
    if (!is.null(nLoadings)) {
      loadings$magnitude <- sqrt(loadings$PCx^2 + loadings$PCy^2)
      loadings <- loadings[order(-loadings$magnitude), ][1:nLoadings, ]
    }
    # Normalize and scale loadings for better visualization
    scale_factor <- min(max(abs(plotData$x)), max(abs(plotData$y))) * 1.5
    loadings$PCx <- loadings$PCx * scale_factor
    loadings$PCy <- loadings$PCy * scale_factor
    
    # Update geom_label_repel to position labels closer to the end of the arrows
    PCAPlot <- PCAPlot +
      ggplot2::geom_segment(data = loadings, 
                            ggplot2::aes(x = 0, y = 0, 
                                         xend = PCx * max(plotData$x) * 0.2, 
                                         yend = PCy * max(plotData$y) * 0.2),
                            arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm")), 
                            color = "black", size = 0.5) +
      ggrepel::geom_label_repel(data = loadings, 
                                ggplot2::aes(x = PCx * max(plotData$x) * 0.2 * 0.9, 
                                             y = PCy * max(plotData$y) * 0.2 * 0.9, 
                                             label = feature),
                                fill = scales::alpha("white", 1), 
                                color = "black", 
                                box.padding = 0.1,
                                point.padding = 0.2)
  }
  if (plotBars) {
    loadings <- as.data.frame(pca$rotation)
    loadings$magnitude <- sqrt(loadings[, paste0("PC", as.character(PCx))]^2 + 
                                 loadings[, paste0("PC", as.character(PCy))]^2)
    if (!is.null(nLoadings)) {
      loadings <- loadings[order(-loadings$magnitude), ][1:nLoadings, ]
    } else {
      loadings <- loadings[order(-loadings[, paste0("PC", as.character(PCx))]), ]
    }
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
    pGrid <- gridExtra::grid.arrange(yBar, PCAPlot, 
                                   ggplot2::ggplot() + ggplot2::theme_void(), 
                                   xBar,
                                   ncol = 2,
                                   heights = c(1, 0.35),
                                   widths = c(0.35, 1))
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
