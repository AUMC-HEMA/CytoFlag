#' Plot PCA
#'
#' @param CF CytoFlag object
#' @param featMethod Which generated features to plot
#' @param flagMethod Which flagging method to plot
#' @param includeRef Whether to visualize reference data
#' @param showFlags Whether to highlight anomalies (TRUE/FALSE)
#' @param arrows Whether to plot loadings (TRUE/FALSE)
#'
#' @return PCA plot
PCAPlot <- function(CF, featMethod, flagMethod, includeRef = "auto",
                    showFlags = "auto", arrows = FALSE){
  testFeatures <- CF$features$test[[featMethod]]
  featNames <- colnames(testFeatures)
  
  if (includeRef == "auto"){
    if ("ref" %in% names(CF$features)){
      includeRef <- TRUE
    }
    else {
      includeRef <- FALSE
    }
  }
  # Detect if novelties and/or outliers have already been flagged
  if (showFlags == "auto"){
    if ("novelties" %in% names(CF)){
      showFlags <- "novelty"
    }
    else if ("outliers" %in% names(CF)){
      showFlags <- "outlier"
    }
    else {
      showFlags <- FALSE
    }
  }

  if (showFlags == "novelty"){
    testAnomaly <- CF$novelties[[flagMethod]]
  }
  else if (showFlags == "outlier"){
    testAnomaly <- CF$outliers[[flagMethod]]
  }
  
  print(testFeatures)
  
  if (includeRef){
    refFeatures <- CF$features$ref[[featMethod]]
    # Reduce to first two principal components
    PCAOutput <- reduceDim(testFeatures = testFeatures[,featNames], 
                           flagStrat = "novelty",
                           refFeatures[,featNames])
    testFeatures <- PCAOutput[["testFeatures"]]
    refFeatures <- PCAOutput[["refFeatures"]]
    refFeatures$category <- "reference"
    if (showFlags != FALSE){
      testFeatures$anomaly <- testAnomaly
      refFeatures$anomaly <- FALSE
    }
    testFeatures$category <- "test"
    features <- rbind(testFeatures, refFeatures)
  }
  else{
    flagStrat <- "outlier"
    # Reduce to first two principal components
    PCAOutput <- reduceDim(testFeatures = testFeatures[,featNames], 
                           flagStrat = "outlier")
    testFeatures <- PCAOutput[["testFeatures"]]
    refFeatures <- PCAOutput[["refFeatures"]]
    testFeatures$category <- "test"
    if (showFlags != FALSE){
      testFeatures$anomaly <- testAnomaly
    }
    features <- testFeatures
  }
  loadings <- PCAOutput[["loadings"]]
  loadings$variable <- rownames(loadings)
  PC1_var <- PCAOutput[["PC1_var"]]
  PC2_var <- PCAOutput[["PC2_var"]]

  print(testFeatures)
  # Plot
  loadScale <- 3
  p <- ggplot2::ggplot(features, ggplot2::aes(x = PC1, y = PC2, colour = category)) +
    ggplot2::xlab(paste0("PC1 (", round(PC1_var, 1), "%)")) +
    ggplot2::ylab(paste0("PC2 (", round(PC2_var, 1), "%)")) +
    ggplot2::theme(panel.background = ggplot2::element_blank(), 
          plot.background = ggplot2::element_blank(), 
          panel.border = ggplot2::element_rect(colour = "black", fill = NA, 
                                      linewidth = 0.5))
  if (arrows){
    p <- p + ggplot2::geom_segment(data = loadings, ggplot2::aes(x = 0, y = 0, xend = (PC1*loadScale), 
                                               yend = (PC2*loadScale)),
                          arrow = ggplot2::arrow(length = ggplot2::unit(1/2, "picas")), 
                          color = "black") +
      ggplot2::annotate("text", x = (loadings$PC1*loadScale), y = (loadings$PC2*loadScale),
               label = loadings$variable)
  }
  if (showFlags == FALSE){
    p <- p + ggplot2::geom_point(size = 3)
  }
  else {
    p <- p + ggplot2::geom_point(size = 3, ggplot2::aes(shape = anomaly))
  }
  return(p)
}
