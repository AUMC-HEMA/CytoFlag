#' Initialize CytoFlag object
#' 
#' Function to initialize CytoFlag workflow. Returns CytoFlab object
#' 
#' @export
CytoFlag <- function(){
  CF <- list(list())
  class(CF) <- "CytoFlag"
  CF[["preprocess_function"]] <- ProcessInput
  CF[["parallel_vars"]] <- c("channels", "CF", "ReadInput")
  CF[["parallel_packages"]] <- c("flowCore", "PeacoQC")
  return(CF)
}


AddLabels <- function(CF, labels, slot){
  CF[["labels"]][[slot]] <- labels
  return(CF)
}
 

AddTestLabels <- function(CF, labels){
  CF <- AddLabels(CF, labels, "test")
  return(CF)
}


AddReferenceLabels <- function(CF, labels){
  CF <- AddLabels(labels, "ref")
  return(CF)
}


ProcessInput <- function(ff){
  ff <- PeacoQC::RemoveMargins(ff, channels)
  ff <- flowCore::compensate(ff, ff@description$SPILL)
  ff <- flowCore::transform(ff, flowCore::transformList(colnames(ff@description$SPILL), 
                                                        flowCore::arcsinhTransform(a = 0, b = 1/150, c = 0)))
  return(ff)
}


#' Read, preprocess and downsample FCS files
#'
#' @param CF CytoFlag object
#' @param path Location of FCS file
#' @param n Number of cells to read from FCS file
#'
#' @return flowFrame
ReadInput <- function(CF, path, n = NULL){
  set.seed(42)
  ff <- flowCore::read.FCS(path, which.lines = n)
  ff <- CF[["preprocess_function"]](ff)
  return(ff)
}


#' Add data to an initialized CytoFlag object
#'
#' @param CF CytoFlag object
#' @param input List of FCS file paths
#' @param slot Slot in the CytoFlag object to store data
#' @param aggSize Total number of cells to read (default = 10000)
#'
#' @return CytoFlag object
AddData <- function(CF, input, slot, aggSize = 10000){
  # Determine how many cells (n) to read per file
  n = round(aggSize / length(input))
  for (path in input){
    ff <- ReadInput(CF, path, n)
    CF[[slot]][[path]] <- data.frame(ff@exprs, check.names = FALSE)
  }
  return(CF)
}


#' Add reference data to initialized CytoFlag object
#'
#' @param CF CytoFlag object
#' @param input List of FCS file paths
#' @param read Whether to read FCS files already (default = FALSE)
#' @param aggSize Total number of cells to read
#'
#' @return CytoFlag object
#' 
#' @export
AddReferenceData <- function(CF, input, read = FALSE, aggSize = 10000){
  # You can force READ to already load the reference data
  # Otherwise this is loaded every time you call FeatureGeneration
  CF$ref_paths <- input
  CF$aggSize <- aggSize
  if (read == TRUE){
    CF <- AddData(CF, input, "ref_data", aggSize)
  }
  return(CF)
}


#' Add reference data to initialized CytoFlag object
#'
#' @param CF CytoFlag object
#' @param input List of FCS file paths
#' @param read Whether to read FCS files already (default = FALSE)
#' @param aggSize Total number of cells to read
#'
#' @return CytoFlag object
#' 
#' @export
AddTestData <- function(CF, input, read = FALSE, aggSize = 10000){
  # You can force READ to already load the reference data
  # Otherwise this is loaded every time you call FeatureGeneration
  CF$test_paths <- input
  CF$aggSize <- aggSize
  if (read == TRUE){
    CF <- AddData(CF, input, "test_data", aggSize)
  }
  return(CF)
}


#' Flag anomalies in CytoFlag object
#'
#' @param CF CytoFlag object
#' @param featMethod Feature generation method to use.
#' Either one of: "summary", "EMD", "binning" or "fingerprint".
#' @param flagStrat Which anomaly detection strategy to use ("outlier" or "novelty").
#' By default, selects intended use based on availability of reference material.
#' @param flagMethod Which anomaly detection algorithm to use ("KDE" or "forest")
#' @param PCA Whether to reduce dimensionality to first two PCs. (default = TRUE)
#'
#' @return CytoFlag object
#' 
#' @export
Flag <- function(CF, featMethod = NULL, flagStrat = "auto", flagMethod = "KDE",
                 PCA = TRUE){
  # Detect  flagging strategy if not supplied based on CytoFlag object
  if (flagStrat == "auto"){
    if ("test" %in% names(CF$features) & "ref" %in% names(CF$features)){
      message("Detected reference and test data, flagging samples using novelty detection.")
      flagStrat <- "novelty"
    }
    else if ("test" %in% names(CF$features)){
      message("Detected only test data, flagging samples using outlier detection.")
      flagStrat <- "outlier"
    }
    else {
      stop("Did not detect generated features in CytoFlag object.")
    }
  }
  
  # Warn of non-supported combinations of flagging and algorithms
  if (flagStrat == "novelty" & flagMethod == "forest"){
    stop("Isolation forest is not supported for novelty detection.")
  }
  
  # Get the features
  if (is.null(featMethod)){
    testFeatures <- CF$features$test[[names(CF$features$test)[1]]]
    if (flagStrat == "novelty"){
      refFeatures <- CF$features$ref[[names(CF$features$ref)[1]]]
    }
  }
  else{
    testFeatures <- CF$features$test[[featMethod]]
    if (flagStrat == "novelty"){
      refFeatures <- CF$features$ref[[featMethod]]
    }
  }
  
  if (PCA){
    message("Reducing features to first two principal components")
    if (flagStrat == "outlier"){
      PCAOutput <- reduceDim(testFeatures, "outlier")
      testFeatures <- PCAOutput[["testFeatures"]]
    }
    if (flagStrat == "novelty"){
      # Fit the PCA on the reference features before reducing the test features
      PCAOutput <- reduceDim(testFeatures, "novelty", refFeatures)
      testFeatures <- PCAOutput[["testFeatures"]]
      refFeatures <- PCAOutput[["refFeatures"]]
    }
  }
  
  if (flagMethod == "KDE"){
    if (flagStrat == "novelty"){
      message("Estimating reference density using KDE")
      # Based on: https://bookdown.org/egarpor/NP-UC3M/kde-ii-mult.html
      ref_scores <- ks::kde(x = refFeatures, eval.points = refFeatures)$estimate
      message("Calculating test density using reference KDE")
      test_scores <- ks::kde(x = refFeatures, eval.points = testFeatures)$estimate
      threshold <- stats::quantile(ref_scores, 0.05)
      message("Predicting novelties using probability cut-off")
      novelties <- as.factor(ifelse(test_scores < threshold, TRUE, FALSE))
      n <- length(novelties[novelties == TRUE])
      perc <- (n / length(novelties)) * 100
      fmt_str <- sprintf("Found %d novelties (%.1f%%)", n, perc)
      message(fmt_str)
      CF$novelties$KDE <- novelties
    }
    else if (flagStrat == "outlier"){
      message("Calculating test density using test KDE")
      message("Please note that this setting always returns X% most extreme values")
      test_scores <- ks::kde(x = testFeatures, eval.points = testFeatures)$estimate
      threshold <- stats::quantile(test_scores, 0.05)
      outliers <- as.factor(ifelse(test_scores >= 0.5, TRUE, FALSE))
      n <- length(outliers[outliers == TRUE])
      perc <- (n / length(outliers)) * 100
      fmt_str <- sprintf("Found %d outliers (%.1f%%)", n, perc)
      message(fmt_str)
      CF$outliers$KDE <- outliers
    }
  }
  
  if (flagStrat == "outlier"){
    if (flagMethod == "forest"){
      message("Building isolation forest")
      forest <- isotree::isolation.forest(testFeatures, seed = 42)
      scores <- stats::predict(forest, testFeatures)
      # Scores above 0.5 are most likely outliers according to the original paper
      # This also the default behavior in the sklearn implementation
      outliers <- as.factor(ifelse(scores >= 0.5, TRUE, FALSE))
      n <- length(outliers[outliers == TRUE])
      perc <- (n / length(outliers)) * 100
      fmt_str <- sprintf("Found %d outliers (%.1f%%)", n, perc)
      message(fmt_str)
      CF$outliers$forest <- outliers
    }
  }
  return(CF)
}