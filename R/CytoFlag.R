#' Initialize CytoFlag object
#' 
#' Function to initialize CytoFlag workflow. Returns CytoFlab object
#' 
#' @return Default initialized CytoFlag object
#' 
#' @export
CytoFlag <- function(){
  CF <- list()
  class(CF) <- "CytoFlag"
  CF[["preprocessFunction"]] <- processInput
  CF[["parallel"]] <- list("parallelVars" = c("channels", "CF", "readInput"),
                           "parallelPackages" = c("flowCore"))
  return(CF)
}

 
#' Add additional annotations for test data 
#' 
#' @param CF CytoFlag object
#' @param labels Vector of numeric or string labels
#'
#' @return CytoFlag object
#'
#' @seealso \code{\link{addReferencelabels}}
#'
#' @export
addTestlabels <- function(CF, labels){
  CF[["labels"]][["test"]] <- labels
  return(CF)
}


#' Add additional annotations for reference data 
#' 
#' @param CF CytoFlag object
#' @param labels Vector of numeric or string labels
#'
#' @return CytoFlag object
#'
#' @seealso \code{\link{addTestlabels}}
#'
#' @export
addReferencelabels <- function(CF, labels){
  CF[["labels"]][["ref"]] <- labels
  return(CF)
}


#' @noRd
processInput <- function(ff){
  spill <- ff@description$SPILL
  ff <- flowCore::compensate(ff, spill)
  ff <- flowCore::transform(ff, flowCore::transformList(colnames(spill), 
                                            flowCore::arcsinhTransform(a = 0, 
                                                                       b = 1/150, 
                                                                       c = 0)))
  return(ff)
}


#' @export
readInput <- function(CF, path, n = NULL){
  set.seed(42)
  ff <- flowCore::read.FCS(path, which.lines = n)
  ff <- CF[["preprocessFunction"]](ff)
  return(ff)
}


#' @noRd
addData <- function(CF, input, type, read, reload, aggSize){
  # Check if the paths are already stored in the CytoFlag object
  if (type %in% names(CF$paths) & !reload){
    for (path in input){
      if (!path %in% CF$paths[[type]]){
        # message(paste("Concatenating additional file path", path))
        CF$paths[[type]] <- c(CF$paths[[type]], path)
        if (read == TRUE){
          ff <- readInput(CF, path, CF$nAgg)
          CF$data[[type]][[path]] <- data.frame(ff@exprs, check.names = FALSE)
          }
        }
      }
    } else {
      CF$paths[[type]] <- input
      if (read == TRUE){
        CF$aggSize <- aggSize
        CF$nAgg <- ceiling(aggSize / length(input))
        for (path in input){
          ff <- readInput(CF, path, CF$nAgg)
          CF$data[[type]][[path]] <- data.frame(ff@exprs, check.names = FALSE)
        }
      }
    }
  return(CF)
}


#' Add test data to CytoFlag object
#'
#' @param CF CytoFlag object
#' @param input List of FCS file paths
#' @param read Whether load subsampled data already (default = FALSE)
#' @param reload Reload all data (default = FALSE)
#' @param aggSize Size of aggregate sample (default = 10000 cells)
#'
#' @return CytoFlag object
#'
#' @seealso \code{\link{addReferencedata}}
#' 
#' @export
addTestdata <- function(CF, input, read = FALSE, reload = FALSE,
                        aggSize = 10000){
  CF <- addData(CF, input, "test", read, reload, aggSize)
  return(CF)
}


#' Add reference data to CytoFlag object
#'
#' @param CF CytoFlag object
#' @param input List of FCS file paths
#' @param read Whether load subsampled data already (default = FALSE)
#' @param reload Reload all data (default = FALSE)
#' @param aggSize Size of aggregate sample (default = 10000 cells)
#'
#' @return CytoFlag object
#'
#' @seealso \code{\link{addTestdata}}
#' 
#' @export
addReferencedata <- function(CF, input, read = FALSE, reload = FALSE,
                             aggSize = 10000){
  CF <- addData(CF, input, "reference", read, reload, aggSize)
  return(CF)
}


#' Flag anomalies based on generated features 
#'
#' @param CF CytoFlag object
#' @param featMethod Which features to use for anomaly detection
#' @param flagStrat Which flaggins strategy to use ("outlier" or "novelty")
#'
#' @return CytoFlag object
#' 
#' @seealso \code{\link{generateFeatures}}
#'
#' @export
Flag <- function(CF, featMethod, flagStrat){
  testFeatures <- CF$features$test[[featMethod]]
  if (flagStrat == "outlier"){
    forest <- isotree::isolation.forest(testFeatures, sample_size = 1,
                                        ntrees = 1000,
                                        ndim = 1, seed = 42)
    scores <- stats::predict(forest, testFeatures)
    # Scores > 0.5 are most likely outliers according to the original paper
    outliers <- as.factor(ifelse(scores >= 0.5, TRUE, FALSE))
    CF$outliers[[featMethod]] <- outliers
  }
  if (flagStrat == "novelty"){
    refFeatures <- CF$features$ref[[featMethod]]
    # Based on: https://bookdown.org/egarpor/NP-UC3M/kde-ii-mult.html
    refScores <- ks::kde(x = refFeatures, eval.points = refFeatures)$estimate
    testScores <- ks::kde(x = refFeatures, eval.points = testFeatures)$estimate
    threshold <- stats::quantile(refScores, 0.05)
    novelties <- as.factor(ifelse(testScores < threshold, TRUE, FALSE))
    CF$novelties[[featMethod]] <- novelties
  }
  return(CF)
}
