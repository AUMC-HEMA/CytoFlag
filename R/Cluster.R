ClusterFeatures <- function(CF, featMethod, k, includeRef = "auto"){
  # Get the features
  testFeatures <- CF$features$test[[featMethod]]
  featNames <- colnames(testFeatures)
  testFeatures$category <- "test"
  if (includeRef == "auto") {
    includeRef <- "ref" %in% names(CF$features)
  }
  if (includeRef){
    refFeatures <- CF$features$ref[[featMethod]]
    refFeatures$category <- "reference"
    features <- rbind(testFeatures, refFeatures)
  } else {
    features <- testFeatures
  }
  
  # Fit GMM for each k and evaluate BIC
  # Search for optimal K
  gmm <- mclust::Mclust(features[,featNames], G=k)
  labels <- gmm$classification

  features <- cbind(features, labels)

  # Add cluster labels to data
  CF$clusters$test[[featMethod]] <- features[features$category == "test", "labels"]
  if (includeRef){
    CF$clusters$ref[[featMethod]] <- features[features$category == "reference", "labels"]
  }
  return(CF)
}
