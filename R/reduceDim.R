#' Reduce features to first two principal components
#'
#' @param testFeatures Dataframe of features from test samples
#' @param flagStrat Which anomaly detection strategy to use ("outlier" or "novelty").
#' @param refFeatures Dataframe of features from reference samples
#'
#' @return List containing PC loadings, variance and reduced features
reduceDim <- function(testFeatures, flagStrat, refFeatures = NULL){
  if (flagStrat == "outlier"){
    pca <- stats::prcomp(testFeatures, center = TRUE, scale = TRUE)
    testFeatures$PC1 <- pca$x[, 1]
    testFeatures$PC2 <- pca$x[, 2]
    testFeatures <- testFeatures[c("PC1", "PC2")]
    return(list("loadings" = as.data.frame(pca$rotation[, 1:2]),
                "PC1_var" = pca$sdev[1]^2,
                "PC2_var" = pca$sdev[2]^2,
                "testFeatures" = testFeatures))
  }
  if (flagStrat == "novelty"){
    # Fit the PCA on the reference features before reducing the test features
    pca <- stats::prcomp(refFeatures, center = TRUE, scale = TRUE)
    test_PCA <- stats::predict(pca, newdata = testFeatures)
    testFeatures$PC1 <- test_PCA[,"PC1"]
    testFeatures$PC2 <- test_PCA[,"PC2"]
    testFeatures <- testFeatures[c("PC1", "PC2")]
    ref_PCA <- stats::predict(pca, newdata = refFeatures)
    refFeatures$PC1 <- ref_PCA[,"PC1"]
    refFeatures$PC2 <- ref_PCA[,"PC2"]
    refFeatures <- refFeatures[c("PC1", "PC2")]
    return(list("loadings" = as.data.frame(pca$rotation[, 1:2]),
                "PC1_var" = pca$sdev[1]^2,
                "PC2_var" = pca$sdev[2]^2,
                "testFeatures" = testFeatures,
                "refFeatures" = refFeatures))
  }
}
