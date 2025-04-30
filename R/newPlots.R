aggregatePlotdata <- function(CF, plotData, channel, n=1000, color=NULL,
                              featMethod = NULL){
  if (plotData == "all"){
    files <- c(CF$paths$test, CF$paths$reference)
  } else {
    files <- CF$paths[[plotData]]
  }
  
  data <- list()
  for (file in files){
    ff <- readInput(CF, file, n=n)
    df <- data.frame(ff@exprs[,channel], check.names=FALSE)
    colnames(df) <- "value"
    df$index <- file
    data[[file]] <- df
  }
  data <- data.frame(dplyr::bind_rows(data), check.names = FALSE)
  
  if (is.null(color)){
    data$group <- plotData
  } else if (color == "input"){
    data$group <- ifelse(data$index %in% CF$paths$test, "test",
                         ifelse(data$index %in% CF$paths$reference, "reference", NA))
  } else if (color == "outlier"){
    data$group <- ifelse(data$index %in% CF$paths$reference, "reference",
                         ifelse(data$index %in% names(CF$outliers[[featMethod]][CF$outliers[[featMethod]] == TRUE]), 
                                "outlier", "inlier"))
  } else if (color == "novelty"){
    data$group <- ifelse(data$index %in% CF$paths$reference, "reference",
                         ifelse(data$index %in% names(CF$novelties[[featMethod]][CF$novelties[[featMethod]] == TRUE]), 
                                "novelty", "non-novelty"))
  }
  return(data)
}


plotBoxplot <- function(CF, plotData, channel, n=1000, color=NULL,
                        featMethod = NULL){
  data <- aggregatePlotdata(CF, plotData="test", channel="A", color="novelty", featMethod = "quantiles")
  
  library(dplyr)
  # Sort the data from low to high median
  data <- data %>%
    group_by(index) %>%
    mutate(median_value = median(value)) %>%
    ungroup() %>%
    mutate(index = reorder(index, median_value))
  
  ggplot2::ggplot(data, ggplot2::aes(x = index, y = value, fill=group)) +
    ggplot2::geom_boxplot(color = "black", outlier.shape = NA) +
    ggplot2::theme_classic() +
    ggplot2::labs(title = paste("Expression of", channel), x = "File", y = "") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))
}