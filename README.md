# CytoFlag



### Loading data

Currently, the package takes a list of FCS files as input. Internally, data is **always** compensated and transformed (arcsinh: co-factor 150) before analysis.

```{r}
CF <- CytoFlag() 
CF <- AddTestData(CF, files) 
```

Loading data only adds the file paths to the object. You can also choose to read the data up front. 

```{r}
CF <- CytoFlag() 
CF <- AddTestData(CF, files, read = TRUE) 
```

If you want to perform novelty detection, you can also add reference data

```{r}
CF <- CytoFlag() 
CF <- AddReferenceData(CF, refFiles)
CF <- AddTestData(CF, testFiles)
```



### Feature generation

Feature generation is performed with the "FeatureGeneration" function.

| **featMethod** | **Description**                                              |
| -------------- | ------------------------------------------------------------ |
| "summary"      | Calculates summary statistics for all files (individually)   |
| "peaks"        | Calculates peak statistics for all files (individually)      |
| "EMD"          | Aggregates data and calculates earth mover's distance to aggregates data |
| "binning"      | Aggregates data, identifies bin boundaries and calculates bin percentages for all files |
| "fingerprint"  | Aggregates data, identifies high-dimensional bin boundaries (FlowFP) and calculates bin percentages for all files |

```R
CF <- CytoFlag() 
CF <- AddTestData(CF, files) 
CF <- FeatureGeneration(CF, channels = channels, featMethod = "peaks")
```

The features are added to the object separately for reference and test data.

Currently, if only test data is added, aggregated data is obtained from the test data. If reference data is available, then the reference data is used for aggregation.

```{r}
# Find generated features
features <- CF$test$peaks
```



### Anomaly detection

Anomaly detection is performed with the "Flag" function.

| **flagMethod** | Description                                                  |
| -------------- | ------------------------------------------------------------ |
| "forest"       | Identifies anomalies using isolation forest. Only works for outlier detection! |
| "KDE"          | Identifies anomalies using kernel density estimation. Only recommended for novelty detection (always filters 5% most extreme samples for outlier detection) |

You can also use PCA to reduce features to 2D. This is recommended for large feature sets in combination with KDE-based anomaly detection due to the curse of dimensionality. 

```{r}
# Detecting outliers based on peak patterns and isolation forest
CF <- CytoFlag() 
CF <- AddTestData(CF, files) 
CF <- FeatureGeneration(CF, channels = channels, featMethod = "peaks")
CF <- Flag(CF, flagMethod = "forest", PCA = FALSE)

# Detecting novelties based on earth mover's distance and KDE
CF <- CytoFlag() 
CF <- AddReferenceData(CF, refFiles)
CF <- AddTestData(CF, testFiles)
CF <- FeatureGeneration(CF, channels = channels, featMethod = "EMD")
CF <- Flag(CF, flagMethod = "KDE", PCA = TRUE)
```



### Visualizations

**WORK IN PROGRESS**
