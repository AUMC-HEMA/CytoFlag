#' @export
generateChannel <- function(nCells){
  proportionPos <- runif(1, min = 0.2, max = 0.8)
  nPos <- round(proportionPos * nCells)
  gauss1 <- stats::rnorm(nPos, mean = -2, sd = 1)
  gauss2 <- stats::rnorm(nCells - nPos, mean = 2, sd = 1)
  channel <- sample(c(gauss1, gauss2))
  return(channel)
}


#' @export
generateFlowframe <- function(nCells, nChannels){
  # Generate expression data
  exprs <- matrix(0, nrow = nCells, ncol = nChannels)
  for (i in 1:nChannels) {
    exprs[, i] <- generateChannel(nCells)
  }
  dimnames(exprs)[[2]] <- LETTERS[1:nChannels]

  # Create metadata
  meta <- data.frame(name = dimnames(exprs)[[2]],
                     desc = dimnames(exprs)[[2]])
  meta$range <- apply(apply(exprs, 2 ,range), 2, diff)
  meta$minRange <- apply(exprs, 2, min)
  meta$maxRange <- apply(exprs, 2, max)
  
  # Create spillover matrix
  spillover <- diag(1, nrow = nChannels, ncol = nChannels)
  dimnames(spillover)[[2]] <- dimnames(exprs)[[2]]
  desc <- list("SPILL" = spillover)
  
  # Create flowframe
  ff <- new("flowFrame",
            exprs = exprs,
            parameters = Biobase::AnnotatedDataFrame(meta),
            description = desc)
  # Apply a reverse logicle transformation
  logicle  <- flowCore::logicleTransform()
  tfList <- flowCore::transformList(colnames(ff), 
                                    flowCore::inverseLogicleTransform(trans = logicle))
  ff <- flowCore::transform(ff, tfList)
  return(ff)
}


#' Generate a folder with a simulated cytometry cohort
#'
#' @param dir Path of directory where FCS files should be saved
#' @param nFiles Number of FCS files to generate
#' @param nCells Number of cells to simulate per FCS file
#' @param nChannels Number of channels to generate per FCS file
#' 
#' @export
generateDemo <- function(dir, nFiles, nCells, nChannels){
  # Create output directory
  if (!dir.exists(dir)){
    dir.create(dir, recursive = TRUE)
  }
  for (i in seq(1, nFiles)){
    ff <- generateFlowframe(nCells, nChannels)
    write.FCS(ff, paste0(dir, "/", as.character(i), ".fcs"))
  }
}