# Libraries needed

library(ggplot2)
library(dplyr)
library(TADkitdev2)
library(BiocGenerics)
library(GenomicRanges)


# Insulation score: Predicted

# Insulation score for predicted matrices : window is given in bin
insScorePredicted <- function(matrix, start, stop, bin.width, window) {
  
  w <- window  # To avoid heavy formula depending on the window
  first <- start / bin.width + w  # = 1st bin + window
  last <- stop / bin.width - w  # = last bin - window
  
  # Creation of a 4-columns data frame: chr, location of the beginning and the end of the bin and insulation score
  IS <- sapply (first:last, function(i) {
    mean(mat[(i-w):(i-1), (i+1):(i+w)])
  })
  
  startBin <- first:last * bin.width - bin.width  # Location of the beginning of the bin
  stopBin <- first:last * bin.width  # Location of the end of the bin
  IS.df <- data.frame(chr <- "1", startBin <- startBin, stopBin <- stopBin, IS <- IS)  # Creation of the data frame
  
  return(IS.df)
}


# Insulation score: Observed

#Insulation score for observed data matrices : window is given in bp and converted in bin within the function
insScoreObserved <- function(matrix, start, stop, bin.width, window) {
  
  w <- window / bin.width  # To avoid heavy formula depending on the window
  first <- start / bin.width + w  # = 1st bin + window
  last <- stop / bin.width - w  # = last bin - window
  
  # Creation of a 4-columns data frame: chr, location of the beginning and the end of the bin and insulation score
  IS <- sapply (first:last, function(i) {
    mean(mat[(i-w):(i-1), (i+1):(i+w)])
  })
  
  startBin <- first:last * bin.width - bin.width  # Location of the beginning of the bin
  stopBin <- first:last * bin.width  # Location of the end of the bin
  IS.df <- data.frame(chr <- "1", startBin <- startBin, stopBin <- stopBin, IS <- IS)  # Creation of the data frame
  
  return(IS.df)
}
