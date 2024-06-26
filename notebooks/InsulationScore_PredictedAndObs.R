# Library needed

library(ggplot2)
library(dplyr)
library(TADkitdev2)
library(BiocGenerics)
library(GenomicRanges)


# Insulation score for predicted matrices

# Insulation score for predicted matrices : window is given in bin
insScorePredicted <- function(matrix, start, stop, bin.width, window) {
  
  w <- window
  first <- start / bin.width + w
  last <- stop / bin.width - w
  
  IS <- sapply (first:last, function(i) {
    mean(mat[(i-w):(i-1), (i+1):(i+w)])
  })
  
  startBin <- first:last * bin.width - bin.width
  stopBin <- first:last * bin.width
  IS.df <- data.frame(chr <- "1", startBin <- startBin, stopBin <- stopBin, IS <- IS)
  
  return(IS.df)
}
