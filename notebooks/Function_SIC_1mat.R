# Libraries

library(ggplot2)
library(dplyr)
library(TADkitdev2)
library(BiocGenerics)
library(GenomicRanges)


# Function: Structural impact score on one matrix

structImpactScore <- function(MUT.mat, WT.mat, start, stop, vpStart, vpStop, binWidth = 64e3) {
  # MUT.mat and WT.mat: obtained by loading matrix thanks to orca2matrix or readOrcaMatrix
  # start and stop: obtained thanks to orca2matrix
  # vpStart, vpStop: obtained thanks to the .bed file with the location of the mutation
  # binWIdth: obtained thanks to readOrcaMatrix
  
  
  # Viewpoint
  
  # Correcting the difference between what was generated with python in the .bed file and the function
  vp_start = vpStart + 1
  vp_stop = vpStop + 1
  
  # Create a list of matrices with the mutated one et the wildtype
  mat.lst <- list(MUT = MUT.mat, WT = WT.mat)
  
  # Set the viewpoint with the viewPointInteract function
  vp = viewPointInteract(
    matrix.lst = mat.lst,
    bin.width = binWidth,
    start = start, stop = stop,
    vp.start = vp_start, vp.stop = vp_stop,
    self_interaction = FALSE
  )
  
  
  # Structural impact score
  
  # Calculate the structural impact score and gives its mean within the viewpoint
  SIC = (vp$MUT$mean_interact/vp$WT$mean_interact) %>% abs %>% mean(., na.rm = TRUE)  #na.rm = remove NA
  
  # Convert the SIC into a data frame to plot a graph
  SIC.df = data.frame(chr = "1", start = vp_start, stop = vp_stop, SIC = SIC)
  
  
  return(SIC.df)
  
}






