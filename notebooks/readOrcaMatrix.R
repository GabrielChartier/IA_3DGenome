# Libraries

library(ggplot2)
library(dplyr)
library(TADkitdev2)
library(BiocGenerics)
library(GenomicRanges)


# Function: read orca matrix

readOrcaMat <- function(path_predicted_mat, output = "OE") {
  
  # Extract the informations of the header
  
  # Reads the header and verify its format
  header <- readLines(path_predicted_mat, n = 1)
  verif <- substr(header, 1, 2)
  
  if (verif == "# ") {
    # Remove the "# " at the beginning of the header
    header <- sub(pattern = "# ", "", header)
    
    # Create a list with the informations in pairs
    pairs.lst <- unlist(strsplit(header, " "))
    
    # Create an empty list to store the previous pairs into a "key-value" format
    kv.lst <- list()
    
    # Build the list with each information (key-value)
    for (s in pairs.lst) {
      kv <- unlist(strsplit(s, "="))
      key <- kv[1]
      value <- kv[2]
      kv.lst[[key]] <- value
    }
    
    
    # Calculation of the missing data
    
    # Calculate the bin width thanks to the wpos form the header
    binWidth <- as.integer(kv.lst$wpos) / 250
    
    # Create a list with all calculated data (from the header)
    calcData <- list(bin_width = binWidth)
    
    # Create a list with all data : those from the header and the calculated
    allData = c(kv.lst, calcData)
    
    # Load the matrix 
    mat <- orca2matrix(
      df_prediction.path = path_predicted_mat,
      mpos = as.integer(kv.lst$mpos),
      scale = as.integer(kv.lst$wpos),
      chromsize = 158534110,  # Va être intégré au header ?
      output = "OE")
    
    # Create a list with: 1: all the data, 2: the matrix
    results <- list(parameters = allData, Matrix = mat)
    
    return(results)
    
    
    # If the format of the header is incorrect, return a message
  } else {
    return("Incorrect header format: it must be: # chr=x region=1:32000001 bin_width=64000 mpos=16000000 wpos=16000000")
  }
}
