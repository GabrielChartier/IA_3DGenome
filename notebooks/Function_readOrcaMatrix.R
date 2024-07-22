# Libraries

library(ggplot2)
library(dplyr)
library(TADkitdev2)
library(BiocGenerics)
library(GenomicRanges)


# Function: read orca matrix

readOrcaMat <- function(path_predicted_mat, path_norm_mat = NULL, output_type = "OE") {
  
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
    
    # Read data we want to return
    data_returned = list(resol = kv.lst$resol, mpos = as.integer(kv.lst$mpos), wpos = as.integer(kv.lst$wpos), start = as.integer(kv.lst$start), end = as.integer(kv.lst$end),
                         nbins = as.integer(kv.lst$nbins), width = as.integer(kv.lst$width), chromlen = as.integer(kv.lst$chromlen))
    
    
    # Calculation of the missing data
    
    # Calculate the bin width thanks to the width from the header
    binWidth <- data_returned$width / data_returned$nbins
    
    # Extract the data about the start and the end of the mutation in bp
    
    # Remove the "1:" before the location
    kv.lst$mutation <- sub(pattern = "1:", "", kv.lst$mutation)
    
    # Remove the "-" between the start and the end of the mutation and create new "keys and values"
    kv.lst$mutation <- unlist(strsplit(kv.lst$mutation, "-"))
    mutation_start <- as.integer(kv.lst$mutation[1]) + 1
    mutation_end <- as.integer(kv.lst$mutation[2]) + 1
    
    # Create a list with all calculated data (from the header)
    calcData <- list(bin_width = binWidth, mutation_start = mutation_start, mutation_end = mutation_end)
    
    # Create a list with all data : those from the header and the calculated
    allData = c(data_returned, calcData)
    
    # Adapt the output depending on the matrix wanted: Obs / Exp or Obs
    if (is.null(path_norm_mat)) {
      # Load the Obs / Exp matrix
      mat <- orca2matrix(
        df_prediction.path = path_predicted_mat,
        mpos = data_returned$mpos,
        scale = data_returned$width,
        chromsize = data_returned$chromlen,
        output = output_type)
    } else {
      # Load the Obs matrix 
      mat <- orca2matrix(
        df_prediction.path = path_predicted_mat,
        df_normmats.path = path_norm_mat,
        mpos = data_returned$mpos,
        scale = data_returned$width,
        chromsize = data_returned$chromlen,
        output = output_type)
    }
    
    # Create a list with: 1: all the data, 2: the matrix
    results <- list(parameters = allData, Matrix = mat)
    
    return(results)
    
    
    # If the format of the header is incorrect, return a message
  } else {
    return("Incorrect header format: it must be: # Orca=predictions resol=16Mb start=7936000 end=23936000 nbins=250 width=16000000 chromlen=158534110")
  }
}