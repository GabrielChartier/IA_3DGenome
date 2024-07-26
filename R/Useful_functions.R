# Libraries

library(ggplot2)
library(dplyr)
library(TADkitdev2)
library(BiocGenerics)
library(GenomicRanges)


#####################################################################################################################################################
                                                               # Read Orca Matrix #
#####################################################################################################################################################

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
    data_returned = list(resol = kv.lst$resol, mpos = (as.integer(kv.lst$mpos) - as.integer(kv.lst$start) + 1), wpos = as.integer(kv.lst$wpos), start = (as.integer(kv.lst$start) - as.integer(kv.lst$start) + 1), 
                         end = (as.integer(kv.lst$end) - as.integer(kv.lst$start) + 1), nbins = as.integer(kv.lst$nbins), width = as.integer(kv.lst$width), chromlen = as.integer(kv.lst$chromlen))
    
    
    # Calculation of the missing data
    
    # Calculate the bin width thanks to the width from the header
    binWidth <- data_returned$width / data_returned$nbins
    
    # Extract the data about the start and the end of the mutation in bp
    
    # Remove the "1:" before the location
    kv.lst$mutation <- sub(pattern = "1:", "", kv.lst$mutation)
    
    # Remove the "-" between the start and the end of the mutation and create new "keys and values"
    kv.lst$mutation <- unlist(strsplit(kv.lst$mutation, "-"))
    mutation_start <- (as.integer(kv.lst$mutation[1]) - as.integer(kv.lst$start))
    mutation_end <- (as.integer(kv.lst$mutation[2]) - as.integer(kv.lst$start))
    
    # Create a list with all calculated data (from the header)
    calcData <- list(bin_width = binWidth, mutation_start = mutation_start, mutation_end = mutation_end)
    
    # Create a list with all data : those from the header and the calculated
    allData = c(data_returned, calcData)
    
    # Adapt the output depending on the matrix wanted: Obs / Exp or Obs
    # Use the function orca2matrix from the TADkitdev2 package
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
    
    mat <- mat[(data_returned$start / binWidth + 1):(data_returned$end / binWidth), (data_returned$start / binWidth + 1):(data_returned$end / binWidth)]
    
    # Create a list with: 1: all the data, 2: the matrix
    results <- list(parameters = allData, Matrix = mat)
    
    return(results)
    
    
    # If the format of the header is incorrect, return a message
  } else {
    return("Incorrect header format: it must be: # Orca=predictions resol=16Mb start=7936000 end=23936000 nbins=250 width=16000000 chromlen=158534110")
  }
}


#####################################################################################################################################################
                                                            # Local Insulation Score #
#####################################################################################################################################################

#viewpoint in bp; window: default value <- 256e3: can be modified; bin.width (<-resolution): default value <- 64e3: can be modified
localInsScore <- function(matrix, view.point, window = 256e3, bin.width = 64e3) {
  
  #conversion of the viewpoint from bp to bin
  vpBin <- round(view.point / bin.width)
  #position of the beginning of the viewpoint (in bin)
  start <- vpBin - round(window / bin.width) + 1
  #position of the end of the viewpoint (in bin)
  stop <- vpBin + round(window / bin.width)
  
  if (start < 1 | stop > 250) {
    
    return(NA)
    
  } else {
    
    #viewpoint converted as a "sub"-matrix
    subMat <- matrix[start:stop, start:stop]
    
    #make the matrix symmetrical
    subMatSym <- as.matrix(subMat)
    subMatSym[lower.tri(subMatSym)] <- t(subMatSym)[lower.tri(subMatSym)]
    
    #definition of the position of A, B & C with (4 zones in the "sub"-matrix):
    #   A   C
    #   C   B
    #then calculate the mean number of interactions within each zone
    
    subStart <- 1
    subStop <- nrow(subMatSym)
    mid1 <- nrow(subMatSym) / 2
    mid2 <- mid1 + 1
    
    A <- subMatSym[subStart:mid1, subStart:mid1] %>% mean(.)
    
    B <- subMatSym[mid2:subStop, mid2:subStop] %>% mean(.)
    
    C <- subMatSym[subStart:mid1, mid2:subStop] %>% mean(.)
    
    IS <- ((A + B) / 2) / C #Insulation Score
    
    return(IS)
  }
  
}


#####################################################################################################################################################
                                                            # Global Insulation Score #
#####################################################################################################################################################

#Insulation score for observed data matrices : window is given in bp and converted in bin within the function
globalInsScore <- function(matrix, start, stop, bin.width, window = 512e3) {
  
  w <- round(window / bin.width)  # To avoid heavy formula depending on the window
  first <- start / bin.width + w  # = 1st bin + window
  last <- stop / bin.width - w  # = last bin - window
  
  # Creation of a 4-columns data frame: chr, location of the beginning and the end of the bin and insulation score
  IS <- sapply (first:last, function(i) {
    mean(matrix[(i-w):(i-1), (i+1):(i+w)])
  })
  
  startBin <- first:last * bin.width - bin.width  # Location of the beginning of the bin
  stopBin <- first:last * bin.width  # Location of the end of the bin
  IS.df <- data.frame(chr <- "1", startBin <- startBin, stopBin <- stopBin, IS <- IS)  # Creation of the data frame
  
  return(IS.df)
}


#####################################################################################################################################################
                                                        # Local Structural Impact Score #
#####################################################################################################################################################

localSIC <- function(mutated.matrix, wildtype.matrix, mutation.start, bin.width, window = 4) {
  # window: window where we calculate the SIC (i.e. SIC at a local scale). Number given in bin.
  
  # Conversion of the location from bp to bin
  i = round(mutation.start/ bin.width) + 1
  
  # Change the name of "window" to avoid heavy formulas
  w <- window
  
  if ((i - w) < 1 | (i + w) > 250) {
    
    return(NA)
    
  } else {
    
    # Extract the data from the matrices to the window
    mut = mutated.matrix[(i-w):(i-1),(i+1):(i+w)]
    wt = wildtype.matrix[(i-w):(i-1),(i+1):(i+w)]
    
    # Calculate the ratio mut/wt. If there is no value, return NA
    ratio = ifelse(is.nan(mut) | is.nan(wt), NaN, mut/wt)
    ratio = matrix(ratio, nrow = nrow(mut), ncol = ncol(mut))
    
    # Calculate the log2 of the ratio for each value in the "ratio matrix", then get its abs, and then calculate the mean of each value
    SIC = mean(abs(log2(ratio)))
    
    return(SIC)
  }
  
}


#####################################################################################################################################################
                                                        # Global Structural Impact Score #
#####################################################################################################################################################

StructImpact <- function(mutated.matrix, wildtype.matrix, mutation.start, bin.width) {
  
  # Conversion of the location from bp to bin
  i = round(mutation.start / bin.width) + 1
  
  # Extract the data from the matrices to the window
  mut = mutated.matrix[i,]
  wt = wildtype.matrix[i,]
  
  # Calculate the ratio mut/wt. If there is no value, return NA
  ratio = ifelse(is.nan(mut) | is.nan(wt), NaN, mut/wt)
  
  # Calculate the log2 of the ratio for each value in the "ratio matrix", then get its abs, and then calculate the mean of each value
  SIC = mean(abs(log2(ratio)), na.rm = TRUE)
  
  return(SIC)
}
