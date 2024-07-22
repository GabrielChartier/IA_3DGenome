localSIC <- function(mutated.matrix, wildtype.matrix, mutation.start, bin.width, window = 4) {
  # window: window where we calculate the SIC (i.e. SIC at a local scale). Number given in bin.
  
  # Conversion of the location from bp to bin
  i = round(mutation.start / bin.width)
  
  # Change the name of "window" to avoid heavy formulas
  w <- window
  
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