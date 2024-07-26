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