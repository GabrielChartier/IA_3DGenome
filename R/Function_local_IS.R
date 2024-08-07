# Insulation score of a viewpoint

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

