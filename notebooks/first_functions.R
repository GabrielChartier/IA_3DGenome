# Library needed

library(ggplot2)
library(dplyr)
library(TADkitdev2)
library(BiocGenerics)
library(GenomicRanges)


# Insulation score of a viewpoint

insScore <- function(matrix, view.point, window = 256e3, bin.width = 64e3) {
#viewpoint in bp; window: default value = 256e3: can be modified; bin.width (=resolution): default value = 64e3: can be modified
  
  vp_bin = view.point/bin.width #conversion of the viewpoint from bp to bin
  start = vp_bin-(window/bin.width)+1 #position of the beginning of the viewpoint (in bin)
  stop = vp_bin+(window/bin.width) #position of the end of the viewpoint (in bin)

  sub_mat = matrix[start:stop, start:stop] #viewpoint converted as a "sub"-matrix

  sub_mat_sym = as.matrix(sub_mat)
  sub_mat_sym[lower.tri(sub_mat_sym)] <- t(sub_mat_sym)[lower.tri(sub_mat_sym)] #make the matrix symmetrical

  #definition of the position of A, B & C with (4 zones in the "sub"-matrix):
  #   A   C
  #   C   B
  #then calculate the mean number of interactions within each zone

  sub_start = 1
  sub_stop = nrow(sub_mat_sym)
  mid1 = nrow(sub_mat_sym)/2
  mid2 = mid1 + 1

  A = sub_mat_sym[sub_start:mid1, sub_start:mid1] %>% mean(.)

  B = sub_mat_sym[mid2:sub_stop, mid2:sub_stop] %>% mean(.)

  C = sub_mat_sym[sub_start:mid1, mid2:sub_stop] %>% mean(.)

  IS = (A+B)/2 - C #Insulation Score

  return(IS)
}


# 
