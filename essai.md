
### Une tentative



### Notes projet

INRAE du 24 au 26 juin, du 8 au 12 juillet, et du 22 au 26 juillet

Structural impact score (ORCA (dans le papier)) : on fait une mutation dans la ligne, et on calcule sur toute la ligne de la matrice : MT[bin mutated,]/WT[bin mutated,] %abs %mean -> valeur absolue puis moyenne du tout
- entrée : mat1 et mat2
- autre paramètres : vpstart et vpstop (viewPoint) : position en bp => il faut donc retrouver la position du bin

1e partie : spécifier le type, le prototype des fonctions
Fonction structural impact score prenant 2 matrices en entrée (+/- MT et WT) et regarder les différences -> "differencial IS"
Créer une fonction où on compare nos 2 IS : (mean(a)+mean(b))/2 - mean(c) = IS avec une matrice partagée comme ça :
  A  C
  C  B
  => le point central risque de poser problème
Faire une fonction pour calculer l'IS
Le but est de balayer toutes les matrices pour savoir quel bin a été modifié dans chacune

# Créer une fonction

=> création d'un fichier <nom>.R

fonction <- function(param1, param2, param3 = TRUE) #en entrée : utiliser 2 matrices, pas la peine d'utiliser des listes pour l'instant ; ne pas oublier bin.width pour savoir où aller chercher le viewPoint ; start et stop ; vpstart et vpstop
  #faire les contrôles (ex : matrices de même taille, etc...)
  if qqc
    stop
  si start = NULL, commencer à l.start=1
  si stop = NULL, finir à l.Stop = nrow
  utiliser la fonction modulo %/%
si on veut appeler une fonction qui n'est pas dans le package de base : librairie::fonction
penser à rendre les matrices symétriques (si on n'a qu'un triangle sur les deux)



### Tuto ORCA

---
title: "Tuto_screeningOrca"
author: "Nicolas MARY"
date: "`r format(Sys.time(), '%d %B, %Y')`"
output:
    html_document:
      number_sections: yes
      toc: yes
      toc_depth: 3
      code_folding: show
---

# Package install

This tutorial requires several packages. All packages are either available on Cran or Bioconductor. Only 1 package created specifically for this tutorial can be installed from GitHub.

First of all, to install packages from GitHub we need to install some packages:
```{r, eval=FALSE}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("biocViews")
install.packages("devtools")
```

Now you should be able to install TADkit_dev2 from GitHub with:
```{r, eval=FALSE}
devtools::install_github("Nico-FR/TADkit_dev2")
```

To learn what you can do with this package and how to do it, you must follow this tutorial:
https://github.com/Nico-FR/TADkit
To be done up to chapter 5 using bovine data.

# Load packages

Here are the packages we need. If they are not installed on your computer, you must install them first.
```{r, include=FALSE}
library(ggplot2)
library(dplyr)
library(TADkitdev2)
library(BiocGenerics)
library(GenomicRanges)
```


# Control matrices

## Load 2 matrices

```{r}
#load control matrices of bovin "3654" and "977" at 64kb resolution for chr 1 as sparse matrix
bin_width = 64e3
mat_norm.lst = list.files("./Datas_tuto/", pattern = "mapq_10.4000.mcool", full.names = T) %>%
  lapply(function(x) {cool2matrix(x, chr = "1", bin.width = bin_width, balance = T)})

# add names of each matrices
names(mat_norm.lst) = list.files("./Datas_tuto/", pattern = "mapq_10.4000.mcool", full.names = F) %>%
  gsub(".ARS-UCD1.2.mapq_10.4000.mcool", "_norm", .)

#get Observed / Expected matrix
mat_norm_OE.lst = lapply(mat_norm.lst, function(MAT) {matObsExp(MAT)})
```

## Plot matrices

```{r}
#plot 1 matrix between 1 to 11Mb
MATplot(matrix = mat_norm.lst[[1]],
        start = 5e6, stop = 15e6,
        bin.width = 64e3,
        log2 = TRUE,
        scale.color = "H",
        tad.upper.tri = "./Datas_tuto/tad3654_10kb.bed",
        tad.chr = 1)+
  ggtitle(names(mat_norm.lst[1]))

MATplot(matrix = mat_norm_OE.lst[[1]],
        start = 5e6, stop = 15e6,
        bin.width = 64e3,
        log2 = TRUE,
        scale.color = "OE",
        tad.upper.tri = "./Datas_tuto/tad3654_10kb.bed",
        tad.chr = 1)+
  ggtitle(paste0(names(mat_norm_OE.lst[1]), " Obs/Exp"))

#plot 2 matrices: 3654 & 0197 (upper and lower part of the matrix respectively)
mMATplot(matrix.upper = mat_norm_OE.lst[[1]],
         matrix.lower = mat_norm_OE.lst[[2]],
         start = 5e6, stop = 15e6, bin.width = 64e3,
         log2 = TRUE, scale.colors = "OE",
         matrix.upper.txt = "3654",
         matrix.lower.txt = "977",
         tad.upper.tri = "./Datas_tuto/tad3654_10kb.bed",
         tad.lower.tri = "./Datas_tuto/tad977_10kb.bed",
         tad.chr = 1)+
  ggtitle("Obs/Exp matrices normalized")
```

## Fold changes between 2 matrices

```{r}
#fold changes between 2 matrices
MATplot(matrix = mat_norm_OE.lst[[1]] / mat_norm_OE.lst[[2]], 
        bin.width = 64e3, start = 5e6, stop = 15e6, scale.colors = "OE2", log2 = T)+
  ggtitle(paste0("ObsExp: ",names(mat_norm_OE.lst[1]), " / ", names(mat_norm_OE.lst[2])))
```


# Compartment analysis

## datas

```{r}
chromsize = read.table("./Datas_tuto/chrom_size.tsv", col.names = c("chr", "bp"))

rnaseq.df = read.table("./Datas_tuto/salmon.merged.gene_tpm.tsv", header = TRUE)

genes.gr = read.table("./Datas_tuto/Bos_taurus.ARS-UCD1.2.104_namecorrected.gr", h = T, sep = "\t", quote = "") %>% dataframes2grange(chromsize = chromsize, strand.col = 5, name.col = 6, metadata.mcols = 8)
```

## PCA

Principal component analysis on correlation matrix with matPCA function.

```{r}
PC_norm_OE.lst = lapply(mat_norm_OE.lst, function(MAT) {
  gr = matPCA(matrix = MAT, bin.width = 64e3)
  suppressWarnings(seqlengths(gr)[1] <- chromsize$bp[1]) #les 2 lignes suivantes servent à rien apparemment (médiane = 0)
  compOrientation(gr, genes.gr, rnaseq.df)[[1]]
})
```

### Plot PC1

```{r}
PC_norm_OE.df = lapply(1:length(PC_norm_OE.lst), function(INT) {
  as.data.frame(PC_norm_OE.lst[[INT]]) %>% 
    dplyr::mutate(ID = names(PC_norm_OE.lst[INT]))
})

ggplot(data = do.call(rbind, PC_norm_OE.df), aes(y=PC1,x=start,color=ID))+geom_step()+
         xlim(1e6,24e6)
```

### Correlation between PC1

```{r}
bgCorr(bedgraph.lst = PC_norm_OE.lst)
```

## Compartment calling

```{r}
comp_norm_OE.lst = lapply(PC_norm_OE.lst, function(MAT) {
  PC1calling(MAT)
})
```

## Plot 

Compartments + PC1s

```{r}
mTADplot(tad.lst = comp_norm_OE.lst, start = 5e6, stop = 15e6, chr = 1,
         bedgraph.lst = PC_norm_OE.lst, tad.id = T)
```



# View point interaction

Now, we would like to visualize the interaction of one compartment with these neighboring regions. To do this, we'll use the function viewPointInteract that calculates the average number of interactions of a region (e.g. a compartment) along the HiC matrix.

Let's take a look at how the compartment B on chromosome 1 at 10Mb interacts :
```{r}
#compartment to analyze
viewPoint = comp_norm_OE.lst[[1]][start(comp_norm_OE.lst[[1]]) > 9.3e6][1]
viewPoint

MATplot(matrix = mat_norm_OE.lst[[1]],
        start = 5e6, stop = 15e6,
        bin.width = 64e3,
        log2 = TRUE,
        scale.color = "OE",
        tad.upper.tri = viewPoint,
        tad.lower.tri = viewPoint,
        tad.chr = 1)+
  ggtitle(paste0(names(mat_norm_OE.lst[1]), " Obs/Exp"))
```

This analysis involves dragging the window (of the compartment above) horizontally and calculating the average number of interactions.

```{r}
viewPointInteract(matrix.lst = mat_norm_OE.lst, bin.width = 64e3,
                  vp.start = 9664000, vp.stop = 10176000, #region of interest
                  start = 5e6, stop = 15e6, #area
                  log2 = TRUE, output = "plot") #return log2 of Obs/Exp count
```

If we want a metric to compare previous distributions (i.e. bedgraphs), we can calculate the correlation between them.

```{r}
viewPointInteract(matrix.lst = mat_norm_OE.lst, bin.width = 64e3,
                  vp.start = 9664000, vp.stop = 10176000, #region of interest
                  start = 5e6, stop = 15e6, #area
                  log2 = TRUE, output = "GRanges") %>%
  bgCorr(bedgraph.lst = .)
```


# Correlations between matrices


## full matrices

```{r}
matCorr(
  matrix.lst = mat_norm.lst,
  log2 = TRUE,
  output = "corr",
  self_interaction = FALSE,  max.distance = 10e6,
  bin.width = 64e3,
  method = "pearson")

matCorr(
  matrix.lst = mat_norm.lst,
  log2 = TRUE,
  output = "plot",
  self_interaction = FALSE,  max.distance = 10e6,
  bin.width = 64e3,
  method = "pearson")
```

## area

```{r}
matCorr(
  matrix.lst = mat_norm.lst,
  log2 = TRUE,
  output = "corr",
  self_interaction = FALSE,  max.distance = 10e6,
  bin.width = 64e3,
  method = "pearson",
  start = 9664000, stop = 10176000) #region of interest

matCorr(
  matrix.lst = mat_norm.lst,
  log2 = TRUE,
  output = "plot",
  self_interaction = FALSE,  max.distance = 10e6,
  bin.width = 64e3,
  method = "pearson",
  start = 9664000, stop = 10176000) #region of interest
```

# Orca matrices

Orca output 250x250 bins log(observed/expected) and log(expected) matrices at different resolutions. Here's the table of possible Orca predictions :
Par défaut : 32Mb de côté (soit 128kb/bin) ; au centre : MPOS à 16Mb
Autour du MPOS : crée une autre matrice de 250*250, cette fois de 16Mb de côté (1bin = 64kb), avec son propre MPOS, autour duquel on a aussi une matrice, etc... jusqu'à 4kb/bin
On peut déplacer le MPOS dans la matrice pour analyser et améliorer la résolution sur la région d'intérêt
Pour analyser les données, on peut créer une matrice de la taille du chr en bins, puis on pose la matrice ORCA de 250*250 sur la zone d'intérêt (le reste de la matrice du chr est nul)
```{r echo=FALSE}
data.frame(
  nb_bins = 250,
  bin_width = c(4e3, 8e3, 16e3, 32e3, 64e3, 128e3),
  predicton_size = c(4e3, 8e3, 16e3, 32e3, 64e3, 128e3) * 250
)
```

We'll therefore use the same resolution of the control matrices (i.e. 64kb) which give us the predicted matrices for a 16Mb sequence. 
```{r}
matOrca = read.table("/home/nmary/mnt/cytogene/Thomas/Orca/bos_taurus/screening_mpos_16Mb/wildtype/orca_predictions_16Mb.txt", header = FALSE) %>% as.matrix()

#plot the 250 bins
MATplot(matrix = matOrca,
        start = 1, stop = 16e6, bin.width = 64e3,
        scale.colors = "OE",
        log2 = FALSE)
```

As mentioned above, Orca output 250x250 log(observed/expected) and log(expected) matrices at different resolutions. On the other hand, our control matrices have a bin number that depends on the resolution. For a resolution of 64kb, the bin number of chromosome 1, which has a size of 158.53411 Mb, is 158.53411e6 / 64e3 = 2478.

To compare matrices, the simplest option is to create an empty matrix (of 2478 bins) into which we'll add Orca matrix at the right position (i.e. between 8 and 24Mb). To do this, we've created a function that returns the observed or observed/expected matrix:

## Observed / expected

```{r}
#load observed/expected matrix of size = 16Mb and centered at 16Mb
matOrca_OE = orca2matrix( #la fonction permettant de poser la matrice ORCA au bon endroit sur la matrice du chr
  df_prediction.path = "./Orca_Thomas/bos_taurus/screening_mpos_16Mb/wildtype/orca_predictions_16Mb.txt", 
  sep = "\t",
  mpos = 16e6, 
  chromsize = 158.53411e6, 
  output = "OE", scale = 16e6)

TADkit::MATplot(matOrca_OE, bin.width = 64e3, start = 8e6, stop = 24e6, scale.colors = "OE", log2 = T)
```

## Observed

```{r}
#load observed matrix of size = 16Mb and centered at 2Mb
matOrca = orca2matrix(
  df_prediction.path = "./Orca_Thomas/bos_taurus/screening_mpos_16Mb/wildtype/orca_predictions_16Mb.txt", 
  sep = "\t",
  mpos = 16e6, 
  chromsize = 158.534110e6, 
  output = "Obs", scale = 16e6,
  df_normmats.path = "./Orca_Thomas/bos_taurus/screening_mpos_16Mb/wildtype/orca_normmats_16Mb.txt")

TADkit::MATplot(matOrca, bin.width = 64e3, start = 8e6, stop = 24e6, scale.colors = "H", log2 = T)
```

# To do

## 1

Perform the same analysis by comparing predicted matrix (Orca wiltype) with the 2 experimental matrices:
  -plots,
  -compartments,
  -correlations,
  -viewpoint...

## 2

Perform the same analysis by comparing Orca wiltype with Orca mutated:

```{r}
list.files("./Orca_Thomas/bos_taurus/screening_mpos_10Mb/64000")
```



### Tests fonctions

---
title: "Tests fonctions"
output: html_document
date: "2024-06-20"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}
library(ggplot2)
library(dplyr)
library(TADkitdev2)
library(BiocGenerics)
library(GenomicRanges)
```

# Test

```{r}
test = function (x, y, z) {
  t = x + y + z + 10
  return(t)
}
```

# Insulation score

=> Insulation Scores are defined for a bin as an average number of interactions within
In Zhou, 2022 p. 23: IS = [mean(A)+mean(B)]/2 - mean(C) with :
  A  C
  C  B
  This square comes from the graph representing f(Predicted IS change (Mut-WT)) = Experimental IS (Mut-WT) : each zone is delimitated by the borders of the graph and straight lines passing through x=0 and y=0. The mean of each zone is calculated with the mean value of every dot inside the zone
  Or different zones of 200kb x 200kb each in the matrix ? Knowing that the size of the matrices is 1, 2, 4, 8, 16, 32, 64 Mb

## Definition of the matrix

```{r}
mat.norm.lst = list.files(path = "./Datas_tuto/", pattern = ".ARS-UCD1.2.mapq_10.4000.mcool", full.names = T) %>%
  lapply(function(MAT) {cool2matrix(MAT, bin.width = 64e3, chr = 1, balance = T)})

mat.norm.OE.lst = lapply(mat.norm.lst, function (MAT) {matObsExp(MAT)})

PC.norm.OE.lst = lapply(mat.norm.OE.lst, function(MAT) {matPCA(matrix = MAT, bin.width = 64e3, seqname = "1")})

comp.norm.OE.lst = lapply(PC.norm.OE.lst, function(MAT) {PC1calling(MAT)})

vp.interactions = viewPointInteract(matrix.lst = mat.norm.OE.lst,
                                    bin.width = 64e3,
                                    vp.start = 4992000, vp.stop = 14976000,
                                    start = 5e6, stop = 15e6,
                                    log2 = T, output = "GRange")

zone.matrix.lst = lapply()
```

