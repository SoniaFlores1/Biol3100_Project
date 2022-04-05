library(tidyverse)
#Install BiocManager Package. Available from Cran.
#Also install the Apes and phangorn package from Cran.

library(BiocManager)
library(ape)
library(phangorn)
#Now download the ggtree, msa, and Biostrings package from Bioconductor.
BiocManager::install("ggtree")
BiocManager::install("Biostrings")
BiocManager::install("msa")

library(Biostrings)
library(ggtree)
library(msa)

read.FASTA(file ="./Data/Raw/archosauria_18s_rrna.fasta")

Biostrings::readDNAStringSet(filepath = "./Data/Raw/archosauria_18s_rrna.fasta")

#Reading in sequences as a DNAStringSet
filepath<- "./Data/Raw/archosauria_18s_rrna.fasta"
mySeqs<- readDNAStringSet(filepath)


#sequence alignment (ClustalW and Muscle)
alignment <- msa::msaClustalW(inputSeqs =mySeqs, cluster = "nj",type = "dna")
alignmentM<- msa::msaMuscle(inputSeqs = mySeqs, cluster = "neighborjoining", type = "dna")

#Conservation scores
data("BLOSUM62")
ClustalWScore<-msa::msaConservationScore(alignment, BLOSUM62)


align_phyDat<-msa::msaConvert(x=alignment, type = "phangorn::phyDat")

#In Phangorn now, making a distance matrix (Maximum Likelihood)
phydist <- dist.ml(align_phyDat)

#applying neighborjoining 
nj <- phangorn::NJ(phydist)

#plotting a tree
plot(nj, main="nj")


