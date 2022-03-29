library(tidyverse)
#Install BiocManager Package. Available from Cran.
#Also install the Apes and phangorn package from Cran.

library(BiocManager)
library(ape)
library(phangorn)
#Now download the ggtree, msa, andBiostrings package from Bioconductor.
BiocManager::install("ggtree")
BiocManager::install("Biostrings")
BiocManager::install("msa")

library(Biostrings)
library(ggtree)
library(msa)

read.FASTA(file ="./Data/Raw/archosauria_18s_rrna.fasta")

Biostrings::readDNAStringSet(filepath = "./Data/Raw/archosauria_18s_rrna.fasta")

