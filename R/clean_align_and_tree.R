#Load Libraries####

library(tidyverse)
library(BiocManager)
library(ape)
library(phangorn)
library(Biostrings)
library(ggtree)
library(msa)
library(ggmsa)
library(ShortRead)
BiocManager::install("ggtree")
BiocManager::install("Biostrings")
BiocManager::install("msa")
BiocManager::install("ggmsa")
BiocManager::install("ShortRead")

# Reading in and combining Fasta files ####

f1<-"./Data/Raw/archosauria_16s.fasta"
f2<-"./Data/Raw/musophagidae_16s.fasta"
f3<-"./Data/Raw/neognathae_16s.fasta"
f4<-"./Data/Raw/other_16s.fasta"

seqs<-ShortRead::readFasta(c(f1,f2,f3,f4))

seqs@id<-
  paste0(
    seqs@id %>% 
      as.character() %>% 
      str_split(" ") %>% 
      map_chr(2),
    " ",
    seqs@id %>% 
      str_split(" ") %>% 
      map_chr(3)
  ) %>% BStringSet()

seqs@id

#Cleaning and saving sequences####

#Aligning sequences using msaMuscle####

#Converting alignments to phyDat####

#Build tree####
