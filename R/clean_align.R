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

f1<-"./Data/Raw/archosauria_18s.fasta"
f2<-"./Data/Raw/other_18s.fasta"


#Cleaning and saving sequences####
seqs<-ShortRead::readFasta(c(f1,f2))

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

ShortRead::writeFasta(object = seqs, file = "./Data/Cleaned/sequences_18s.fasta", mode ="w")
sequences<- Biostrings::readDNAStringSet("./Data/Cleaned/sequences_18s.fasta")


sequences_s<- ShortRead::readFasta("./Data/Cleaned/sequences_18s.fasta")
newids <- BStringSet(paste0(sequences_s@id, "_",1:length(sequences_s)))
sequences_s@id <- newids
sequences_s@id

writeFasta(sequences_s,"./Data/Cleaned/sequences_nonredundant18s.fasta")

#Aligning sequences using msaMuscle####
cleanseq<- Biostrings::readDNAStringSet("./Data/Cleaned/sequences_nonredundant18s.fasta")

alignmentM<-msa::msaMuscle(inputSeqs=cleanseq, cluster="neighborjoining", type="dna", order="input")
print(alignmentM, show="complete")

ShortRead::writeFasta(DNAStringSet(alignmentM),"./Data/Cleaned/alignmentM_18s.fasta")

#creating alignment graphic####
alignM18s<- read.FASTA("./Data/Cleaned/alignmentM_18s.fasta")

alignM18sSet<- readDNAMultipleAlignment("./Data/Cleaned/alignmentM_18s.fasta")
print(alignM18sSet)

ggmsa(alignM18s, 225, 350, color = "Chemistry_NT", font = "DroidSansMono", 
      char_width = 0.7, seq_name = TRUE, border=NA, by_conservation=TRUE)

ggsave(filename="18smusclealign_img.png", path ="./Images", dpi="print" ,
       width=11, height=4)


#All that again, but for the COI sequences####
f1<-"./Data/Raw/archosauria_COI.fasta"
f2<-"./Data/Raw/other_COI.fasta"

seqs<-ShortRead::readFasta(c(f1,f2))

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

ShortRead::writeFasta(object = seqs, file = "./Data/Cleaned/sequences_COI.fasta", mode ="w")
sequences<- Biostrings::readDNAStringSet("./Data/Cleaned/sequences_COI.fasta")


sequences_s<- ShortRead::readFasta("./Data/Cleaned/sequences_COI.fasta")
newids <- BStringSet(paste0(sequences_s@id, "_",1:length(sequences_s)))
sequences_s@id <- newids
sequences_s@id

writeFasta(sequences_s,"./Data/Cleaned/sequences_nonredundantCOI.fasta")

cleanseq<- Biostrings::readDNAStringSet("./Data/Cleaned/sequences_nonredundantCOI.fasta")

alignmentM<-msa::msaMuscle(inputSeqs=cleanseq, cluster="neighborjoining", type="dna", order="input")
print(alignmentM, show="complete")

ShortRead::writeFasta(DNAStringSet(alignmentM),"./Data/Cleaned/alignmentM_COI.fasta")

alignMCOI<- read.FASTA("./Data/Cleaned/alignmentM_COI.fasta")

alignMCOISet<- readDNAMultipleAlignment("./Data/Cleaned/alignmentM_COI.fasta")
print(alignMCOISet)

ggmsa(alignMCOI, 250, 375, color = "Chemistry_NT", font = "DroidSansMono", 
      char_width = 0.7, seq_name = TRUE, border=NA, by_conservation=TRUE)

ggsave(filename="COImusclealign_img.png", path ="./Images", dpi="print" ,
       width=11, height=4)

#16s Alignment####
alignM16s<- read.FASTA("./Data/Cleaned/alignmentM_16s.fasta")

alignM16sSet<- readDNAMultipleAlignment("./Data/Cleaned/alignmentM_16s.fasta")
print(alignM16sSet)

ggmsa(alignM16s, 250, 375, color = "Chemistry_NT", font = "DroidSansMono", 
      char_width = 0.7, seq_name = TRUE, border=NA, by_conservation=TRUE)

ggsave(filename="16smusclealign_img.png", path ="./Images", dpi="print" ,
       width=11, height=4)
