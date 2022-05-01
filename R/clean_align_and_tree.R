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

#Converting alignments to phyDat####
align_phyDatM<-msa::msaConvert(x=alignmentM, type= "phangorn::phyDat")

phydist<-dist.ml(align_phyDatM)


nj<- phangorn::NJ(phydist)

ape::is.rooted(nj)
root_nj<-ape::root(phy = nj,"Didelphis virginiana_37",
                   resolve.root=TRUE)

plot(root_nj)
#Build tree####
gtt2<-ggtree(root_nj, mapping=NULL, branch.length="none" ,root.position = -1)+
  geom_tiplab(alignt=FALSE, size=3)+
  geom_treescale(fontsize=3)+
  ggtree::geom_nodelab(node= 'internal')
gtt2


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

align_phyDatM<-msa::msaConvert(x=alignmentM, type= "phangorn::phyDat")

phydist<-dist.ml(align_phyDatM)


nj<- phangorn::NJ(phydist)

ape::is.rooted(nj)
root_nj<-ape::root(phy = nj,"Didelphis virginiana_41",
                   resolve.root=TRUE)

plot(root_nj)

gtt2<-ggtree(root_nj, mapping=NULL, branch.length="none" ,root.position = -1)+
  geom_tiplab(alignt=FALSE, size=3)+
  geom_treescale(fontsize=3)+
  ggtree::geom_nodelab(node= 'internal')

gtt2

