library(tidyverse)
#Install BiocManager Package. Available from Cran.
#Also install the Apes and phangorn package from Cran.

library(BiocManager)
library(ape)
library(phangorn)
library(taxonomizr)
#Now download the ggtree, msa, and Biostrings package from Bioconductor.
BiocManager::install("ggtree")
BiocManager::install("Biostrings")
BiocManager::install("msa")
BiocManager::install("ggmsa")

library(Biostrings)
library(ggtree)
library(msa)
library(ggmsa)


read.FASTA(file ="./Data/Raw/archosauria_18s_rrna.fasta")
read.FASTA(file ="./Data/Raw/musophagidae_12s_trnaval_16s.fasta")
read.FASTA(file ="./Data/Raw/neognathae_12s_trnaval_16s_trnaleu.fasta")

Biostrings::readDNAStringSet(filepath = "./Data/Raw/archosauria_18s_rrna.fasta")

#Reading in sequences as a DNAStringSet####
filepath<- "./Data/Raw/archosauria_18s_rrna.fasta"
filepath2<- "./Data/Raw/musophagidae_12s_trnaval_16s.fasta"
filepath3<- "./Data/Raw/neognathae_12s_trnaval_16s_trnaleu.fasta"
#mySeqs<- readDNAStringSet(filepath)

mySeqs2<-readDNAStringSet(c(filepath, filepath2, filepath3))
mySeqs2

#obtaining taxonomic names from accession numbers
prepareDatabase('accessionTaxa.sql')


accessions_nums<- mySeqs2 %>% names() %>% 
str_split(" ") %>% 
  map_chr(1)

taxonomizr::accessionToTaxa(accessions = accession_nums, 'accessionTaxa.sql')

#sequence alignment (ClustalW and Muscle)####
#alignmentW <- msa::msaClustalW(inputSeqs =mySeqs, cluster = "nj",type = "dna")
#alignmentM<- msa::msaMuscle(inputSeqs = mySeqs, cluster = "neighborjoining", type = "dna")

alignmentW2<- msaClustalW(inputSeqs = mySeqs2, cluster= "nj",type="dna")
#alignmentM2<- msaMuscle(inputSeqs= mySeqs2, cluster = "upgmb", type= "dna",
                        #maxiters=5)
print(alignmentW2, show="complete")
ape_alignW2<-msa::msaConvert(x = alignmentW2,type="ape::DNAbin")



tidyW2<-tidy_msa(msa=ape_alignW2, start = 10, end= 70)
ggplot()+ggmsa::geom_msa(data=tidy, font=NULL)

#Conservation scores####
data("BLOSUM62")
ClustalWScore<-msa::msaConservationScore(alignment, BLOSUM62)

#converting the msa alignments to other files####
align_phyDat<-msa::msaConvert(x=alignmentW, type = "phangorn::phyDat")
align_phyApe<- msa::msaConvert(x = alignmentW, type="ape::DNAbin")

#In Phangorn now, making a distance matrix (Maximum Likelihood)####
phydist <- dist.ml(align_phyDat)

#applying neighborjoining ####
nj <- phangorn::NJ(phydist)

#plotting a tree####
plot(nj, main="nj")

#messing with ggtrees####
ggt<-ggtree(nj, cex = 1, aes(color=branch.length))+
  scale_color_continuous(high='lightskyblue1',low='coral4')+
  geom_tiplab(align=TRUE, size=2)+
  geom_treescale(y = - 4, color = "coral4", fontsize = 5)

njmsaplot<-msaplot(ggt, align_phyApe, offset = 0.009, width=1, height = 0.5)
njmsaplot

#experimenting with ape.####
ape::dist.dna(align_phyApe, model="JC69", as.matrix = TRUE )


