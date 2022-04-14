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
BiocManager::install("ShortRead")

library(Biostrings)
library(ggtree)
library(msa)
library(ggmsa)
library(ShortRead)


read.FASTA(file ="./Data/Raw/archosauria_18s_rrna.fasta")
read.FASTA(file ="./Data/Raw/musophagidae_12s_trnaval_16s.fasta")
read.FASTA(file ="./Data/Raw/neognathae_12s_trnaval_16s_trnaleu.fasta")

Biostrings::readDNAStringSet(filepath = "./Data/Raw/archosauria_18s_rrna.fasta")

#Reading in sequences as a DNAStringSet####
filepath<- "./Data/Raw/archosauria_18s_rrna.fasta"
filepath2<- "./Data/Raw/musophagidae_12s_trnaval_16s.fasta"
filepath3<- "./Data/Raw/neognathae_12s_trnaval_16s_trnaleu.fasta"
#Reading in with ShortRead
seqs<- ShortRead::readFasta(c(filepath,filepath2,filepath3))
seqs@id
#changing the sequence IDs with ShortRead
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
seqs@sread

View(seqs@id)

#Make a new Fasta file with ShortRead that includes the changes made to IDs
ShortRead::writeFasta(object = seqs, file = "./Data/Cleaned/sequences.fasta", mode ="w")
#NOTE: Mode 'w' is used to create a new Fasta file. Mode 'a" is used to append an existing Fasta.
#Now read that new file in for alignment
sequences<- Biostrings::readDNAStringSet("./Data/Cleaned/sequences.fasta")

duplicated(sequences)
#We'll need to give some unique identifiers to each sequence to avoid errors from duplicate names
sequences_s<- ShortRead::readFasta("./Data/Cleaned/sequences.fasta")

newids <- BStringSet(paste0(sequences_s@id, "_",1:length(sequences_s)))
sequences_s@id <- newids
writeFasta(sequences_s,"./Data/Cleaned/sequences_nonredundant.fasta")

cleanseq<- Biostrings::readDNAStringSet("./Data/Cleaned/sequences_nonredundant.fasta")
#mySeqs<- readDNAStringSet(filepath)

#mySeqs2<-readDNAStringSet(c(filepath, filepath2, filepath3))

#obtaining taxonomic names from accession numbers with taxonomizr####
#NOTE! This will make a database and will take some time to make
#Please make sure you have the time and internet bandwidth to do this.
#also try to save this in a more central location instead of copying it for
#each project
#prepareDatabase('accessionTaxa.sql')


#accessions_nums<- mySeqs2 %>% names() %>% 
#str_split(" ") %>% 
  #map_chr(1)

#taxonomizr::accessionToTaxa(accessions = accession_nums, 'accessionTaxa.sql')

#sequence alignment (ClustalW and Muscle)####
#alignmentW <- msa::msaClustalW(inputSeqs =mySeqs, cluster = "nj",type = "dna")
#alignmentM<- msa::msaMuscle(inputSeqs = mySeqs, cluster = "neighborjoining", type = "dna")

alignmentW<- msa::msaClustalW(inputSeqs =cleanseq, cluster = "nj", type = "dna")

alignmentM<- msa::msaMuscle(inputSeqs=cleanseq, cluster="neighborjoining", type="dna")
#alignmentW2<- msaClustalW(inputSeqs = mySeqs2, cluster= "nj",type="dna")
#alignmentM2<- msaMuscle(inputSeqs= mySeqs2, cluster = "upgmb", type= "dna",
                        #maxiters=5)

print(alignmentW, show="complete")
msa::msaConsensusSequence(x = alignmentW)

print(alignmentM, show="complete")


#Save alignment in new Fasta as a DNAStringset
ShortRead::writeFasta(DNAStringSet(alignmentW),"./Data/Cleaned/alignmentW.fasta")

ShortRead::writeFasta(DNAStringSet(alignmentM),"./Data/Cleaned/alignmentM.fasta")

#load the Fasta

alignmentWset <- DNAStringSet(alignmentW)
alignmentWset

alignmentMset<- DNAStringSet(alignmentM)
alignmentMset

#using ggmsa to display alignment in a graphical format
ggmsa::ggmsa(alignmentset, char_width=0.9, border="Blue",start = 10,end=60)
ggsave("")
?ggmsa

#ape_alignW<-msa::msaConvert(x = alignmentW,type="ape::DNAbin")



#tidyW<-tidy_msa(msa=ape_alignW, start = 15, end= 55) 
                                                    
#ggplot()+ggmsa::geom_msa(data=tidyW, font=NULL)

#Conservation scores####
#data("BLOSUM62")
#ClustalWScore<-msa::msaConservationScore(alignment, BLOSUM62)

#converting the msa alignments to other files####
align_phyDat<-msa::msaConvert(x=alignmentW, type = "phangorn::phyDat")

align_phyDatM<-msa::msaConvert(x=alignmentM, type= "phangorn::phyDat")

#In Phangorn now, making a distance matrix (Maximum Likelihood)####
phydist <- dist.ml(align_phyDatM)

#applying neighborjoining ####
nj <- phangorn::NJ(phydist)

njtrees<- bootstrap.phyDat(align_phyDatM, FUN=function(x)nj(dist.ml(x)),bs=20)
treesnj<- plotBS(tree=nj, njtrees, type = "phylogram" )

#maximum likelihood
fit<-pml(nj, align_phyDatM)
fit<- optim.pml(fit, rearrangement = "NNI")
bs<- bootstrap.pml(fit, bs=20, optNni=TRUE)
treeBS<- plotBS(tree=fit$tree, bs, type = "phylogram")

#Maximum Parsimony
treeMP<- pratchet(align_phyDatM)
treeMP<- acctran(treeMP, align_phyDatM)
BStrees<- bootstrap.phyDat(align_phyDatM, pratchet, bs=20)
treeMP2<- plotBS(treeMP, BStrees, type="phylogram")

#You can also do modelTest and bootstrap to find the best models to use for the tree
#But DO NOT RUN on SurfacePro!!!!
#mt<- modeltest(align_phyDatM)


#plotting a tree####
plot(nj, main="nj")

#messing with ggtrees####
ggt<-ggtree(nj, cex = 1, aes(color=branch.length))+
  scale_color_continuous(high='lightskyblue1',low='coral4')+
  geom_tiplab(align=TRUE, size=2)+
  geom_treescale(y = - 4, color = "coral4", fontsize = 5)

gtt2<-ggtree(nj, mapping=NULL, layout='rectangular',root.position = -4)+
  geom_tiplab(alignt=FALSE, size=3)+
  geom_treescale(y=6, fontsize=3,x = 0.5)+
  ggtree::geom_nodelab(node= 'internal')
options(ignore.negative.edge=TRUE)



#NOTE: LOOK at DECIPHER Package for multiple sequence alignments, too



