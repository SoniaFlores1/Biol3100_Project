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
read.FASTA(file="./Data/Raw/neoganathae_12s_16s_18s.fasta")
Biostrings::readDNAStringSet(filepath = "./Data/Raw/archosauria_18s_rrna.fasta")

#Reading in sequences as a DNAStringSet####
filepath<- "./Data/Raw/archosauria_18s_rrna.fasta"
filepath2<- "./Data/Raw/musophagidae_12s_trnaval_16s.fasta"
filepath3<- "./Data/Raw/neognathae_12s_trnaval_16s_trnaleu.fasta"
filepath4<- "./Data/Raw/other_12s_16s_18s.fasta"
#Reading in with ShortRead
seqs<- ShortRead::readFasta(c(filepath,filepath2,filepath3, filepath4))
seqs@id

#Cleaning up the sequences####
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

#Sequence alignment with msa####


alignmentM<- msa::msaMuscle(inputSeqs=cleanseq, cluster="neighborjoining", type="dna", order="input")
   #alignmentM_upgma<- msa:msaMuscle(inputSeqs=cleanseq, cluster ="upgma", type="dna")
   #alignmentW<- msa::msaClustalW(inputSeqs =cleanseq, cluster = "nj", type = "dna")


print(alignmentM, show="complete")


#Save alignment in new Fasta as a DNAStringset

ShortRead::writeFasta(DNAStringSet(alignmentM),"./Data/Cleaned/alignmentM.fasta")

#load the Fasta

#as a DNAStringSet
alignmentMset<- DNAStringSet(alignmentM)
alignmentMset
#Just as a Fasta
alignM<- read.FASTA("./Data/Cleaned/alignmentM.fasta")


#using ggmsa to display alignment in a graphical format####
ggmsa::ggmsa(alignmentMset, start = 25,end=70, color='Clustal') #For some reason the names aren't read this way.

#this method uses the alignment read as just a Fasta
ggmsa(alignM, 30, 70, color = "Chemistry_NT", font = "DroidSansMono", char_width = 0.6, 
      seq_name = TRUE, border=NA, by_conservation=TRUE)

ggsave(filename="musclealign_img.png", path ="./Images", dpi="print" , width=11, height=9,)


#converting the msa alignments to other files####

align_phyDatM<-msa::msaConvert(x=alignmentM, type= "phangorn::phyDat")


#In Phangorn now, making a distance matrix (Maximum Likelihood)####
phydist <- dist.ml(align_phyDatM)

#making some trees####
nj <- phangorn::NJ(phydist)
upgma<- phangorn::upgma(phydist)

#adding roots
ape::is.rooted(nj)
ape::unroot(nj)
root_nj<-ape::root(phy = nj,"Didelphis virginiana_79",
                   resolve.root=TRUE)

ape::is.rooted(upgma)
ape::unroot(upgma)
root_upgma<-ape::root(phy = upgma,"Didelphis virginiana_79",
                   resolve.root=TRUE)

is.rooted(root_upgma)



#simple plotting
plot(nj, main="nj")
plot(upgma, main="UPGMA" )
plot(root_upgma)
plot(root_upgma, main="upgma")

#Model testing for trees (do not run on Surfcace Pro!!)
mt<- modelTest(align_phyDatM, root_upgma)

mt[order(mt$AICc),]
# choose best model from the table according to AICc
bestmodel <- mt$Model[which.min(mt$AICc)]

env <- attr(mt, "env")
fitStart <- eval(get("GTR+G+I", env), env)

fit <- optim.pml(fitStart, rearrangement = "stochastic",
                 optGamma=TRUE, optInv=TRUE, model="GTR")
#using best fit model to plot the tree, with bootstrap values
bs <- bootstrap.pml(fit, bs=50, optNni=TRUE)
treeBS <- plotBS(fit$tree,bs, type ="phylogram")

root_bs<- ape::root(phy = treeBS,"Didelphis virginiana_79", resolve.root=TRUE)
plotBS(root_bs, type="phylogram")


#messing with ggtrees####



gtt2<-ggtree(root_bs, mapping=NULL, branch.length="none" ,root.position = -1)+
  geom_tiplab(alignt=FALSE, size=3)+
  geom_treescale(fontsize=3)+
  ggtree::geom_nodelab(node= 'internal')

gtt2

tree2 <- groupClade(root_bs, c(103,101))
p <- ggtree(tree2, aes(color=group), branch.length="none") + 
  scale_color_manual(values=c("firebrick", "steelblue","black"))+
  geom_tiplab(alignt=FALSE, size=3)+
  geom_treescale(fontsize=3)+
  ggtree::geom_nodelab(node= 'internal')

p
scaleClade(p, node=101, scale=.1)
  

#NOTE: LOOK at DECIPHER Package for multiple sequence alignments, too



