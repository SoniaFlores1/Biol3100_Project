library(tidyverse)
library(ape)
library(phangorn)
library(Biostrings)
library(ggtree)
library(ShortRead)


#Converting alignments to phyDat####
alignM16sSet<- readDNAMultipleAlignment("./Data/Cleaned/alignmentM_16s.fasta")
print(alignM16sSet)

align_phyDatM<-msa::msaConvert(x=alignM16sSet, type= "phangorn::phyDat")

phydist<-dist.ml(align_phyDatM)


nj<- phangorn::NJ(phydist)

ape::is.rooted(nj)
root_nj<-ape::root(phy = nj,"Didelphis virginiana_51",
                   resolve.root=TRUE)

plot(root_nj)
#Build tree####
fit <- pml(root_nj, align_phyDatM)
fit <- optim.pml(fit, rearrangement="NNI")

bs <- bootstrap.pml(fit, bs=100, optNni=TRUE)
treeBS <- plotBS(fit$tree, bs, type ="phylogram")

root_bs<- ape::root(phy = treeBS,"Didelphis virginiana_51", resolve.root=TRUE)
plotBS(root_bs, type="phylogram")


gtt2<-ggtree(root_bs, mapping=NULL, branch.length="none" ,root.position = -5)+
  geom_tiplab(alignt=FALSE, size=3)+
  ggtree::geom_nodelab(node= 'internal')
gtt2

tree2 <- groupClade(root_bs, c(83,55,90))
p <- ggtree(tree2, aes(color=group), branch.length="none") + 
  scale_color_manual(values=c("black","firebrick", "steelblue", "orchid"))+
  geom_tiplab(alignt=FALSE, size=3)+
  geom_treescale(fontsize=3)+
  ggtree::geom_nodelab(node= 'internal')
p

ggsave(filename="16s_tree.png", path ="./Images", dpi="print" ,
       width=20, height=15)


#Build COI tree####
alignMCOISet<- readDNAMultipleAlignment("./Data/Cleaned/alignmentM_COI.fasta")
print(alignMCOISet)

align_phyDatM<-msa::msaConvert(x=alignMCOISet, type= "phangorn::phyDat")

phydist<-dist.ml(align_phyDatM)


nj<- phangorn::NJ(phydist)

ape::is.rooted(nj)
root_nj<-ape::root(phy = nj,"Didelphis virginiana_41",
                   resolve.root=TRUE)

plot(root_nj)

fit <- pml(root_nj, align_phyDatM)
fit <- optim.pml(fit, rearrangement="NNI")

bs <- bootstrap.pml(fit, bs=100, optNni=TRUE)
treeBS <- plotBS(fit$tree, bs, type ="phylogram")

root_bs<- ape::root(phy = treeBS,"Didelphis virginiana_41", resolve.root=TRUE)
plotBS(root_bs, type="phylogram")

gtt2<-ggtree(root_bs, mapping=NULL, branch.length="none" ,root.position = -1)+
  geom_tiplab(alignt=FALSE, size=3)+
  geom_treescale(fontsize=3)+
  ggtree::geom_nodelab(node= 'internal')

gtt2

tree2 <- groupClade(root_bs, c(70, 58, 75,78, 57))
p <- ggtree(tree2, aes(color=group), branch.length="none") + 
  scale_color_manual(values=c("black","firebrick", "steelblue", "orchid", "black","firebrick"))+
  geom_tiplab(alignt=FALSE, size=3)+
  geom_treescale(fontsize=3)+
  ggtree::geom_nodelab(node= 'internal')
p

ggsave(filename="COI_tree.png", path ="./Images", dpi="print" ,
       width=20, height=15)

#Build 18s Tree####


alignM18sSet<- readDNAMultipleAlignment("./Data/Cleaned/alignmentM_18s.fasta")
print(alignM18sSet)

align_phyDatM<-msa::msaConvert(x=alignM18sSet, type= "phangorn::phyDat")

phydist<-dist.ml(align_phyDatM)


nj<- phangorn::NJ(phydist)

ape::is.rooted(nj)
root_nj<-ape::root(phy = nj,"Didelphis virginiana_37",
                   resolve.root=TRUE)

plot(root_nj)

fit <- pml(nj, align_phyDatM)
fit <- optim.pml(fit, rearrangement="NNI")

bs <- bootstrap.pml(fit, bs=100, optNni=TRUE)
treeBS <- plotBS(fit$tree, bs, type ="phylogram")

root_bs<- ape::root(phy = treeBS,"Didelphis virginiana_37", resolve.root=TRUE)
plotBS(root_bs, type="phylogram")

gtt2<-ggtree(root_bs, mapping=NULL, branch.length="none" ,root.position = -1)+
  geom_tiplab(alignt=FALSE, size=3)+
  geom_treescale(fontsize=3)+
  ggtree::geom_nodelab(node= 'internal')

gtt2

tree2 <- groupClade(root_bs, c(66, 51, 62))
p <- ggtree(tree2, aes(color=group), branch.length="none") + 
  scale_color_manual(values=c("black", "black", "steelblue", "firebrick","orchid","firebrick"))+
  geom_tiplab(alignt=FALSE, size=3)+
  geom_treescale(fontsize=3)+
  ggtree::geom_nodelab(node= 'internal')
p

ggsave(filename="18s_tree.png", path ="./Images", dpi="print" ,
       width=20, height=15)


