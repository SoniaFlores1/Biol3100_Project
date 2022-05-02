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
gtt2<-ggtree(root_nj, mapping=NULL, branch.length="none" ,root.position = -1)+
  geom_tiplab(alignt=FALSE, size=3)+
  geom_treescale(fontsize=3)+
  ggtree::geom_nodelab(node= 'internal')
gtt2

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

gtt2<-ggtree(root_nj, mapping=NULL, branch.length="none" ,root.position = -1)+
  geom_tiplab(alignt=FALSE, size=3)+
  geom_treescale(fontsize=3)+
  ggtree::geom_nodelab(node= 'internal')

gtt2
