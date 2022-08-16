#R
#load the package ape
require(ape)
#load the package gplots
require(gplots)
#load the package RColorBrewer
require(RColorBrewer)

#have made this script by looking at various forums and other code. So the heatmap is not mine but is a one that is used by me for my genetic distance calculations

#now read the DNA sequence
#x <- read.dna("CPCOIF.fasta",format="fasta")
y <- read.dna("/home/siby/workspace/works/papers/TOR/AUG_17/lanka/SLanka_tor.fasta",format="fasta")

#calculate the raw distances with pairwise deletion
d <- dist.dna(y,model="JC69",pairwise.deletion=TRUE,as.matrix=TRUE) #substiture x with "y" for cytb

#write the distances to a csv file for future reference
#write.csv(d,file="COI_dist_raw.csv")
write.csv(d,file="CYTB_dist_JC69wwht_1.csv")

#convert the above matrix to a matrix with only 3 decimal places
D <- (matrix(as.numeric(sprintf("%.3f",d)),nrow=5))

#make the distance matrix's upper triangle NA's so that we get a plesant lowerdiagonal heatmap
D[upper.tri(D)]<-NA
#now create the heatmap

#save as svg in the respective folder
svg("/home/siby/workspace/works/papers/TOR/AUG_17/lanka/SLanka_tor.svg")

heatmap.2(D,Rowv=F,Colv=F,distfun=dist,hclustfun=hclust,dendrogram="none",key=T,keysize=1,trace="none",margins=c(3,3),cellnote=D,notecol="black",sepwidth=c(0.5,0.5),sepcol="black",col=brewer.pal(8,"RdBu"))
dev.off()

#y <- as.matrix(data)#convert the dataframe as matrix
#Y[upper.tri(y)]<-NA #this is to convert the upper diagonal to NA's useful to draw heatmap
