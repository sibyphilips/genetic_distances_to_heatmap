#R
#load the package ape
require(ape)
#load the package gplots
require(gplots)
#load the package RColorBrewer
require(RColorBrewer)

#have made this script by looking at various forums and other code. So the heatmap is not mine but is a one that is used by me for my genetic distance calculations

#now read the DNA sequence

y <- read.dna("provide_your_fasta_alignment_here.fasta",format="fasta")

#calculate the raw distances with pairwise deletion
d <- dist.dna(y,model="JC69",pairwise.deletion=TRUE,as.matrix=TRUE) #I am using raw genetic distances if you want K2P please arrange accordingly
#write the distances to a csv file for future reference
#write.csv(d,file="COI_dist_raw.csv")
write.csv(d,file="CYTB_dist_JC69wwht_1.csv")#or any suitable name for your csv file

#convert the above matrix to a matrix with only 3 decimal places
D <- (matrix(as.numeric(sprintf("%.3f",d)),nrow=5))#here in the place of nrow replace the digit 5 with the number of sequences that you have in your alignment

#make the distance matrix's upper triangle NA's so that we get a plesant lowerdiagonal heatmap
D[upper.tri(D)]<-NA
#now create the heatmap

#save as svg in the respective folder
svg("provide_the_name_for_your_svg_file_here.svg")

heatmap.2(D,Rowv=F,Colv=F,distfun=dist,hclustfun=hclust,dendrogram="none",key=T,keysize=1,trace="none",margins=c(3,3),cellnote=D,notecol="black",sepwidth=c(0.5,0.5),sepcol="black",col=brewer.pal(8,"RdBu"))
#see RColorBrewer palettes to choose the needed colour spectrum for your genetic distance heatmap
dev.off()


