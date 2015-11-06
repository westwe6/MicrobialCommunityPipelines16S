#Post Process 18S and 16S microbial data

rm(list=ls())

#setwd("/Users/MacbookPro/Desktop/16S")
#s16otu=read.table("stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.pick.tx.shared",header=TRUE, sep="\t",fill=T,stringsAsFactors=F)
#s16tax=read.table("stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.pick.tx.1.cons.taxonomy",header=TRUE, sep="\t",fill=T,stringsAsFactors=F)
#s1=t(s16otu)
#s1=s1[-c(1,2,3),]
#colnames(s1)=s16otu[,2]
#s16=cbind(s16tax[,c(3,1)],s1)
#colnames(s16)[c(1,2)]=c("Taxonomy","OTU")


setwd("/Users/MacbookPro/Desktop/18S")
s18otu=read.table("stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v119.wang.pick.tx.shared",header=TRUE, sep="\t",fill=T,stringsAsFactors=F)
s18tax=read.table("stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v119.wang.pick.tx.1.cons.taxonomy",header=TRUE, sep="\t",fill=T,stringsAsFactors=F)
s2=t(s18otu)
s2=s2[-c(1,2,3,239),]
colnames(s2)=s18otu[,2]
s18=cbind(s18tax[,c(3,1)],s2)
colnames(s18)[c(1,2)]=c("Taxonomy","OTU")

#write.table(s16,"Nambia_Iowa_16S_9-28-15.txt",sep="\t",row.names=FALSE)
write.table(s18,"Nambia_Iowa_18S_9-28-15.txt",sep="\t",row.names=FALSE)

setwd("/Users/MacbookPro/Desktop/16S")

# Load and view simple matrix
data16=as.matrix(read.table("Nambia_Iowa_16S_wo_taxa_92815.txt",header=TRUE,sep="\t",row.names=1))

# look at dimensions of matrix and example contents
dim(data16)
data16[1,]
data16[,3]

# Create a presence-absence transformed version of the toy matrix
dataPA=(data16>0)*1
dataPA

rich=colSums(dataPA)
rich

# Faster method using a for loop
dataREL2=data16[,1:115]*0
for(i in 1:115){
  dataREL2[,i]=as.numeric(as.character(data16[,i]))/sum(data16[,i])
}
dataREL2

colSums(dataREL2)

#transposing a matrix can be useful if running a function that requires the matrix in a particular orientation
dataREL2
t(dataREL2)



library(vegan)

# calculating pairwise distance matrix
samplePA.dist=vegdist(t(dataPA),method="jaccard")
samplePA.dist

otuPA.dist=vegdist(dataPA,method="jaccard")
otuPA.dist

sampleREL.dist=vegdist(t(dataREL2),method="bray")
sampleREL.dist

# visualization of toy matrix
samplePA.pcoa=cmdscale(samplePA.dist)
samplePA.pcoa

sampleREL.pcoa=cmdscale(sampleREL.dist)
sampleREL.pcoa

#hierarchical clustering
samplePA.clust=hclust(samplePA.dist)
sampleREL.clust=hclust(sampleREL.dist)

plot(samplePA.pcoa[,1],samplePA.pcoa[,2],cex=0,main="16S sample PCOA")
text(samplePA.pcoa[,1],samplePA.pcoa[,2],colnames(dataREL2),cex=1)
plot(sampleREL.pcoa[,1],sampleREL.pcoa[,2],cex=0,main="standardized sample PCOA")
text(sampleREL.pcoa[,1],sampleREL.pcoa[,2],colnames(dataREL2),cex=1)
plot(samplePA.clust,main="16S samples", hang=-1)
samplePA.dist
plot(sampleREL.clust,main="standardized samples", hang=-1)
sampleREL.dist

# generate a heatmap of data
heatmap(dataREL2,scale="none")








