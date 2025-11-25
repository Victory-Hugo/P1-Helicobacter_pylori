# input
args <- commandArgs()
input_prefix = args[6]
##################################################################
## Finestructure R Example
## Author: Daniel Lawson (dan.lawson@bristol.ac.uk)
## For more details see www.paintmychromosomes.com ("R Library" page)
## Date: 14/02/2012
## Notes:
##    These functions are provided for help working with fineSTRUCTURE output files
## but are not a fully fledged R package for a reason: they are not robust
## and may be expected to work only in some very specific cases! USE WITH CAUTION!
## SEE FinestrictureLibrary.R FOR DETAILS OF THE FUNCTIONS
##
## Licence: GPL V3
## 
##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.

##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.

##    You should have received a copy of the GNU General Public License
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.
fslibrary = "/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/4-fineSTRCTURE/python/FinestructureLibrary.R"
source(fslibrary) # read in the R functions, which also calls the needed packages

## make some colours
some.colors<-MakeColorYRP() # these are yellow-red-purple
some.colorsEnd<-MakeColorYRP(final=c(0.2,0.2,0.2)) # as above, but with a dark grey final for capped values

### Define our input files
chunkfile<-paste(input_prefix, ".chunkcounts.out", sep='') ## chromopainter chunkcounts file
mcmcfile<-paste(input_prefix, "_mcmc.xml", sep='') ## finestructure mcmc file
treefile<-paste(input_prefix, "_tree.xml", sep='') ## finestructure tree file
treefile_wh<-paste(input_prefix, "_tree_wh.xml", sep='') ## finestructure tree file

## Additional files that you can extract from finestructure
mappopchunkfile<-paste(input_prefix, ".mutationprobs.out", sep='') # population-by-population chunkcount file for the populations used in the MAP (i.e tree)
system( paste("fs fs -X -Y -e X2",chunkfile,treefile,mappopchunkfile) )
meancoincidencefile<-paste(input_prefix, ".meancoincidence.csv", sep='') # pairwise coincidence, .i.e. proportion of MCMC files where individuals are found in the same 
system( paste("fs fs -X -Y -e meancoincidence",chunkfile,mcmcfile,meancoincidencefile) )
## there are ways of generating these within R but are either slower or more annoying - its your call how you do it

###### READ IN THE CHUNKCOUNT FILE
dataraw<-as.matrix(read.table(chunkfile,row.names=1,header=T,skip=1)) # read in the pairwise coincidence 

###### READ IN THE MCMC FILES
mcmcxml<-xmlTreeParse(mcmcfile) ## read into xml format
mcmcdata<-as.data.frame.myres(mcmcxml) ## convert this into a data frame

###### READ IN THE TREE FILES

treexml<-xmlTreeParse(treefile) ## read the tree as xml format
ttree<-extractTree(treexml) ## extract the tree into ape's phylo format
## If you dont want to plot internal node labels (i.e. MCMC posterior assignment probabilities)
## now is a good time to remove them via:
#     ttree$node.label<-NULL
## Will will instead remove "perfect" node labels
ttree$node.label[ttree$node.label=="1"] <-""
## And reduce the amount of significant digits printed:
ttree$node.label[ttree$node.label!=""] <-format(as.numeric(ttree$node.label[ttree$node.label!=""]),digits=2)

tdend<-myapetodend(ttree,factor=1) # convert to dendrogram format

####################################
## PLOT 1: RAW DENDROGRAM PLOT
pdf(file=paste(input_prefix, "FullDendrogram.pdf", sep=''),height=6,width=14)
par(mar=c(6,0,2,0),mfrow=c(1,1))
fs.plot.dendrogram(tdend,horiz=FALSE,nodePar=list(cex=0,lab.cex=0.6),edgePar=list(p.lwd=0,t.srt=90,t.off=-0.5),axes=F)
dev.off()


###########m目前不能工作
## Now we work on the MAP state
mapstate<-extractValue(treexml,"Pop") # map state as a finestructure clustering
mapstatelist<-popAsList(mapstate) # .. and as a list of individuals in populations
#mapstateunsorted<-popAsList(extractValue(treexml,"Pop")) # map state as a finestructure clustering

popnames<-lapply(mapstatelist,NameSummary) # population names IN A REVERSIBLE FORMAT (I.E LOSSLESS)
#write.csv(popnames,file=paste(input_prefix, "_pcarespopnames.csv", sep=''))
## NOTE: if your population labels don't correspond to the format we used (NAME<number>) YOU MAY HAVE TROUBLE HERE. YOU MAY NEED TO RENAME THEM INTO THIS FORM AND DEFINE YOUR POPULATION NAMES IN popnamesplot BELOW
popnamesplot<-lapply(mapstatelist,NameMoreSummary) # a nicer summary of the populations

#write.csv(popnamesplot,file=paste(input_prefix, "_pcarespopnamespopnamesplot.csv", sep=''))

names(popnames)<-popnamesplot # for nicety only
names(popnamesplot)<-popnamesplot # for nicety only


popdend<-makemydend(tdend,mapstatelist) # use NameSummary to make popdend
#popdend<-fixMidpointsComplete(popdend) # needed for obscure dendrogram reasons不用跑

popdendclear<-makemydend(tdend,mapstatelist,"NameMoreSummary")# use NameMoreSummary to make popdend
####popdendclear<-fixMidpointsComplete(popdendclear) # needed for obscure dendrogram reasons不用跑

########################
## PLOT 2: population tree
pdf(file=paste(input_prefix, "PopulationDendrogram.pdf", sep=''),height=6,width=12)
par(mar=c(8,2,2,2),cex=0.8)
fs.plot.dendrogram(popdendclear,horiz=FALSE,nodePar=list(cex=0,lab.cex=1.2,las=2),edgePar=list(p.lwd=0,t.srt=0,t.off=0.3),yaxt="n",height=0.5,dLeaf=0.2)
dev.off()

########################
## PAIRWISE COINCIDENCES

fullorder<-labels(tdend) # the order according to the tree
mcmcmatrixraw<-as.matrix(read.csv(meancoincidencefile,row.names=1)) # read in the pairwise coincidence file we created earlier
mcmcmatrix<-mcmcmatrixraw[fullorder,fullorder] 
mapstatematrix<-groupingAsMatrix(mapstatelist)[fullorder,fullorder] # map state for reference

#########################
## PLOT 3: Pairwise coincidence, showing the MAP state

source(fslibrary)
pdf(file=paste(input_prefix, "PairwiseCoincidence.pdf", sep=''),height=12,width=12)
plotFinestructure(mcmcmatrix,dimnames(mcmcmatrix)[[1]],dend=tdend,optpts=mapstatematrix,cex.axis=0.6,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=0.8))
dev.off()

pdf(file=paste(input_prefix, "PairwiseCoincidence_nooptpts.pdf", sep=''),height=12,width=12)
plotFinestructure(mcmcmatrix,dimnames(mcmcmatrix)[[1]],dend=tdend,cex.axis=0.6,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=0.8))
dev.off()

########################
## COANCESTRY MATRIX

datamatrix<-dataraw[fullorder,fullorder] # reorder the data matrix

tmatmax<-500 # cap the heatmap
tmpmat<-datamatrix 
tmpmat[tmpmat>tmatmax]<-tmatmax # 
pdf(file=paste(input_prefix, "Coancestry.pdf", sep=''),height=12,width=12)
plotFinestructure(tmpmat,dimnames(tmpmat)[[1]],dend=tdend,cols=some.colorsEnd,cex.axis=0.6,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=0.8))
dev.off()

## Population averages
popmeanmatrix<-getPopMeanMatrix(datamatrix,mapstatelist)

tmatmax<-500 # cap the heatmap
tmpmat<-popmeanmatrix
tmpmat[tmpmat>tmatmax]<-tmatmax # 
pdf(file=paste(input_prefix, "AveragedCoancestry.pdf", sep=''),height=12,width=12)
plotFinestructure(tmpmat,dimnames(tmpmat)[[1]],dend=tdend,cols=some.colorsEnd,cex.axis=0.6,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=0.8))
dev.off()


### Useful tricks with labels

mappopcorrectorder<-NameExpand(labels(popdend))
mappopsizes<-sapply(mappopcorrectorder,length)
labellocs<-PopCenters(mappopsizes)
labelcols<-c(2,2,3,3,4,4,1,1,1,5,6,6,6,7,8,8,2) # different label colours allow clearer identification of individuals, too
labelcrt=45

pdf(file=paste(input_prefix, "Coancestry2.pdf", sep=''),height=12,width=12)
plotFinestructure(tmpmat,labelsx=labels(popdendclear),labelsatx=labellocs,crt=labelcrt,dend=tdend,text.col=labelcols,cex.axis=1.0,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=0.8))
dev.off()


####################################
## PCA Principal Components Analysis
pcares<-mypca(dataraw)
# For figuring out how many PCs are important; see Lawson & Falush 2012
# You need packages GPArotation and paran
tmap<-optimalMap(dataraw)
thorn<-optimalHorn(dataraw)
c(tmap,thorn) # 11 and 5. Horn typically underestimates, Map is usually better
pcapops<-getPopIndices(rownames(dataraw),mapstatelist)
pcanames<-rownames(dataraw)
rcols<-rainbow(max(pcapops))
library(ggplot2)



write.csv(pcares$vectors,file=paste(input_prefix, "_pcares.csv", sep=''))
#mydata<-read.table("./input_prefix_pcares.csv", sep=''),header=TRUE,sep=",")
#mydata<-read.table("./", paste(input_prefix, "_pcares.csv", sep=''), sep=''),header=TRUE,sep=",")
#paste(input_prefix, ".chunkcounts.out", sep='')

#paste("mydata<-read.table("./", input_prefix, _pcares.csv", sep=''),header=TRUE,sep=","), sep='')

pdf(paste(input_prefix, "_PCA.pdf", sep=''),height=16,width=12)
par(mfrow=c(4,3))
for(i in 1:4) for(j in (i+1):5) {
  plot(pcares$vectors[,i],pcares$vectors[,j],col=rcols[pcapops],xlab=paste("PC",i),ylab=paste("PC",j),main=paste("PC",i,"vs",j),pch=rcols)
  text(pcares$vectors[,i],pcares$vectors[,j],labels=pcanames,col=rcols[pcapops],cex=0.5,pos=1)
}
dev.off()

pdf(paste(input_prefix, "_PCA2.pdf", sep=''),height=16,width=12)
par(mfrow=c(4,3))
for(i in 1:4) for(j in (i+1):5) {
  plot(pcares$vectors[,i],pcares$vectors[,j],col=rcols[pcapops],xlab=paste("PC",i),ylab=paste("PC",j),main=paste("PC",i,"vs",j))
  text(pcares$vectors[,i],pcares$vectors[,j],labels=pcanames,col=rcols[pcapops],cex=0.5,pos=1)
}
dev.off()

#write.csv(pcares$vectors,file=paste(input_prefix, "_pcares.csv", sep=''))
#mydata<-read.table("./", paste(input_prefix, "_pcares.csv", sep=''), sep=''),header=TRUE,sep=",")

##pdf(paste(input_prefix, "_PCA1.pdf", sep=''),height=16,width=12)
#ggsave(paste(input_prefix, "_PCA1.pdf", sep=''),height=16,width=12)
#par(mfrow=c(4,3))
#for(i in 1:4) for(j in (i+1):5) {
#ggplot(mydata, aes(x = pcares$vectors[,i], y = pcares$vectors[,j]))+ geom_point(size = 2, aes(shape = rcols[pcapops], color = rcols[pcapops]))+ 
#  scale_shape_manual(values = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18))+
#  scale_colour_manual(values=c("DarkRed", "DarkSalmon", "DarkSlateBlue", "DarkTurquoise", "DarkViolet", "DarkBlue", "DarkCyan", "DarkGoldenrod", "DarkGreen", "DarkKhaki", "DarkMagenta",  "DarkOrange", "DeepPink", "DimGrey", "DodgerBlue", "Firebrick", "Cyan", "DarkSlateGray", "DarkOliveGreen","DarkSeaGreen", "DeepSkyBlue", "DarkGrey", "DarkOrchid","Cornsilk"))+
 # geom_hline(yintercept=0,linetype="dashed", size=1,colour="black") + 
 # geom_vline(xintercept=0,linetype="dashed", size=1,colour="black")+theme(legend.position = "left",legend.background = element_rect(fill="white",  size=0.5, linetype="solid", colour ="black")) + 
 # theme_bw()}
 #dev.off()


#########################
## CHROMOPAINTER
## NOTE: Requires downloading "ChromoPainterExampleHGDPdata.zip" containing the additional Chromosome 1 example files
#system("ChromoPainterv2 -a 1 1 -b -in -iM -i 10 -g input_prefix_ibd.chr1.phase, sep='') -r genetic_map_Finalversion_GRCh37_chr1.recombfile -o input_prefix.chrom1.linked.hap1")

#fs pops152_4005.cp -hpc 1 -idfile pops152_4005_ibd.ids -phasefiles pops152_4005_ibd.chr{1..22}.phase -recombfiles genetic_map_Finalversion_GRCh37_chr{1..22}.recombfile -s3iters 100000 -s4iters 50000 -s1minsnps 1000 -s1indfrac 0.1 -go 
#for x in {1..22}; do
#nohup ChromoPainterv2 -g ${sample}_ibd.chr$x.phase -r genetic_map_Finalversion_GRCh37_chr$x.recombfile -t pops152_4005_ibd.ids -f pops152_4005_ibd.poplistfull 0 0 -o  ${sample}_ChromoPainterchrChr$x &
#done

#copyprobsfile<-paste(input_prefix, ".chrom1.linked.hap1.copyprobsperlocus.out.gz", sep='')
## file contains only haplotypes from individual 1 (it was run with "ChromoPainterv2 -a 1 1 -b -in -iM -i 10 -g input_prefix_ibd.chr1.phase -r genetic_map_Finalversion_GRCh37_chr1.recombfile -o input_prefix.chrom1.linked.hap1")
#myhap<-getHap(1,copyprobsfile,verbose=TRUE) # read in first haplotpe (it takes a minute or two)
#myhap2<-getHap(2,copyprobsfile,nlines=length(myhap$snps),verbose=TRUE) # second haplotype (takes half the time as provided with SNP count)
#simplecollist<-MakeColorYRP(0.1) # construct a list of colours
#cpdensityplot(myhap$snps[1:1000],myhap$probs[1:1000,],simplecollist) # plot the first 1000 SNPs
###########
# Now we will use the finestructure run to cluster the individuals

#dataraw<-as.matrix(read.table(chunkfile,row.names=1,header=T,skip=1)) # read in the pairwise coincidence 
#treexml<-xmlTreeParse(treefile) ## read the tree as xml format
#mapstate<-extractValue(treexml,"Pop") # map state as a finestructure clustering
#mapstatelist<-popAsList(mapstate) # .. and as a list of individuals in populations
## note: of course, you can just cluster by labels here if you prefer

#names(mapstatelist)<-sapply(mapstatelist,NameMoreSummary) # choose how we will name populations
#colnames(myhap$probs)<-dimnames(dataraw)[[1]] # name the probability matrix (not named otherwise, and needed for summing)
#colnames(myhap2$probs)<-dimnames(dataraw)[[1]] # name the probability matrix (not named otherwise, and needed for summing)

#popsnpmathap1<-matColSums(myhap$probs,mapstatelist) # construct a population level SNP matrix
#popsnpmathap2<-matColSums(myhap2$probs,mapstatelist) # construct a population level SNP matrix
## CAREFUL with matColSums: if the names don't match, you'll miss individuals
#collist2<-c(MakeColorYRP(0.2),rgb(0.3,0.3,0.3)) # contstruct a list of colours the correct length
#collist2<-collist2[c(1,4,7,10,13,16,2,5,8,11,14,17,3,6,9,12,15)] # reorder for higher contrast
## FINALLY MAKE THE PLOT
#png(file=paste(input_prefix, "Haplotype.png", sep=''),height=768,width=1024)
#par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)
#layout(matrix(1:2,nrow=2,ncol=1),height=c(3,1))
#cpdensityplot(myhap2$snps[1:1000],popsnpmathap2[1:1000,],collist2) # plot the first 1000 SNPs
## Make a legend
#par(mar=c(0,4,2,2)+0.1)
#plot(c(0,1),c(0,1),type="n",axes=F,xlab="",ylab="")
#legend("topleft",legend=names(mapstatelist)[1:6],col=collist2[1:6],lty=1,lwd=2,bty="n")
#legend("top",legend=names(mapstatelist)[7:12],col=collist2[7:12],lty=1,lwd=2,bty="n")
#legend("topright",legend=names(mapstatelist)[13:17],col=collist2[13:17],lty=1,lwd=2,bty="n")
#dev.off()






##################################################################
## Finestructure R Example
## Author: Daniel Lawson (dan.lawson@bristol.ac.uk)
## For more details see www.paintmychromosomes.com ("R Library" page)
## Date: 14/02/2012
## Notes:
##    These functions are provided for help working with fineSTRUCTURE output files
## but are not a fully fledged R package for a reason: they are not robust
## and may be expected to work only in some very specific cases! USE WITH CAUTION!
## SEE FinestrictureLibrary.R FOR DETAILS OF THE FUNCTIONS
##
## Licence: GPL V3
## 
##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.

##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.

##    You should have received a copy of the GNU General Public License
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.

source(fslibrary) # read in the R functions, which also calls the needed packages

## make some colours
some.colors<-MakeColorYRP() # these are yellow-red-purple
some.colorsEnd<-MakeColorYRP(final=c(0.2,0.2,0.2)) # as above, but with a dark grey final for capped values

### Define our input files
chunkfile<-paste(input_prefix, ".chunkcounts.out", sep='') ## chromopainter chunkcounts file
mcmcfile<-paste(input_prefix, "_mcmc.xml", sep='') ## finestructure mcmc file
treefile<-paste(input_prefix, "_tree.xml", sep='') ## finestructure tree file
treefile_wh<-paste(input_prefix, "_tree_wh.xml", sep='') ## finestructure tree file

###### READ IN THE CHUNKCOUNT FILE
dataraw<-as.matrix(read.table(chunkfile,row.names=1,header=T,skip=1)) # read in the pairwise coincidence 

###### READ IN THE MCMC FILES
mcmcxml<-xmlTreeParse(mcmcfile) ## read into xml format
mcmcdata<-as.data.frame.myres(mcmcxml) ## convert this into a data frame

###### DO A 'LIKELIHOOD-BASED' TREE EXTRACTION, SO THAT CUTTING THE TREE IS MEANINGFUL
## (This is similar to that done in the Leslie 2015 Peopling of the British Isles paper)
#system("fs fs -m T -k 2 paste(input_prefix, .chunkcounts.out, sep='') paste(input_prefix, _mcmc.xml, sep='') paste(input_prefix, _tree_wh.xml, sep='')")

#system("fs fs -m T -k 2 chunkfile 44pops361_mcmc.xml 44pops361_tree_wh.xml")

system( paste("fs fs -m T -k 2",chunkfile,mcmcfile,treefile_wh) )

#system( paste("fs fs -X -Y -e X2",chunkfile,treefile,mappopchunkfile) )
###### READ IN THE TREE FILES

## The whole likelihood version:
treexml_wh<-xmlTreeParse(treefile_wh) ## read the tree as xml format
ttree_wh<-extractTree(treexml_wh) ## extract the tree into ape's phylo format

ttree_wh$node.label[ttree_wh$node.label=="1"] 
ttree_wh$node.label[ttree_wh$node.label!=""] <-format(as.numeric(ttree_wh$node.label[ttree_wh$node.label!=""]),digits=2)
tdend_wh<-myapetodend(ttree_wh,factor=1) # convert to dendrogram format

## Cut the tree at the height of 200 log-likelihood units
tdend2<-cutdend(tdend_wh,200)

## And the normal version
treexml<-xmlTreeParse(treefile) ## read the tree as xml format
ttree<-extractTree(treexml) ## extract the tree into ape's phylo format
## If you dont want to plot internal node labels (i.e. MCMC posterior assignment probabilities)
## now is a good time to remove them via:
#     ttree$node.label<-NULL
## Will will instead remove "perfect" node labels
ttree$node.label[ttree$node.label=="1"] <-""
## And reduce the amount of significant digits printed:
ttree$node.label[ttree$node.label!=""] <-format(as.numeric(ttree$node.label[ttree$node.label!=""]),digits=2)

tdend<-myapetodend(ttree,factor=1) # convert to dendrogram format

####################################
## PLOT 1: RAW DENDROGRAM PLOT
pdf(file=paste(input_prefix, "_CuttingDendrograms.pdf", sep=''),height=10,width=10)
par(mfrow=c(3,1))
fs.plot.dendrogram(tdend,main="Basic Tree",nodePar=list(cex=0,lab.cex=1),edgePar=list(p.lwd=0,t.srt=0,t.off=-0.5))
fs.plot.dendrogram(tdend_wh,main="Marginal Likelihood Tree",nodePar=list(cex=0,lab.cex=1),
                   edgePar=list(p.lwd=0,t.srt=0,t.off=c(-0.5,1)))
abline(h=200,col="grey")
fs.plot.dendrogram(tdend2,main="Cut tree",nodePar=list(cex=0,lab.cex=1),edgePar=list(p.lwd=0,t.srt=0,t.off=-0.3))
dev.off()


####################################
## Now we work on the MAP state
mapstateunsorted<-popAsList(extractValue(treexml,"Pop")) # map state as a finestructure clustering

popdend<-makemydend(tdend,mapstateunsorted)
mapstatelist<-lapply(labels(popdend),getIndivsFromSummary)

popnames<-lapply(mapstatelist,NameSummary) # population names IN A REVERSIBLE FORMAT (I.E LOSSLESS)
## NOTE: if your population labels don't correspond to the format we used (NAME<number>) YOU MAY HAVE TROUBLE HERE. YOU MAY NEED TO RENAME THEM INTO THIS FORM AND DEFINE YOUR POPULATION NAMES IN popnamesplot BELOW
popnamesplot<-lapply(mapstatelist,NameMoreSummary) # a nicer summary of the populations
names(popnames)<-popnamesplot # for nicety only
names(popnamesplot)<-popnamesplot # for nicety only

## Same things for our cut state
cutstatelist<-lapply(labels(tdend2),getIndivsFromSummary)
cutnamesplot<-lapply(cutstatelist,NameMoreSummary) # a nicer summary of the populations
names(cutstatelist)<-cutnamesplot # for nicety only
names(cutnamesplot)<-cutnamesplot # for nicety only


popdend<-makemydend(tdend,mapstatelist) # use NameSummary to make popdend
popdendclear<-makemydend(tdend,mapstatelist,"NameMoreSummary")# use NameMoreSummary to make popdend

cutdendclear<-cutdend(tdend_wh,200,"NameMoreSummary")

##################################################
## That is it, now we just plot

popmatrixcut<-getPopMatrix(dataraw,cutstatelist)
popmatrix<-getPopMatrix(dataraw,mapstatelist)

pdf(file=paste(input_prefix, "_CoancestryCut.pdf", sep=''),height=12,width=12)
plotFinestructure(popmatrixcut,dimnames(popmatrixcut)[[1]],dend=cutdendclear,cols=some.colorsEnd,cex.axis=2,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=2),labmargin=15,main="Cut Dendrogram reporting average chunk counts",cex.main=2)
dev.off()

pdf(file=paste(input_prefix, "_CoancestryPop.pdf", sep=''),height=12,width=12)
plotFinestructure(popmatrix,popnamesplot,dend=popdendclear,cols=some.colorsEnd,cex.axis=2,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=0.8),labmargin=15,main="Population Dendrogram reporting average chunk counts",cex.main=2)
dev.off()

########################
## COANCESTRY MATRIX

fullorder<-labels(tdend) # the order according to the tree
datamatrix<-dataraw[fullorder,fullorder] # reorder the data matrix

tmatmax<-500 # cap the heatmap
tmpmat<-datamatrix 
tmpmat[tmpmat>tmatmax]<-tmatmax # 
pdf(file=paste(input_prefix, "_CoancestryRaw.pdf", sep=''),height=12,width=12)
plotFinestructure(tmpmat,dimnames(tmpmat)[[1]],dend=tdend,cols=some.colorsEnd,cex.axis=0.8,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=0.8),main="Raw Coancestry chunk count matrix",cex.main=2)
dev.off()

## Population averages
popmeanmatrix<-getPopMeanMatrix(datamatrix,mapstatelist)

tmatmax<-500 # cap the heatmap
tmpmat<-popmeanmatrix
tmpmat[tmpmat>tmatmax]<-tmatmax # 
pdf(file=paste(input_prefix, "_CoancestryPopAveFull.pdf", sep=''),height=12,width=12)
plotFinestructure(tmpmat,dimnames(tmpmat)[[1]],dend=tdend,cols=some.colorsEnd,cex.axis=0.6,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=0.8))
dev.off()

