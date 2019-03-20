#############
# R script to make Figure S10 of Chang et al. 2019 PLoS Biology
# Figure S10 is modified from the R plot generated from this script in Illustrator.
# Relationship of IGS in D. melanogaster and closely related species of the simulans 
# clade (D. simulans and D. sechellia) and D. yakuba. Maximum likelihood phylogenetic 
# tree of all individual IGS sequences found in the D. melanogaster genome with related
# outgroups. Node support is only shown for key nodes in the tree (complete tree is in
# supplemental file 14). All centromeric IGS sequences appear to have a single origin:
# they duplicated from sex-linked IGS interspersed at the rDNA loci at some time near 
#  the divergence of the simulans clade and D. melanogaster. IGS repeats in blue (extra)
# are similar to the IGS at the cen 3 island, but are on small contigs of unknown 
# origin, one of which is moderately enriched in CENP-A. We refer to all black and 
# blue IGS repeats as IGS3cen. 
#
# By Ching-Ho Chang
# cchang45@ur.rochester.edu
#############
library(ape)
MyTree <- read.tree("RAxML_bipartitionsBranchLabels.igs.automre.tree")
tree.bionj <- root(MyTree, "yak_240bp")
plot(tree.bionj,show.tip.label=FALSE,no.margin=TRUE,type="fan",open.angle=10)
clust <-read.table("IGS_for_R",fill=T)
tips <- clust[,3]
names <- clust[,4]
tiplabels(cex=1,pch=16,col=tips)
 tiplabels(names,frame="none",cex=1,adj=0)
 add.scale.bar(-0.1873992,-0.3075182,cex = 0.7, font = 2, col = "black")
 legend("topright", legend=c("sechellia","simulans","Sex-linked","3cen","3 minor islands"), col=c(3,5,2,1,4), lwd=0,cex=0.8,pch=c(16,16,16,16,16,16),bty='n')
