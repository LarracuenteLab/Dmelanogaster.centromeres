#############
#	R script to make Figure 3B of Chang et al. 2019 PLoS Biology
# Maximum likelihood phylogenetic tree based on the entire sequence of all
# G2/Jockey-3 copies in D. melanogaster inside (squares) and outside (circles) of 
# centromeric contigs, and on the consensus repeat in its sister species D. sechellia 
# and D. simulans, and a more distantly related species (D. yakuba). 
# The tree shows that centromeric G2/Jockey-3 elements do not have a single origin.
# By Xiaolu Wei
# Xiaolu_Wei@URMC.Rochester.edu
############## 

library(ape)

clust <- read.table("G2_and_Jockey-3_for_R.txt",sep='\t',header=T,fill=T)
tips <- clust[,3]
names <- clust[,4]
category <- clust[,5]

pdf("G2_Jockey-3_tree.pdf", width=9, height=6, useDingbats=FALSE)

MyTree <- read.tree("RAxML_bipartitionsBranchLabels.consensus_1_bigelements1000_with_outgroup_Jockey-3_from_Yak_and_Sec.automre")
tree.bionj <- root(MyTree, "Jockey_3_Dyak")
plot(tree.bionj,show.tip.label=FALSE,no.margin=TRUE,type="fan",open.angle=10)
tiplabels(cex=2, pch=category, col=tips, bg="white")
tiplabels(names,frame="none",cex=1.6,adj=0)
add.scale.bar(cex = 1, font = 2, col = "black")
legend("topright", legend=c("G2","Jockey-3","simulans","sechellia"), col=c(4,2,5,3), lwd=0, cex=1.5, pch=c(16,16,16,16), bty='n')

dev.off()
