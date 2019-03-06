#############
# R script to make Figure S3B of Chang et al. 2019 PLoS Biology
# Figure S3B is modified from the R plot generated from this script in Illustrator.
# Plot of ChIP-seq data from S2 cells (this paper, Talbert et al. 2018, and Chen et al.
# 2015) and an independent embryo CID-GFP ChIP-seq dataset (see details in Table S4; 
# Talbert et al. 2018; 5m and 15m represent different MNase treatments). The centromeric
# contigs are also CENP-A enriched in these independent datasets, with the exception of
# the X chromosome centromere contig. S2 cells lack a Y and are therefore not expected
# to have peaks on the Y candidate centromere contig.
#
#
# By Ching-Ho Chang
# cchang45@ur.rochester.edu
#############
library(ggplot2)
idr.peak.supp<-read.table("idr_test_for_fig_supp",header=T)

idr.peak.supp$Contig <- factor(idr.peak.supp$Contig, levels=c("Y_Contig26","Contig79","Contig119","tig00057289","3R_5","Others"))
idr.peak.supp$identity <- factor(idr.peak.supp$identity, levels=c("Y_Contig26","Contig79","Contig119","tig00057289","3R_5","Others","Total"))
p<-ggplot(idr.peak.supp,aes(y=num,x=Contig,fill=identity))+ geom_bar(stat = "identity",aes(),position = position_stack(reverse = TRUE))+ facet_grid(samp1~samp2, switch = "both",drop = FALSE)+scale_x_discrete(breaks=NULL)+ylab("# of peaks")+xlab("")+scale_fill_manual(guide=FALSE,values=c("#999999","#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2","#000000"))+ theme_light()

