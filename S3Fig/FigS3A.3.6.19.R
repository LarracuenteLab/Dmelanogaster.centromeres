#############
# R script to make Figure S3A of Chang et al. 2019 PLoS Biology
#   Reproducibility of CENP-A ChIP enrichment among replicates in embryos and S2 cells.
# Locations of the top 100 strongest peaks for each ChIP experiment. A) Plot of the 
# location of top 100 strongest peaks for each ChIP experiment on the diagonal
# (see details in Table S4). For the four replicate ChIP experiment in our OreR embryos,
# we examined the reproducibility of our experiments by first applying the IDR 
# (Irreproducible Discovery Rate) test and only keeping peaks with IDR≤0.05. 
# The number of these peaks are plotted below the diagonal. 
# Between Replicates 2 and 3, we found a total of 16,870 overlapping peaks, 
# but 16,833 were weakly enriched relative to the overlapping peaks between other 
# datasets because they are technical repeats with a shared library bias 
# (Accel, see Supplemental Methods). We therefore only report the 37 strongest peaks 
# (the average peak number of other comparisons between replicates). 
# The IDR dataset comparisons are in Table S5. We show the correlation between the 
# CENP-A ChIP replicates above the diagonal. Plotted are the signal strength after
# IDR tests (normalized ChIP over input ratio from 1-1,000 on a log10 scale) with 
# Spearman’s rho. The five contigs with the most consistent peaks within and among 
# replicates correspond to the five centromeric candidates.  
#
# By Ching-Ho Chang
# cchang45@ur.rochester.edu
#############
library(ggplot2)
idr.peak<-read.table("idr_test_for_fig",header=T)

idr.peak$Contig <- factor(idr.peak$Contig, levels=c("Y_Contig26","Contig79","Contig119","tig00057289","3R_5","Others"))
idr.peak$identity <- factor(idr.peak$identity, levels=c("Y_Contig26","Contig79","Contig119","tig00057289","3R_5","Others","Total"))
p<-ggplot(idr.peak,aes(y=num,x=Contig,fill=identity))+ geom_bar(stat = "identity",aes(),position = position_stack(reverse = TRUE))+ facet_grid(samp2~samp1, switch = "both",drop = FALSE)+scale_x_discrete(breaks=NULL)+ylab("# of peaks")+xlab("")+scale_fill_manual(guide=FALSE,values=c("#999999","#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2","#000000"))+ theme_light()

