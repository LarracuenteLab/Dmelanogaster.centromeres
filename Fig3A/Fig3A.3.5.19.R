#############
#	R script to make Figure 3A of Chang et al. 2019 PLoS Biology
# Density of all repetitive elements on each candidate centromere 
# contig and the entire genome (minus the centromeres) grouped by type: 
# non-LTR retroelements, LTR retroelements, rDNA-related sequences, and 
# simple satellites. G2/Jockey-3 is present on all centromeres.
#   
# By Ching-Ho Chang
# cchang45@ur.rochester.edu
#############
library(ggplot2)
repeat.composition<-read.table("Component_for_plot",header=T)
repeat.composition$Component <- factor(repeat.composition$Component, levels=rev(c("G2/Jockey-3","Jockey-1_dsi","Doc2","Doc","G6_dm","G5_dm","G_dm","BS","Dmrt1B","TART_B1","R2","Copia","Gypsy-2_dsim","Gypsy-24_dy","Gypsy-27_dya","Gypsy-7_dse","Nomad","NTS","rDNA.ETS","rDNA.ITS2","5.8s_rDNA","28s_rDNA","18s_rDNA","AAGAT","AAGAG","AATAG","Prodsat","Dodeca")))

repeat.composition$location <- factor(repeat.composition$location, levels=c("Contig79","Contig119","Y_Contig26","3R_5","tig00057289","whole_genome"))
repeat.composition$density=round(log10(repeat.composition$density),2)
b<-ggplot(repeat.composition, aes(location,Component)) + geom_tile(aes(fill = density))+ scale_fill_gradient(low ="white", high = "blue", labels=c("1e-5","1e-4","1e-3","1e-2","1e-1"))+   scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) 
ggsave("plots.pdf",width = 8, height = 10,device = "pdf")
