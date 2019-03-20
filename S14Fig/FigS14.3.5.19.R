
#############
#	R script to make Figure S14 of Chang et al. 2019 PLoS Biology
# Quantification of interactions between centromeres and different genomic regions 
# by Hi-C. Plots showing intra- and inter-chromosomal interactions between regions 
# in Hi-C data from: A) stage 16 embryos (end of embryogenesis) and B) embryonic cycles
# 1–8 (before zygotic genome activation; data from Ogiyama et al. 2018). The different
# colors indicate interactions with individual centromeres of all chromosomes.
# Centromere-centromere interactions are significantly more frequent than interactions
# between centromeres and distal heterochromatin, inter distal heterochromatin and 
# euchromatin, and marginally more significant than centromere-inter proximal 
# heterochromatin interactions. **** adjusted P<0.0001; * adjusted P<0.02, 
# Pairwise Wilcoxon rank sum test with false discovery rate (FDR) correction; 
# Kruskal-Wallis test by ranks with Dunn’s test for post-hoc analysis. 
#
# By Xiaolu Wei
# Xiaolu_Wei@URMC.Rochester.edu
#
#############
library(ggplot2)

count=read.table(file.choose(),header=FALSE)
summary(count)

#count=read.csv(file.choose(),header=FALSE)
#summary(count)

pdf("plot_50kbwindow_boxplot_log2_jitter_cycle1-8.pdf",width=15,height=5)
ggplot(count, aes(x=factor(V1, levels=c('centromere','proximal_heterochromatin','distal_heterochromatin','inter-proximal_heterochromatin','inter-distal_heterochromatin','euchromatin')), y=log2(V6+1) ))+
  geom_boxplot() + 
  xlab("category")+ 
  ylab("log2(count per 100kb)")+
  geom_jitter(alpha = 0.6, shape=16, size=2.5, position=position_jitter(0), aes(colour = factor(V2,levels=c("centromere_2", "centromere_3", "centromere_4", "centromere_X", "centromere_Y") )))+
  #geom_point(alpha = 0.8, shape=16, position=position_jitter(0), aes(colour = factor(V2,levels=c("centromere_2", "centromere_3", "centromere_4", "centromere_X", "centromere_Y") )))+
  scale_colour_manual(name="Centromere", values=c("orange", "green", "blue", "red", "orchid"))+
  #plot background
  theme(panel.background = element_rect(colour = "white", fill = "white")) +
  #legend key background
  theme(legend.key=element_rect(colour = "white", fill = "white")) +
  #legend text size
  theme(legend.text=element_text(size=rel(1))) +
  #size of text
  theme(axis.title=element_text(size=15))+
  theme(axis.text=element_text(colour = "black", size=10))+
  #axis line
  theme(axis.line=element_line(size=0.5, color="black"))

dev.off()

