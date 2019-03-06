#############
#	R script to make Figure 1B of Chang et al. 2019 PLoS Biology
#	Kseek plot showing the relative abundance of simple repeat sequences in CENP-A ChIP compared to the input. Plotted on the x-axis is the median of CENP-A ChIP reads normalized over total mapped CENP-A 
#	 ChIP reads across four ChIP replicates. Plotted on the y-axis is the median of input reads normalized over total mapped input reads across four replicates. The top 7 kmers in the ChIP read abundance are labeled.
#	 The line represents the enrichment of CENP-A ChIP/Input for AATAC, a non-centromeric simple repeat. Repeats to the right of the line are putatively enriched in CENP-A (see also Fig S1; Table S1).  
#	mapping to a 'comprehensive repeat library' with bowtie.
#   by Ching-Ho Chang
# cchang45@ur.rochester.edu
#
#############
kseek<-read.table("kseek_all_table_Fig1B.txt",header=T)
plot(kseek$IP_med,kseek$input_med,cex=0.5,pch=19,xlim=c(0,0.1),ylab=c('Input_median'),xlab=c('IP_median'))
#plot(kseek$IP_ave,kseek$input_ave,cex=0.5,pch=19,xlim=c(0,0.15),ylab=c('Input_average'),xlab=c('IP_average'))

lines(c(0,0.013102895*10),c(0,0.133418771*10)) #medians
with(kseek[1:7,], text(kseek[1:7,]$IP_med,kseek[1:7,]$input_med, labels = c("Prodsat","AAGAG","AATAG","AATAT","AATAC","AAGAT","Dodeca"), pos = 4,cex=0.8))
with(kseek[1:7,], points(kseek[1:7,]$IP_med,kseek[1:7,]$input_med, cex=1,col="blue",pch=16))
