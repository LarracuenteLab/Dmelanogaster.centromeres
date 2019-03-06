#############
# R script to make Figure S9A of Chang et al. 2019 PLoS Biology
# Figure S9 is modified from the R plot generated from this script in Illustrator. 
# Transcription of G2/Jockey-3 elements. A) Shown is the plot of the normalized reads 
# depth from uniquely mapped reads (mapping quality ≥ 10) across the G2/Jockey-3 
# consensus element obtained from mapping total and poly-A RNA-seq data from testes 
# [15, 16] to our repeat library.  
#
# By Ching-Ho Chang
# cchang45@ur.rochester.edu
#############
data<-read.table("Expression_G2.out",header=T)
# normalization by mapped reads after rRNA in silico removement
data$bdgp_1=data$bdgp_1/392566684 * 10^6
data$bdgp_2=data$bdgp_2/440701989 * 10^6
data$wt_1=data$wt_1/10316907 * 10^6
data$wt_2=data$wt_2/10999567 * 10^6
pdf("FigS9_G2_read_distribution.pdf",8,4)
plot(data$position,data$bdgp_1,col=1,type='l',xlab='G2 position',ylab='Depth/Total # of reads * 10^6',ylim=c(0,0.5))
lines(data$position,data$bdgp_2,col="orange",lwd=2)
lines(data$position,data$wt_1,col="blue",lwd=2)
lines(data$position,data$wt_2,col="forestgreen",lwd=2)
legend("topright", legend=c("total_RNA1","total_RNA2","polyA_RNA1","polyA_RNA2"), col=c(1,"orange","blue","forestgreen"),lty=c(1,1,1,1),lwd=2,cex=0.8,bty='n')
dev.off()