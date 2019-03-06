#############
# R script to make Figure S6 of Chang et al. 2019 PLoS Biology
# Figure S6 is modified from the R plot generated from this script in Illustrator. 
#
# Relative depth of Pacbio reads across centromeric contigs. Pacbio reads were mapped
# to the genome using Minimap (v 2.11) and the parameter “-ax map-pb.”. Shown are A) 
# X centromere, B) centromere 4; C) Y centromere; D) centromere 3; E) centromere 2.
# The depth of only the high-quality mapped reads (mapped Q ≥ 30) was estimated for 
# each position and normalized by the median depth of other genomic regions (98.32x for
# autosomes and 49.16x for sex chromosomes) to get relative depth. The relative depth 
# of the TE-rich islands are close to 1, whereas the depth of the flanking simple 
# satellites are uneven with some regions >1 and some <1. We therefore exclude simple
# repeats from any assembly-based analyses and color these regions gray in Fig. 2 and 
# Fig. S4 to indicate that caution should be used in interpreting these regions of the
# assembly.
#
# By Ching-Ho Chang
# cchang45@ur.rochester.edu
#############
#https://bioconductor.org/packages/devel/bioc/manuals/karyoploteR/man/karyoploteR.pdf
#https://bernatgel.github.io/karyoploter_tutorial/
library(karyoploteR)
library(regioneR)
library(GenomicRanges) 
library(rtracklayer) 
library(IRanges) 

#enrichment only based on mapq>30
all<-read.table("dmel_V5_pacbio.bam_q30.out",header=T)

#subset data frames for CENP-A and H3
c119t<-subset(all,all$Contig =='Contig119')
c79t<-subset(all,all$Contig =='Contig79')
c142t<-subset(all,all$Contig =='Contig142')
c3t<-subset(all,all$Contig =='3R_5')
cY26t<-subset(all,all$Contig =='Y_Contig26')
c22795t<-subset(all,all$Contig =='tig00022795')
c57289t<-subset(all,all$Contig =='tig00057289')

custom.genome <- toGRanges("mygenome.txt")
custom.cytobands <- toGRanges("TE_masker_cen_noAATAT.new.colors")

cen3oligo <- toGRanges("oligopaints.info/Dmel_Cen3_Oligopaints.rename.bed")
cenXoligo <- toGRanges("oligopaints.info/Dmel_CenX_30nt_Oligopaints.bed")
cen4oligo <- toGRanges("oligopaints.info/Dmel_CenUnk_Oligopaints.bed")
cenYoligo <- toGRanges("oligopaints.info/Dmel_CenY_Oligopaints.bed")


# Set colors based on annotations no AATAT and AATATAT
v = c("darkcyan" , "blue"    , "springgreen3",  "cyan",  "blueviolet" , "goldenrod1" ,  "red", "yellow" , "maroon","dodgerblue3","forestgreen","grey","grey" )
names(v) =c( "10bp",  "AAGAG", "AAGAT" , "AATAG", "Dodeca", "Gypsy",  "Other_DNA", "IGS",   "Jockey","Other_LTR","Other_Non-LTR","Other_simple","Other")


#depth
kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands,plot.type=2,ideogram.plotter = NULL,cex=.7,chromosomes=c("Contig79","Contig119","Y_Contig26","3R_5","tig00057289"))
kpAddCytobandsAsLine(kp,color.table=v,lwd=5)
kpAddBaseNumbers(kp,tick.dist=10000)
kpAxis(kp, r0=0, r1=1,tick.pos = c(0,1,2,3),ymin=0,ymax=3,cex=.5)
kpBars(kp, chr="Contig119",border=NA,col=1, x0=c119t$pos,x1=c119t$pos+1, y1=c119t$depth/98.32,data.panel = 1, r0=0,r1=1/3)
kpLines(kp,chr="Contig119",x=1:93914, y=1/3,col="blue")
kpBars(kp, chr="Contig79",border=NA,col=1, x0=c79t$pos,x1=c79t$pos+1, y1=c79t$depth/49.16,data.panel = 1, r0=0,r1=1/3)
kpLines(kp,chr="Contig79",x=1:70181, y=1/3,col="blue")
kpBars(kp, chr="Y_Contig26",border=NA,col=1, x0=cY26t$pos,x1=cY26t$pos+1, y1=cY26t$depth/49.16,data.panel = 1, r0=0,r1=1/3)
kpLines(kp,chr="Y_Contig26",x=1:139957, y=1/3,col="blue")
kpBars(kp, chr="3R_5",border=NA,col=1,  x0=c3t$pos,x1=c3t$pos+1, y1=c3t$depth/98.32,data.panel = 1, r0=0,r1=1/3)
kpLines(kp,chr="3R_5",x=1:103827, y=1/3,col="blue")
kpBars(kp, chr="tig00057289",border=NA,col=1,  x0=c57289t$pos,x1=c57289t$pos+1, y1=c57289t$depth/98.32,data.panel = 1, r0=0,r1=1/3)
kpLines(kp,chr="tig00057289",x=1:24561, y=1/3,col="blue")

dev.copy(pdf,"Figure_2_grey_R1_no_AATAT_small.pdf", width=12, height=10)
dev.off()

legend("topright",c( "10bp",  "AAGAG", "AAGAT" , "AATAG", "AATAT",  "AATATAT", "Dodeca", "Gypsy",  "Other_DNA", "IGS",   "Jockey","Other_LTR","Other_Non-LTR","Other_simple","Other"),col=c("darkcyan" , "blue"    , "springgreen3",  "cyan",    "black"  ,  "black", "blueviolet" , "goldenrod1" ,  "red", "yellow" , "maroon","dodgerblue3","forestgreen","grey","grey" ),pch=15)
