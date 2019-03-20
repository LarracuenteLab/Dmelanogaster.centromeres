
#############
#	R script to make Figure 2, Figure 5 and S4 of Chang et al. 2019 PLoS Biology
# CENP-A occupies DNA sequences within putative centromere contigs. 
# Organization of each CENP-A enriched island corresponding to centromere candidates:
# A) X centromere, B) centromere 4; C) Y centromere; D) centromere 3; E) centromere 2.
# Different repeat families are color coded (see legend; note that Jockey elements 
# are shown in one color even though they are distinct elements). Shown are the 
# normalized CENP-A enrichment over input (plotted on a log scale) from one 
# replicate (Replicate 2, other replicates are in Fig. S4) colored in gray for simple
# repeats and black for complex island sequences. While the mapping quality scores 
# are high in simple repeat regions, we do not use these data to make inferences 
# about CENP-A distribution (see text for details). The coordinates of the 
# significantly CENP-A-enriched ChIPtigs mapped to these contigs (black) and the 
# predicted ChIP peaks (orange) are shown below each plot. 
#
# By Amanda Larracuente, Ching-Ho Chang
# alarracu@bio.rochester.edu; cchang45@ur.rochester.edu
############## 
#https://bioconductor.org/packages/devel/bioc/manuals/karyoploteR/man/karyoploteR.pdf
#https://bernatgel.github.io/karyoploter_tutorial/
#source("https://bioconductor.org/biocLite.R") 
#biocLite("karyoploteR") 
#biocLite("regioneR") 
library(karyoploteR)
library(regioneR)
library(GenomicRanges) 
library(rtracklayer) 
library(IRanges) 
#library(devtools)
#install_github('davetang/bedr')

all<-read.table("/enriched_all_depth_plus_md_q30.3.12.out.txt",header=T)

#subset data frames for CENP-A and H3
c119t<-subset(all,all$Contig =='Contig119')
c79t<-subset(all,all$Contig =='Contig79')
c142t<-subset(all,all$Contig =='Contig142')
c3t<-subset(all,all$Contig =='3R_5')
cY26t<-subset(all,all$Contig =='Y_Contig26')
c22795t<-subset(all,all$Contig =='tig00022795')
c57289t<-subset(all,all$Contig =='tig00057289')

cid119.R1<-subset(c119t,c119t$R1+c119t$R1_Input>4)
cid79.R1<-subset(c79t,c79t$R1+c79t$R1_Input>4)
cid142.R1<-subset(c142t,c142t$R1+c142t$R1_Input>4)
cid26.R1<-subset(cY26t,cY26t$R1+cY26t$R1_Input>4)
cid3.R1<-subset(c3t,c3t$R1+c3t$R1_Input>4)
cid22795.R1<-subset(c22795t,c22795t$R1+c22795t$R1_Input>4)
cid57289.R1<-subset(c57289t,c57289t$R1+c57289t$R1_Input>4)

cid119.R2<-subset(c119t,c119t$R2+c119t$R2_Input>4)
cid79.R2<-subset(c79t,c79t$R2+c79t$R2_Input>4)
cid142.R2<-subset(c142t,c142t$R2+c142t$R2_Input>4)
cid26.R2<-subset(cY26t,cY26t$R2+cY26t$R2_Input>4)
cid3.R2<-subset(c3t,c3t$R2+c3t$R2_Input>4)
cid22795.R2<-subset(c22795t,c22795t$R2+c22795t$R2_Input>4)
cid57289.R2<-subset(c57289t,c57289t$R2+c57289t$R2_Input>4)

cid119.R3<-subset(c119t,c119t$R3+c119t$R3_Input>4)
cid79.R3<-subset(c79t,c79t$R3+c79t$R3_Input>4)
cid142.R3<-subset(c142t,c142t$R3+c142t$R3_Input>4)
cid26.R3<-subset(cY26t,cY26t$R3+cY26t$R3_Input>4)
cid3.R3<-subset(c3t,c3t$R3+c3t$R3_Input>4)
cid22795.R3<-subset(c22795t,c22795t$R3+c22795t$R3_Input>4)
cid57289.R3<-subset(c57289t,c57289t$R3+c57289t$R3_Input>4)

cid119.R4<-subset(c119t,c119t$R4+c119t$R4_Input>4)
cid79.R4<-subset(c79t,c79t$R4+c79t$R4_Input>4)
cid142.R4<-subset(c142t,c142t$R4+c142t$R4_Input>4)
cid26.R4<-subset(cY26t,cY26t$R4+cY26t$R4_Input>4)
cid3.R4<-subset(c3t,c3t$R4+c3t$R4_Input>4)
cid22795.R4<-subset(c22795t,c22795t$R4+c22795t$R4_Input>4)
cid57289.R4<-subset(c57289t,c57289t$R4+c57289t$R4_Input>4)


R1.chip=49512826
R1.input=101043744
R2.chip=393640001
R2.input=495457706
R3.chip=63807709
R3.input=115587082
R4.chip=78211153
R4.input=43451898

R1.norm=R1.chip/R1.input
R2.norm=R2.chip/R2.input
R3.norm=R3.chip/R3.input
R4.norm=R4.chip/R4.input

custom.genome <- toGRanges("mygenome.txt")
custom.cytobands <- toGRanges("TE_masker_cen_noAATAT.new.colors")
R1_chiptigs  <- toGRanges("reformatted.chiptig/R1_denovo_peaks2.out.txt")
R1_peaks  <- toGRanges("post_filter/R1_plus_md_q30_peaks.narrowPeak")
R2_chiptigs  <- toGRanges("reformatted.chiptig/R2_denovo_peaks2.out.txt")
R2_peaks  <- toGRanges("post_filter/R2_plus_md_q30_peaks.narrowPeak")
R3_chiptigs  <- toGRanges("reformatted.chiptig/R3_denovo_peaks2.out.txt")
R3_peaks  <- toGRanges("post_filter/R3_plus_md_q30_peaks.narrowPeak")
R4_chiptigs  <- toGRanges("reformatted.chiptig/R4_denovo_peaks2.out.txt")
R4_peaks  <- toGRanges("post_filter/R4_plus_md_q30_peaks.narrowPeak")

cen3oligo <- toGRanges("~/Dropbox/mel.centromere/oligopaints.info/Dmel_Cen3_Oligopaints.rename.bed")
cenXoligo <- toGRanges("~/Dropbox/mel.centromere/oligopaints.info/Dmel_CenX_30nt_Oligopaints.bed")
cen4oligo <- toGRanges("~/Dropbox/mel.centromere/oligopaints.info/Dmel_CenUnk_Oligopaints.bed")
cenYoligo <- toGRanges("~/Dropbox/mel.centromere/oligopaints.info/Dmel_CenY_Oligopaints.bed")


# Set colors based on annotations
#no AATAT and AATATAT
v = c("darkcyan" , "blue"    , "springgreen3",  "cyan",  "blueviolet" , "goldenrod1" ,  "red", "yellow" , "maroon","dodgerblue3","forestgreen","grey","grey" )
names(v) =c( "10bp",  "AAGAG", "AAGAT" , "AATAG", "Dodeca", "Gypsy",  "Other_DNA", "IGS",   "Jockey","Other_LTR","Other_Non-LTR","Other_simple","Other")


#qPCR primer coords
primers <- toGRanges("~/Dropbox/mel.centromere/Final_assembly_0311/KaryoploteR_file/Primer.coords.5.23.18.txt")


#R1 (in Fig S4)
cid79_island.R1=subset(cid79.R1,cid79.R1$position>5655&cid79.R1$position<48432)
cid79_rep.R1=subset(cid79.R1,cid79.R1$position<=5655 | cid79.R1$position>=48432)
cid119_island.R1=subset(cid119.R1,cid119.R1$position>29252&cid119.R1$position<69146)
cid119_rep.R1=subset(cid119.R1,cid119.R1$position<=29252 | cid119.R1$position>=69146)
cid26_island.R1=subset(cid26.R1,cid26.R1$position>5655&cid26.R1$position<48432)
cid26_rep.R1=subset(cid26.R1,cid26.R1$position<=5655 | cid26.R1$position>=48432)
cid3_island.R1=subset(cid3.R1,cid3.R1$position>19152&cid3.R1$position<69541)
cid3_rep.R1=subset(cid3.R1,cid3.R1$position<=19152 | cid3.R1$position>=69541)
cid57289_island.R1=subset(cid57289.R1,cid57289.R1$position>5561&cid57289.R1$position<7864)
cid57289_rep.R1=subset(cid57289.R1,cid57289.R1$position<=5561 | cid57289.R1$position>=7864)

kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands,plot.type=2,ideogram.plotter = NULL,cex=.7,chromosomes=c("Contig79","Contig119","Y_Contig26","3R_5","tig00057289"))
kpAddCytobandsAsLine(kp,color.table=v,lwd=5)
kpAddBaseNumbers(kp,tick.dist=10000)
kpAxis(kp, r0=0.5, r1=1,tick.pos = c(0,1.5,3.0),ymin=0,ymax=3.5,cex=.5)
kpAxis(kp, r0=0.5, r1=0,tick.pos = c(0,-1.5,-3.0),ymin=0,ymax=-3.5,cex=.5)
kpBars(kp, chr="Contig119",col=1,border=NA, x0=cid119_island.R1$position,x1=cid119_island.R1$position+1, y1=log(((cid119_island.R1$R1)/(cid119_island.R1$R1_Input))/R1.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="Contig119",col="grey",border=NA, x0=cid119_rep.R1$position,x1=cid119_rep.R1$position+1, y1=log(((cid119_rep.R1$R1)/(cid119_rep.R1$R1_Input))/R1.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="Contig79",col=1,border=NA, x0=cid79_island.R1$position,x1=cid79_island.R1$position+1, y1=log(((cid79_island.R1$R1)/(cid79_island.R1$R1_Input))/R1.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="Contig79",col="grey",border=NA, x0=cid79_rep.R1$position,x1=cid79_rep.R1$position+1, y1=log(((cid79_rep.R1$R1)/(cid79_rep.R1$R1_Input))/R1.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="Y_Contig26",col=1,border=NA, x0=cid26.R1$position,x1=cid26.R1$position+1, y1=log(((cid26.R1$R1)/(cid26.R1$R1_Input))/R1.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="3R_5",col=1,border=NA,  x0=cid3_island.R1$position,x1=cid3_island.R1$position+1, y1=log(((cid3_island.R1$R1)/(cid3_island.R1$R1_Input))/R1.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="3R_5",col="grey",border=NA,  x0=cid3_rep.R1$position,x1=cid3_rep.R1$position+1, y1=log(((cid3_rep.R1$R1)/(cid3_rep.R1$R1_Input))/R1.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="tig00057289",col=1,border=NA,  x0=cid57289_island.R1$position,x1=cid57289_island.R1$position+1, y1=log(((cid57289_island.R1$R1)/(cid57289_island.R1$R1_Input))/R1.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="tig00057289",col="grey",border=NA,  x0=cid57289_rep.R1$position,x1=cid57289_rep.R1$position+1, y1=log(((cid57289_rep.R1$R1)/(cid57289_rep.R1$R1_Input))/R1.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpPlotRegions(kp, data=R1_chiptigs,data.panel = 2, r0=0.2,r1=0.3,avoid.overlapping=FALSE,border=NA)
kpPlotRegions(kp, data=R1_peaks,data.panel = 2, r0=0.4,r1=0.5,avoid.overlapping=FALSE,col="orange",border=NA)
dev.copy(pdf,"Figure_2_grey_R1_no_AATAT_small.pdf", width=12, height=10)
dev.off()

kpBars(kp, chr="Contig119",col=1,border=NA, x0=cid119.R1$position,x1=cid119.R1$position+1, y1=log(((cid119.R1$R1)/(cid119.R1$R1_Input))/R1.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="Contig79",col=1,border=NA, x0=cid79.R1$position,x1=cid79.R1$position+1, y1=log(((cid79.R1$R1)/(cid79.R1$R1_Input))/R1.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="Y_Contig26",col=1,border=NA, x0=cid26.R1$position,x1=cid26.R1$position+1, y1=log(((cid26.R1$R1)/(cid26.R1$R1_Input))/R1.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="3R_5",col=1,border=NA,  x0=cid3.R1$position,x1=cid3.R1$position+1, y1=log(((cid3.R1$R1)/(cid3.R1$R1_Input))/R1.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="tig00057289",col=1,border=NA,  x0=cid57289.R1$position,x1=cid57289.R1$position+1, y1=log(((cid57289.R1$R1)/(cid57289.R1$R1_Input))/R1.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpPlotRegions(kp, data=R1_chiptigs,data.panel = 2, r0=0.2,r1=0.3,avoid.overlapping=FALSE,border=NA)
kpPlotRegions(kp, data=R1_peaks,data.panel = 2, r0=0.4,r1=0.5,avoid.overlapping=FALSE,col="orange",border=NA)


legend("topright",c( "10bp",  "AAGAG", "AAGAT" , "AATAG", "AATAT",  "AATATAT", "Dodeca", "Gypsy",  "Other_DNA", "IGS",   "Jockey","Other_LTR","Other_Non-LTR","Other_simple","Other"),col=c("darkcyan" , "blue"    , "springgreen3",  "cyan",    "black"  ,  "black", "blueviolet" , "goldenrod1" ,  "red", "yellow" , "maroon","dodgerblue3","forestgreen","grey","grey" ),pch=15)


#R2 (In Fig 2)
cid79_island.R2=subset(cid79.R2,cid79.R2$position>5655&cid79.R2$position<48432)
cid79_rep.R2=subset(cid79.R2,cid79.R2$position<=5655 | cid79.R2$position>=48432)
cid119_island.R2=subset(cid119.R2,cid119.R2$position>29252&cid119.R2$position<69146)
cid119_rep.R2=subset(cid119.R2,cid119.R2$position<=29252 | cid119.R2$position>=69146)
cid26_island.R2=subset(cid26.R2,cid26.R2$position>5655&cid26.R2$position<48432)
cid26_rep.R2=subset(cid26.R2,cid26.R2$position<=5655 | cid26.R2$position>=48432)
cid3_island.R2=subset(cid3.R2,cid3.R2$position>19152&cid3.R2$position<69541)
cid3_rep.R2=subset(cid3.R2,cid3.R2$position<=19152 | cid3.R2$position>=69541)
cid57289_island.R2=subset(cid57289.R2,cid57289.R2$position>5561&cid57289.R2$position<7864)
cid57289_rep.R2=subset(cid57289.R2,cid57289.R2$position<=5561 | cid57289.R2$position>=7864)


kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands,plot.type=2,ideogram.plotter = NULL,cex=.7,chromosomes=c("Contig79","Contig119","Y_Contig26","3R_5","tig00057289"))
kpAddCytobandsAsLine(kp,color.table=v,border=NA,lwd=5)
kpAddBaseNumbers(kp,tick.dist=10000)
kpAxis(kp, r0=0.5, r1=1, tick.pos = c(0,2.5, 5),ymax=5,cex=.5)
kpAxis(kp, r0=0.5, r1=0,  tick.pos = c(0,-2.5, -5),ymax=-5,cex=.5)
kpBars(kp, chr="Contig119",col=1,border=NA, x0=cid119_island.R2$position,x1=cid119_island.R2$position+1, y1=log(((cid119_island.R2$R2)/(cid119_island.R2$R2_Input))/R2.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="Contig119",col="grey",border=NA, x0=cid119_rep.R2$position,x1=cid119_rep.R2$position+1, y1=log(((cid119_rep.R2$R2)/(cid119_rep.R2$R2_Input))/R2.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="Contig79",col=1,border=NA, x0=cid79_island.R2$position,x1=cid79_island.R2$position+1, y1=log(((cid79_island.R2$R2)/(cid79_island.R2$R2_Input))/R2.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="Contig79",col="grey",border=NA, x0=cid79_rep.R2$position,x1=cid79_rep.R2$position+1, y1=log(((cid79_rep.R2$R2)/(cid79_rep.R2$R2_Input))/R2.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="Y_Contig26",col=1,border=NA, x0=cid26.R2$position,x1=cid26.R2$position+1, y1=log(((cid26.R2$R2)/(cid26.R2$R2_Input))/R2.norm,2)/5,data.panel = 1, r0=0.5,r1=1)

kpBars(kp, chr="3R_5",col=1,border=NA,  x0=cid3_island.R2$position,x1=cid3_island.R2$position+1, y1=log(((cid3_island.R2$R2)/(cid3_island.R2$R2_Input))/R2.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="3R_5",col="grey",border=NA,  x0=cid3_rep.R2$position,x1=cid3_rep.R2$position+1, y1=log(((cid3_rep.R2$R2)/(cid3_rep.R2$R2_Input))/R2.norm,2)/5,data.panel = 1, r0=0.5,r1=1)

kpBars(kp, chr="tig00057289",col=1,border=NA,  x0=cid57289_island.R2$position,x1=cid57289_island.R2$position+1, y1=log(((cid57289_island.R2$R2)/(cid57289_island.R2$R2_Input))/R2.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="tig00057289",col="grey",border=NA,  x0=cid57289_rep.R2$position,x1=cid57289_rep.R2$position+1, y1=log(((cid57289_rep.R2$R2)/(cid57289_rep.R2$R2_Input))/R2.norm,2)/5,data.panel = 1, r0=0.5,r1=1)

kpPlotRegions(kp, data=R2_chiptigs,data.panel = 2, r0=0.2,r1=0.3,avoid.overlapping=FALSE,border=NA)
kpPlotRegions(kp, data=R2_peaks,data.panel = 2, r0=0.4,r1=0.5,avoid.overlapping=FALSE,col="orange",border=NA)
#kpPlotRegions(kp, data=primers,data.panel = 2, r0=0.6,r1=0.7,avoid.overlapping=FALSE,col="green")
dev.copy(pdf,"myplot.pdf", width=12, height=10)
dev.off()

#R3 (in Fig S4)
cid79_island.R3=subset(cid79.R3,cid79.R3$position>5655&cid79.R3$position<48432)
cid79_rep.R3=subset(cid79.R3,cid79.R3$position<=5655 | cid79.R3$position>=48432)
cid119_island.R3=subset(cid119.R3,cid119.R3$position>29252&cid119.R3$position<69146)
cid119_rep.R3=subset(cid119.R3,cid119.R3$position<=29252 | cid119.R3$position>=69146)
cid26_island.R3=subset(cid26.R3,cid26.R3$position>5655&cid26.R3$position<48432)
cid26_rep.R3=subset(cid26.R3,cid26.R3$position<=5655 | cid26.R3$position>=48432)
cid3_island.R3=subset(cid3.R3,cid3.R3$position>19152&cid3.R3$position<69541)
cid3_rep.R3=subset(cid3.R3,cid3.R3$position<=19152 | cid3.R3$position>=69541)
cid57289_island.R3=subset(cid57289.R3,cid57289.R3$position>5561&cid57289.R3$position<7864)
cid57289_rep.R3=subset(cid57289.R3,cid57289.R3$position<=5561 | cid57289.R3$position>=7864)

kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands,plot.type=2,ideogram.plotter = NULL,cex=.7,chromosomes=c("Contig79","Contig119","Y_Contig26","3R_5","tig00057289"))
kpAddCytobandsAsLine(kp,color.table=v,lwd=5)
kpAddBaseNumbers(kp,tick.dist=10000)
kpAxis(kp, r0=0.5, r1=1, tick.pos = c(0,2.5, 5),ymax=5,cex=.5)
kpAxis(kp, r0=0.5, r1=0,  tick.pos = c(0,-2.5, -5),ymax=-5,cex=.5)
kpBars(kp, chr="Contig119",col=1,border=NA, x0=cid119_island.R3$position,x1=cid119_island.R3$position+1, y1=log(((cid119_island.R3$R3)/(cid119_island.R3$R3_Input))/R3.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="Contig119",col="grey",border=NA, x0=cid119_rep.R3$position,x1=cid119_rep.R3$position+1, y1=log(((cid119_rep.R3$R3)/(cid119_rep.R3$R3_Input))/R3.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="Contig79",col=1,border=NA, x0=cid79_island.R3$position,x1=cid79_island.R3$position+1, y1=log(((cid79_island.R3$R3)/(cid79_island.R3$R3_Input))/R3.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="Contig79",col="grey",border=NA, x0=cid79_rep.R3$position,x1=cid79_rep.R3$position+1, y1=log(((cid79_rep.R3$R3)/(cid79_rep.R3$R3_Input))/R3.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="Y_Contig26",col=1,border=NA, x0=cid26.R3$position,x1=cid26.R3$position+1, y1=log(((cid26.R3$R3)/(cid26.R3$R3_Input))/R3.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="3R_5",col=1,border=NA,  x0=cid3_island.R3$position,x1=cid3_island.R3$position+1, y1=log(((cid3_island.R3$R3)/(cid3_island.R3$R3_Input))/R3.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="3R_5",col="grey",border=NA,  x0=cid3_rep.R3$position,x1=cid3_rep.R3$position+1, y1=log(((cid3_rep.R3$R3)/(cid3_rep.R3$R3_Input))/R3.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="tig00057289",col=1,border=NA,  x0=cid57289_island.R3$position,x1=cid57289_island.R3$position+1, y1=log(((cid57289_island.R3$R3)/(cid57289_island.R3$R3_Input))/R3.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="tig00057289",col="grey",border=NA,  x0=cid57289_rep.R3$position,x1=cid57289_rep.R3$position+1, y1=log(((cid57289_rep.R3$R3)/(cid57289_rep.R3$R3_Input))/R3.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpPlotRegions(kp, data=R3_chiptigs,data.panel = 2, r0=0.2,r1=0.3,avoid.overlapping=FALSE,border=NA)
kpPlotRegions(kp, data=R3_peaks,data.panel = 2, r0=0.4,r1=0.5,avoid.overlapping=FALSE,col="orange",border=NA)
dev.copy(pdf,"Figure_2_grey_R3_no_AATAT_small.pdf", width=12, height=10)
dev.off()


kpBars(kp, chr="Contig119", x0=cid119.R3$position,x1=cid119.R3$position+1, y1=log(((cid119.R3$R3)/(cid119.R3$R3_Input))/R3.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="Contig79", x0=cid79.R3$position,x1=cid79.R3$position+1, y1=log(((cid79.R3$R3)/(cid79.R3$R3_Input))/R3.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="Y_Contig26", x0=cid26.R3$position,x1=cid26.R3$position+1, y1=log(((cid26.R3$R3)/(cid26.R3$R3_Input))/R3.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="3R_5",  x0=cid3.R3$position,x1=cid3.R3$position+1, y1=log(((cid3.R3$R3)/(cid3.R3$R3_Input))/R3.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="tig00057289",  x0=cid57289.R3$position,x1=cid57289.R3$position+1, y1=log(((cid57289.R3$R3)/(cid57289.R3$R3_Input))/R3.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpPlotRegions(kp, data=R3_chiptigs,data.panel = 2, r0=0.2,r1=0.3,avoid.overlapping=FALSE)
kpPlotRegions(kp, data=R3_peaks,data.panel = 2, r0=0.4,r1=0.5,avoid.overlapping=FALSE,col="orange")

#R4 (in Fig S4)
cid79_island.R4=subset(cid79.R4,cid79.R4$position>5655&cid79.R4$position<48432)
cid79_rep.R4=subset(cid79.R4,cid79.R4$position<=5655 | cid79.R4$position>=48432)
cid119_island.R4=subset(cid119.R4,cid119.R4$position>29252&cid119.R4$position<69146)
cid119_rep.R4=subset(cid119.R4,cid119.R4$position<=29252 | cid119.R4$position>=69146)
cid26_island.R4=subset(cid26.R4,cid26.R4$position>5655&cid26.R4$position<48432)
cid26_rep.R4=subset(cid26.R4,cid26.R4$position<=5655 | cid26.R4$position>=48432)
cid3_island.R4=subset(cid3.R4,cid3.R4$position>19152&cid3.R4$position<69541)
cid3_rep.R4=subset(cid3.R4,cid3.R4$position<=19152 | cid3.R4$position>=69541)
cid57289_island.R4=subset(cid57289.R4,cid57289.R4$position>5561&cid57289.R4$position<7864)
cid57289_rep.R4=subset(cid57289.R4,cid57289.R4$position<=5561 | cid57289.R4$position>=7864)

kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands,plot.type=2,ideogram.plotter = NULL,cex=.7,chromosomes=c("Contig79","Contig119","Y_Contig26","3R_5","tig00057289"))
kpAddCytobandsAsLine(kp,color.table=v,lwd=5)
kpAddBaseNumbers(kp,tick.dist=10000)
kpAxis(kp, r0=0.5, r1=1, tick.pos = c(0,2.5, 5),ymax=5,cex=.5)
kpAxis(kp, r0=0.5, r1=0,  tick.pos = c(0,-2.5, -5),ymax=-5,cex=.5)
kpBars(kp, chr="Contig119",col=1,border=NA, x0=cid119_island.R4$position,x1=cid119_island.R4$position+1, y1=log(((cid119_island.R4$R4)/(cid119_island.R4$R4_Input))/R4.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="Contig119",col="grey",border=NA, x0=cid119_rep.R4$position,x1=cid119_rep.R4$position+1, y1=log(((cid119_rep.R4$R4)/(cid119_rep.R4$R4_Input))/R4.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="Contig79",col=1,border=NA, x0=cid79_island.R4$position,x1=cid79_island.R4$position+1, y1=log(((cid79_island.R4$R4)/(cid79_island.R4$R4_Input))/R4.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="Contig79",col="grey",border=NA, x0=cid79_rep.R4$position,x1=cid79_rep.R4$position+1, y1=log(((cid79_rep.R4$R4)/(cid79_rep.R4$R4_Input))/R4.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="Y_Contig26",col=1,border=NA, x0=cid26.R4$position,x1=cid26.R4$position+1, y1=log(((cid26.R4$R4)/(cid26.R4$R4_Input))/R4.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="3R_5",col=1,border=NA,  x0=cid3_island.R4$position,x1=cid3_island.R4$position+1, y1=log(((cid3_island.R4$R4)/(cid3_island.R4$R4_Input))/R4.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="3R_5",col="grey",border=NA,  x0=cid3_rep.R4$position,x1=cid3_rep.R4$position+1, y1=log(((cid3_rep.R4$R4)/(cid3_rep.R4$R4_Input))/R4.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="tig00057289",col=1,border=NA,  x0=cid57289_island.R4$position,x1=cid57289_island.R4$position+1, y1=log(((cid57289_island.R4$R4)/(cid57289_island.R4$R4_Input))/R4.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="tig00057289",col="grey",border=NA,  x0=cid57289_rep.R4$position,x1=cid57289_rep.R4$position+1, y1=log(((cid57289_rep.R4$R4)/(cid57289_rep.R4$R4_Input))/R4.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpPlotRegions(kp, data=R4_chiptigs,data.panel = 2, r0=0.2,r1=0.3,avoid.overlapping=FALSE,border=NA)
kpPlotRegions(kp, data=R4_peaks,data.panel = 2, r0=0.4,r1=0.5,avoid.overlapping=FALSE,col="orange",border=NA)
dev.copy(pdf,"Figure_2_grey_R4_no_AATAT_small.pdf", width=12, height=10)
dev.off()




kpBars(kp, chr="Contig119",border=NA, x0=cid119.R4$position,x1=cid119.R4$position+1, y1=log(((cid119.R4$R4)/(cid119.R4$R4_Input))/R4.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="Contig79",border=NA, x0=cid79.R4$position,x1=cid79.R4$position+1, y1=log(((cid79.R4$R4)/(cid79.R4$R4_Input))/R4.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="Y_Contig26",border=NA, x0=cid26.R4$position,x1=cid26.R4$position+1, y1=log(((cid26.R4$R4)/(cid26.R4$R4_Input))/R4.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="3R_5",border=NA,  x0=cid3.R4$position,x1=cid3.R4$position+1, y1=log(((cid3.R4$R4)/(cid3.R4$R4_Input))/R4.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpBars(kp, chr="tig00057289",border=NA,  x0=cid57289.R4$position,x1=cid57289.R4$position+1, y1=log(((cid57289.R4$R4)/(cid57289.R4$R4_Input))/R4.norm,2)/5,data.panel = 1, r0=0.5,r1=1)
kpPlotRegions(kp, data=R4_chiptigs,data.panel = 2, r0=0.2,r1=0.3,avoid.overlapping=FALSE,border=NA)
kpPlotRegions(kp, data=R4_peaks,data.panel = 2, r0=0.4,r1=0.5,avoid.overlapping=FALSE,col="orange",border=NA)
dev.copy(pdf,"Figure_2_grey_R4_no_AATAT_small.pdf", width=12, height=10)
dev.off()

#For Figure 5 oligopaint plots
#cen3
kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands,plot.type=2,ideogram.plotter = NULL,cex=.7,chromosomes="3R_5")
kpAddCytobandsAsLine(kp,color.table=v,lwd=5)
kpAddBaseNumbers(kp,tick.dist=10000)
kpPlotRegions(kp, data=cen3oligo,data.panel = 2, r0=0.2,r1=0.3,avoid.overlapping=FALSE,col="magenta")

#cenX
kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands,plot.type=2,ideogram.plotter = NULL,cex=.7,chromosomes="Contig79")
kpAddCytobandsAsLine(kp,color.table=v,lwd=5)
kpAddBaseNumbers(kp,tick.dist=10000)
kpPlotRegions(kp, data=cenXoligo,data.panel = 2, r0=0.2,r1=0.3,avoid.overlapping=FALSE,col="magenta")

#cenY
kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands,plot.type=2,ideogram.plotter = NULL,cex=.7,chromosomes="Y_Contig26")
kpAddCytobandsAsLine(kp,color.table=v,lwd=5)
kpAddBaseNumbers(kp,tick.dist=10000)
kpPlotRegions(kp, data=cenYoligo,data.panel = 2, r0=0.2,r1=0.3,avoid.overlapping=FALSE,col="magenta")

#cen4
kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands,plot.type=2,ideogram.plotter = NULL,cex=.7,chromosomes="Contig119")
kpAddCytobandsAsLine(kp,color.table=v,lwd=5)
kpAddBaseNumbers(kp,tick.dist=10000)
kpPlotRegions(kp, data=cen4oligo,data.panel = 2, r0=0.2,r1=0.3,avoid.overlapping=FALSE,col="magenta")
