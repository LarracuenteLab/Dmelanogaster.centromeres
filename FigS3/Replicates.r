#############
# R script to make Figure S3A of Chang et al. 2019 PLoS Biology
# The part of Figure S3A is modified from the R plot generated from this script by Illustrator.  
# By Ching-Ho Chang
# cchang45@ur.rochester.edu
#############

data<-read.table("idr_peak_comparision",na.strings = "NA",fill=T,head=T)

par(lwd=4,oma=c(3, 3, 7, 7),mar=c(1,1,1,1))
layout(matrix(1:16,4,4))
plot(data$R12_R1,data$R12_R2,pch=20,cex=1,ann=FALSE,axes=FALSE,frame.plot=TRUE,lwd=3,xlim=c(0,3000),ylim=c(0,3000),type="n")
text(1500,1500,"R1",cex=6)
plot(data$R12_R1,data$R12_R2,pch=20,cex=1,ann=FALSE,axes=FALSE,frame.plot=TRUE,lwd=3,xlim=c(0,3000),ylim=c(0,3000),type="n")
text(1500,1500,formatC(cor.test(data$R12_R1,data$R12_R2,method = c("spearman"))$estimate,digits=3, format="f"),cex=2)
plot(data$R12_R1,data$R12_R2,pch=20,cex=1,ann=FALSE,axes=FALSE,frame.plot=TRUE,lwd=3,xlim=c(0,3000),ylim=c(0,3000),type="n")
text(1500,1500,formatC(cor.test(data$R13_R1,data$R13_R3,method = c("spearman"))$estimate,digits=3, format="f"),cex=2)
plot(data$R12_R1,data$R12_R2,pch=20,cex=1,ann=FALSE,axes=FALSE,frame.plot=TRUE,lwd=3,xlim=c(0,3000),ylim=c(0,3000),type="n")
text(1500,1500,formatC(cor.test(data$R14_R1,data$R14_R4,method = c("spearman"))$estimate,digits=3, format="f"),cex=2)

plot(data$R12_R1,data$R12_R2,pch=20,cex=0.5,ann=FALSE,axes=FALSE,frame.plot=TRUE,lwd=1,xlim=c(0,3000),ylim=c(0,3000))
axis(3,labels=c("0","1","2","3","4","5"),at=c(0,1000,2000,3000,4000,3000),line=NA, cex.axis=1.5, lwd.ticks=5,tck=-0.02, mgp=c(1,0.5,0))
plot(data$R12_R2,data$R12_R1,pch=20,cex=0.5,ann=FALSE,axes=FALSE,frame.plot=TRUE,lwd=3,xlim=c(0,3000),ylim=c(0,3000),type="n")
text(1500,1500,"R2",cex=6)
plot(data$R12_R1,data$R12_R2,pch=20,cex=1,ann=FALSE,axes=FALSE,frame.plot=TRUE,lwd=3,xlim=c(0,3000),ylim=c(0,3000),type="n")
text(1500,1500,formatC(cor.test(data$R23_R2,data$R23_R3,method = c("spearman"))$estimate,digits=3, format="f"),cex=2)
plot(data$R12_R1,data$R12_R2,pch=20,cex=1,ann=FALSE,axes=FALSE,frame.plot=TRUE,lwd=3,xlim=c(0,3000),ylim=c(0,3000),type="n")
text(1500,1500,formatC(cor.test(data$R24_R2,data$R24_R4,method = c("spearman"))$estimate,digits=3, format="f"),cex=2)

plot(data$R13_R1,data$R13_R3,pch=20,cex=0.5,ann=FALSE,axes=FALSE,frame.plot=TRUE,lwd=1,xlim=c(0,3000),ylim=c(0,3000))
axis(3,labels=c("0","1","2","3","4","5"),at=c(0,1000,2000,3000,4000,3000),line=NA, cex.axis=1.5, lwd.ticks=5,tck=-0.02, mgp=c(1,0.5,0))
plot(data$R23_R2,data$R23_R3,pch=20,cex=0.5,ann=FALSE,axes=FALSE,frame.plot=TRUE,lwd=1,xlim=c(0,3000),ylim=c(0,3000))
plot(data$R12_R1,data$R12_R2,pch=20,cex=0.5,ann=FALSE,axes=FALSE,frame.plot=TRUE,lwd=1,xlim=c(0,3000),ylim=c(0,3000),type="n")
text(1500,1500,"R3",cex=6)
plot(data$R34_R3,data$R34_R4,pch=20,cex=0.5,ann=FALSE,axes=FALSE,frame.plot=TRUE,lwd=3,xlim=c(0,3000),ylim=c(0,3000),type="n")
text(1500,1500,formatC(cor.test(data$R34_R3,data$R34_R4,method = c("spearman"))$estimate,digits=3, format="f"),cex=2)

plot(data$R14_R1,data$R14_R4,pch=20,cex=0.5,ann=FALSE,axes=FALSE,frame.plot=TRUE,lwd=1,xlim=c(0,3000),ylim=c(0,3000))
axis(3,labels=c("0","1","2","3","4","5"),at=c(0,1000,2000,3000,4000,3000),line=NA, cex.axis=1.5, lwd.ticks=5,tck=-0.02, mgp=c(1,0.5,0))
axis(4,labels=c("0","1","2","3","4","5"),at=c(0,1000,2000,3000,4000,3000),line=NA, cex.axis=1.5, lwd.ticks=5,tck=-0.02, mgp=c(1,0.5,0))
plot(data$R24_R2,data$R24_R4,pch=20,cex=0.5,ann=FALSE,axes=FALSE,frame.plot=TRUE,lwd=1,xlim=c(0,3000),ylim=c(0,3000))
axis(4,labels=c("0","1","2","3","4","5"),at=c(0,1000,2000,3000,4000,3000),line=NA, cex.axis=1.5, lwd.ticks=5,tck=-0.02, mgp=c(1,0.5,0))
plot(data$R34_R3,data$R34_R4,pch=20,cex=0.5,ann=FALSE,axes=FALSE,frame.plot=TRUE,lwd=1,xlim=c(0,3000),ylim=c(0,3000))
axis(4,labels=c("0","1","2","3","4","5"),at=c(0,1000,2000,3000,4000,3000),line=NA, cex.axis=1.5, lwd.ticks=5,tck=-0.02, mgp=c(1,0.5,0))
plot(data$R12_R1,data$R12_R2,pch=20,cex=0.5,ann=FALSE,axes=FALSE,frame.plot=TRUE,lwd=3,xlim=c(0,3000),ylim=c(0,3000),type="n")
text(1500,1500,"R4",cex=6)

jpeg("Peakvalue_correlation.jpg",width = 4800, height = 4800)



dev.off()

 