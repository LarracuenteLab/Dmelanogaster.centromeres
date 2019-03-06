
#############
#	R script to make Figure 6 of Chang et al. 2019 PLoS Biology
# The association between G2/Jockey-3 and centromeres is conserved in D. simulans. 
# A) Plot of the normalized CENP-A enrichment over input across the D. simulans 
# G2/Jockey-3 consensus sequence using CENP-A ChIP-seq data from D. simulans 
# ML82-19a cells (Talbert et al. 2018) showing that G2/Jockey-3 is enriched in 
# CENP-A in D. simulans. 15m and 5m indicate duration of MNase digestion and IP 
# and IP2 are technical replicates. Note that the first 487bp of D. simulans 
# G2/Jockey-3 consensus sequence, which are homologous to the 500bp satellite, 
# are not included in this figure; the 500bp satellite was previously reported 
# as enriched in CENP-A in D. simulans [16]. B) Plot of the normalized CENP-A 
# enrichment over input across the D. melanogaster G2/Jockey-3 consensus sequence 
# using our CENP-A ChIP-seq replicates (R1-R4) and ChIP-seq from CENP-A-GFP transgenic
# flies from Talbert et al. 2018. 
#
# By Ching-Ho Chang
# cchang45@ur.rochester.edu
#############

#Fig 6A D. simulans
pdf("Jockey3_dsim_enrichement.pdf",8,4)
data<-read.table("Fig6A.All_repeat_Jockey_sim.out",header=T)
#Transposon pos CidM_15m_IP CidM_15m_IN CidM_5m_IP CidM_5m_IN CidM_5m_IP2 CidM_5m_IN2
#normalized by total reads...
data$CidM_15m_IP=data$CidM_15m_IP/39368106
data$CidM_5m_IP=data$CidM_5m_IP/31610618
data$CidM_5m_IP2=data$CidM_5m_IP2/17139506
data$CidM_15m_IN=(data$CidM_15m_IN+1)/62603992
data$CidM_5m_IN=(data$CidM_5m_IN+1)/56951406
data$CidM_5m_IN2=(data$CidM_5m_IN2+1)/31712194
plot(data$pos,log(data$CidM_15m_IP/(data$CidM_15m_IN),2),col=1,type='l',xlab='Jockey3 position',ylab=expression('log'[2]*'(ChIP/Input)'),ylim=c(0,10))
lines(data$pos,log(data$CidM_5m_IP/(data$CidM_5m_IN),2),col="orange")
lines(data$pos,log(data$CidM_5m_IP2/(data$CidM_5m_IN2),2),col="blue")
legend("topright", legend=c("15m_IP","5m_IP","5m_IP2"), col=c(1,"orange","blue"), lwd=1,lty=c(1,1,1),cex=0.8,bty='n')
#abline(v=487,col = "lightgray", lty = 3)
dev.off()
#normalized by mapped reads...
#data$CidM_15m_IP=data$CidM_15m_IP/27780864
#data$CidM_5m_IP=data$CidM_5m_IP/16424658
#data$CidM_5m_IP2=data$CidM_5m_IP2/10841216
#data$CidM_15m_IN=(data$CidM_15m_IN+1)/38744532
#data$CidM_5m_IN=(data$CidM_5m_IN+1)/35037652
#data$CidM_5m_IN2=(data$CidM_5m_IN2+1)/16062990


####### Fig 6B D. melanogaster
pdf("All_repeat_Jockey_DM.pdf",8,5)
data<-read.table("Fig6B.All_repeat_Jockey_DM.out",header=T)

data$ChIP_R1=data$ChIP_R1/49512826
data$ChIP_R2=data$ChIP_R2/393640001
data$ChIP_R3=data$ChIP_R3/63807709
data$ChIP_R5=data$ChIP_R5/78211153
data$ChIP_CidGFP=data$ChIP_CidGFP/34698260
data$Input_R1=(data$Input_R1+1)/101043744
data$Input_R2=(data$Input_R2+1)/495457706
data$Input_R3=(data$Input_R3+1)/115587082
data$Input_R5=(data$Input_R5+1)/43451898
data$Input_CidGFP=data$Input_CidGFP/30648608
plot(data$pos,log(data$ChIP_R1/(data$Input_R1),2),col=1,type='l',xlab='G2 position',ylab=expression('log'[2]*'(ChIP/Input)'),ylim=c(0,14))
lines(data$pos,log(data$ChIP_R2/(data$Input_R2),2),col="orange",lwd=2)
lines(data$pos,log(data$ChIP_R3/(data$Input_R3),2),col="blue",lwd=2)
lines(data$pos,log(data$ChIP_R5/(data$Input_R5),2),col="forestgreen",lwd=2)
lines(data$pos,log(data$ChIP_CidGFP/(data$Input_CidGFP),2),col="grey", lty=2,lwd=2)
legend("topright", legend=c("R1","R2","R3","R5","CidGFP"), col=c(1,"orange","blue","forestgreen","grey"),lty=c(1,1,1,1,2),lwd=2,cex=0.8,bty='n')
dev.off()