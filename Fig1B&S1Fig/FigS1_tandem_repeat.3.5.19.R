
#############
#	R script to make Figure S1 of Chang et al. 2019 PLoS Biology
# Figure S1. Enrichment of simple tandem repeats in CENP-A ChIP-seq across four 
# replicates. Plot of normalized CENP-A/Input for simple tandem repeats for each 
# ChIP-seq replicate, sorted by median (red lines). Shown are only the simple tandem 
# repeats with median CENP-A/Input > 1 in all four CENP-A ChIP replicates 
# (see details in Table S1). The simple tandem repeats with less than 10 counts of 
# input reads in any one replicate are not shown. 
#   
# By Ching-Ho Chang
# cchang45@ur.rochester.edu
#
#############
library(ggplot2)

total<-read.table("kseek_all_enrichment_FigS1",header=T)

plot.median <- function(x) {
  m <- median(x)
  c(y = m, ymin = m, ymax = m)
}
plot.average <- function(x) {
  m <- mean(x)
  c(y = m, ymin = m, ymax = m)
}

total$Repeat <- factor(total$Repeat, levels=c("ACCAGTACGGGACCGGG","AGC","ACCGAGTACGGG","AACCGAGTACGG","ACCAGTACGGG","AAGGCAATGC","ACGAG","AAGGC","AATGG","AGG","ACAGG","AATAG","AGAGG","AAGGAC","AG","AATAGACGAC","AAGACAAGAGAC","AC","AAGACAAGACAC","AACACAACACC","ACACC"))
p <- ggplot(total, aes(x=Repeat, y=Enrichment)) +
labs(list(title = "", x = "Repeat", y = "ChIP/Input")) +
theme(axis.title.x = element_text(face="bold"), axis.text.x = element_text(size = 8)) +
theme(axis.title.y = element_text(face="bold"), axis.text.y = element_text(size = 10))+
theme(panel.background = element_rect(fill = "white",  colour = "grey", linetype = "solid"))+
coord_cartesian(ylim=c(0,10))+ylab("ChIP/Input")
p1 <- p + stat_summary(fun.data="plot.median", geom="errorbar", colour="red", width=0.5, size=1) +
#stat_summary(fun.data="plot.average", geom="errorbar", colour="blue", width=0.5, size=1) +
geom_hline(yintercept=1,linetype="dashed", , size=1) + 
geom_dotplot(binaxis='y', stackdir='center', method="histodot", binwidth=0.3)+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p1)


