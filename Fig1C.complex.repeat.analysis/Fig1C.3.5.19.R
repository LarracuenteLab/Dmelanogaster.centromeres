
#############
#	R script to make Figure 1C of Chang et al. 2019 PLoS Biology
# Plot of the normalized CENP-A/Input reads on a log scale for each replicate, 
# sorted by median (red lines) for complex repeat families. Shown are only the 
# complex repeats in the top 20% across all four CENP-A ChIP replicates.  
#   
# By Ching-Ho Chang
# cchang45@ur.rochester.edu
# 
#############
library(ggplot2)


total<-read.table("repeat.plots_complex_repeat_data",header=T)
plot.median <- function(x) {
  m <- median(x)
  c(y = m, ymin = m, ymax = m)
}
plot.average <- function(x) {
  m <- mean(x)
  c(y = m, ymin = m, ymax = m)
}

#p <- ggplot(total, aes(x=Repeat, y=Enrichment)) +
#labs(list(title = "", x = "Repeat", y = "ChIP/Input")) +

total$Repeat <- factor(total$Repeat, levels=c("Jockey-3_Dmel","Jockey-1_DSi","NTS_3cen_DM","G6_DM","DMRT1B","Gypsy8_LTR","DOC2_DM","DM1731_I","TART-A","G_DM","BLASTOPIA_LTR","Chimpo_I","Bica_LTR","R2_DM","BEL-6_DYa-I"))
p <- ggplot(total, aes(x=Repeat, y=log(Enrichment,2))) +
coord_cartesian(ylim=c(0,10))+ylab("log2 ChIP/Input")+
theme(axis.title.x = element_text(face="bold"), axis.text.x = element_text(size = 10)) +
theme(axis.title.y = element_text(face="bold"), axis.text.y = element_text(size = 10))+
theme(panel.background = element_rect(fill = "white",  colour = "grey", linetype = "solid"))
ylab("ChIP/Input")+coord_cartesian(ylim=c(0,10))

p1 <- p + stat_summary(fun.data="plot.median", geom="errorbar", colour="red", width=0.5, size=1) +
geom_dotplot(binaxis='y', stackdir='center', method="histodot", binwidth=0.3) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p1)
