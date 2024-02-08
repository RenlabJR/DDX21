####FigureS3F
type="m6A"
type1=${type}_fwd
type2=${type}_rev
computeMatrix reference-point -p 20 --scoreFileName $Rloop_fwd -R ${type}.fwd.bed --referencePoint center -b 2000 -a 2000 --binSize 10 --skipZeros -o ${type1}.point.plot_matrix_mRNA.gz --outFileNameMatrix ${type1}.point.loadR_matrix.tab --averageTypeBins mean --samplesLabel shCTRL
computeMatrix reference-point -p 20 --scoreFileName $Rloop_rev -R ${type}.rev.bed --referencePoint center -b 2000 -a 2000 --binSize 10 --skipZeros -o ${type2}.point.plot_matrix_mRNA.gz --outFileNameMatrix ${type2}.point.loadR_matrix.tab --averageTypeBins mean --samplesLabel shCTRL
computeMatrixOperations rbind -m ${type1}.point.plot_matrix_mRNA.gz ${type2}.point.plot_matrix_mRNA.gz -o ${type}.merged.point.plot_matrix_mRNA.gz
plotHeatmap -m ${type}.merged.point.plot_matrix_mRNA.gz -out median.hot_spectral.${type}.point.merged_heatmap.pdf --outFileSortedRegions sorted.${type}.bed --sortRegions descend --sortUsing mean --heatmapWidth 1.3 --heatmapHeight 5 --colorList '#0A0000, #4C0000, #900000, #FF1700, #FF5B00, #FF9D00, #FFE100, #FFFF36, #FFFF9C, #A9DCA4, #66C2A5, #3286BC' --averageTypeSummaryPlot median --zMax 80


####FigureS3H
multiBamSummary BED-file -p 12 --BED m6A_Rloop.gene.bed -o m6A_Rloop.gene.DDX21_METTL3.npz --outRawCounts m6A_Rloop.gene.DDX21_METTL3.tab --bamfiles $DDX21 $METTL3
multiBamSummary BED-file -p 12 --BED no_m6A_Rloop.gene.bed -o no_m6A_Rloop.gene.DDX21_METTL3.npz --outRawCounts no_m6A_Rloop.gene.DDX21_METTL3.tab --bamfiles $DDX21 $METTL3

sed '1d' m6A_Rloop.gene.DDX21_METTL3.tab | awk -v OFS="\t" '{n_len=($3-$2); print $1, $2, $3, $4, $5, ($4*1e9)/(n_len*40344146), ($5*1e9)/(n_len*"'$DDX21'")}' > m6A_Rloop.gene.DDX21_METTL3.RPKM
sed '1d' no_m6A_Rloop.gene.DDX21_METTL3.tab | awk -v OFS="\t" '{n_len=($3-$2); print $1, $2, $3, $4, $5, ($4*1e9)/(n_len*40344146), ($5*1e9)/(n_len*"'$METTL3'")}' > no_m6A_Rloop.gene.DDX21_METTL3.RPKM

library(ggplot2)
library(RColorBrewer)
m6A_Rloop <- read.table("m6A_Rloop.gene.DDX21_METTL3.RPKM",header=F)
no_m6A_Rloop <- read.table("no_m6A_Rloop.gene.DDX21_METTL3.RPKM",header=F)
wilcox.test(m6A_Rloop$V6, no_m6A_Rloop$V6, alternative="greater")
wilcox.test(m6A_Rloop$V7, no_m6A_Rloop$V7, alternative="greater")

library(ggsci)
cols = pal_npg("nrc", alpha = 0.7)(9)
#DDX21
data <- data.frame(RPKM=c(m6A_Rloop$V6, no_m6A_Rloop$V6), type=c(rep("m6A+ Rloop+ gene", nrow(m6A_Rloop)), rep("m6A- Rloop+ gene", nrow(no_m6A_Rloop))))
data$type <- factor(data$type, level=c("m6A- Rloop+ gene", "m6A+ Rloop+ gene"))
pdf("DDX21 binding in m6A R_loop genes boxplot.pdf",width=2.1,height=2.3)
ggplot(data, aes(type, log2(RPKM)))+ 
  geom_boxplot(fill=cols[1:2], width = 0.3, outlier.shape = NA)+
  labs(x="",y="log2(RPKM)",title="") + 
  scale_fill_manual(values=c(brewer.pal(5,"Set2")[c(1,3)]))+
  theme_classic() + theme(
  legend.title=element_blank(),
  axis.title=element_text(size=12, color="black"),
  axis.text.x=element_text(size=12, color="black"),
  axis.text.y=element_text(size=12, color="black"),
  legend.text=element_text(size=12, color="black"),
  axis.ticks.x=element_blank())
dev.off()
#METTL3
data <- data.frame(RPKM=c(m6A_Rloop$V7, no_m6A_Rloop$V7), type=c(rep("m6A+ Rloop+ gene", nrow(m6A_Rloop)), rep("m6A- Rloop+ gene", nrow(no_m6A_Rloop))))
data$type <- factor(data$type, level=c("m6A- Rloop+ gene", "m6A+ Rloop+ gene"))
pdf("METTL3 binding in m6A R_loop genes boxplot.pdf",width=2.1,height=2.3)
ggplot(data, aes(type, log2(RPKM)))+ 
  geom_boxplot(fill=cols[1:2], width = 0.3, outlier.shape = NA)+
  labs(x="",y="log2(RPKM)",title="") + 
  scale_fill_manual(values=c(brewer.pal(5,"Set2")[c(1,3)]))+
  theme_classic() + theme(
  legend.title=element_blank(),
  axis.title=element_text(size=12, color="black"),
  axis.text.x=element_text(size=12, color="black"),
  axis.text.y=element_text(size=12, color="black"),
  legend.text=element_text(size=12, color="black"),
  axis.ticks.x=element_blank())
dev.off()


####FigureS3I
awk -v OFS="\t" '{print $0, "Rloop_"NR}' $Rloop_macs2 | intersectBed -a - -b m6A_Rloop.gene.bed -wa -wb -f 0.51 | awk '$6==$14' | grep mRNA | cut -f7,8,12 | awk '!a[$0]++' | awk -v OFS="\t" '{print $3,$2,$1}' > m6A_Rloop.gene.Rloop.peak.bed
awk -v OFS="\t" '{print $0, "Rloop_"NR}' $Rloop_macs2 | intersectBed -a - -b no_m6A_Rloop.gene.bed -wa -wb -f 0.51 | awk '$6==$14' | grep mRNA | cut -f7,8,12 | awk '!a[$0]++' | awk -v OFS="\t" '{print $3,$2,$1}' > no_m6A_Rloop.gene.Rloop.peak.bed

options(stringsAsFactors=F)
library(plyr)
library(RColorBrewer)	
cols <- brewer.pal(9,"Set1")[c(2,4)]
m6A_Rloop <- read.table("m6A_Rloop.gene.Rloop.peak.bed",header=F)
no_m6A_Rloop <- read.table("no_m6A_Rloop.gene.Rloop.peak.bed",header=F)
m6A_Rloop <- ddply(m6A_Rloop, .(V1), sum_Rloop=sum(V3), summarise)
no_m6A_Rloop <- ddply(no_m6A_Rloop, .(V1), sum_Rloop=sum(V3), summarise)

library(ggplot2)
data <- rbind(m6A_Rloop, no_m6A_Rloop)
data$type <- c(rep("m6A+ R-loop genes",nrow(m6A_Rloop)), rep("m6A- R-loop genes",nrow(no_m6A_Rloop)))
data$type <- factor(data$type, level=c("m6A- R-loop genes", "m6A+ R-loop genes"))

p1 <- wilcox.test(m6A_Rloop$sum_Rloop,no_m6A_Rloop$sum_Rloop,alternative="greater")$p.value
pdf("Rloop_m6A_genes_Rloop_level.pdf", width=4.1,height=4.5)
ggplot(data, aes(type,log2(sum_Rloop)))+
  geom_boxplot(aes(fill = type), outlier.shape=NA,width=0.42)+
  labs(x="",y="log2(R-loop level)",title="") +
  scale_fill_manual(values=cols)+
  theme_classic() + theme(
  legend.title=element_blank(),
  axis.title=element_text(size=16),
  axis.text.x=element_text(size=14,angle=30,hjust=1,vjust=1,color="Black"),
  axis.text.y=element_text(size=14,color="Black"),
  legend.position="none",
  axis.ticks.x=element_blank()) +
  annotate("text", x = 1, y = -1, label = paste("p=",p1))  #+
dev.off()








