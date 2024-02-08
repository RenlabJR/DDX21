####Figure2A
multiBamSummary BED-file -p 12 --BED $genebody -o RNA_signal.npz --outRawCounts RNA_signal.tab --bamfiles \
./shCTRL_DDX21/no_MT.bam \
./shCTRL_METTL3/no_MT.bam 

options(stringsAsFactors=F)
data <- read.table("RNA_signal.tab",header=F,comment.char="#")
colnames(data) <- c("chr","start","end","DDX21","METTL3")
data$DDX21 <- data$DDX21*1e6/DDX21_reads
data$METTL3 <- data$METTL3*1e6/METTL3_reads
cor.test(data$DDX21,data$METTL3)

pdf("DDX21_METTL3_RNA_signal_smoothScatter.pdf",width=4.2,height=4.5)
smoothScatter(log2(data$DDX21+1),log2(data$METTL3+1),xlab="DDX21",ylab="METTL3",main="PAR-CLIP") 
abline(lm(log2(data$METTL3+1)~log2(data$DDX21+1)),lty=2,col="red",lwd=2)
text(2,10,"r = 0.98")
text(2,8,"p-value < 2.2e-16")
dev.off()


####Figure2B
awk -v OFS="\t" '$6=="+"' $anno_DDX21 > DDX21.fwd.bed
awk -v OFS="\t" '$6=="-"' $anno_DDX21 > DDX21.rev.bed
awk -v OFS="\t" '$6=="+"' $anno_METTL3 > METTL3.fwd.bed
awk -v OFS="\t" '$6=="-"' $anno_METTL3 > METTL3.rev.bed
type="DDX21"
type1=${type}_fwd
type2=${type}_rev
computeMatrix reference-point -p 20 --scoreFileName $Rloop_fwd -R ${type}.fwd.bed --referencePoint center -b 2000 -a 2000 --binSize 10 --skipZeros -o ${type1}.point.plot_matrix_mRNA.gz --outFileNameMatrix ${type1}.point.loadR_matrix.tab --averageTypeBins mean
computeMatrix reference-point -p 20 --scoreFileName $Rloop_rev -R ${type}.rev.bed --referencePoint center -b 2000 -a 2000 --binSize 10 --skipZeros -o ${type2}.point.plot_matrix_mRNA.gz --outFileNameMatrix ${type2}.point.loadR_matrix.tab --averageTypeBins mean
computeMatrixOperations rbind -m ${type1}.point.plot_matrix_mRNA.gz ${type2}.point.plot_matrix_mRNA.gz -o ${type}.merged.point.plot_matrix_mRNA.gz
plotProfile -m ${type}.merged.point.plot_matrix_mRNA.gz -out mean.${type}.point.plot_matrix_mRNA_profile.pdf --plotHeight 7 --plotWidth 7 --averageType mean    #--perGroup 
plotProfile -m ${type}.merged.point.plot_matrix_mRNA.gz -out median.${type}.point.plot_matrix_mRNA_profile.pdf --plotHeight 7 --plotWidth 7 --averageType median    #--perGroup 
plotHeatmap -m ${type}.merged.point.plot_matrix_mRNA.gz -out Spectral.${type}.point.merged_heatmap.pdf --outFileSortedRegions sorted.${type}.bed --sortRegions descend --sortUsing mean --heatmapWidth 1.3 --heatmapHeight 5 --colorMap Spectral --sortUsing mean --zMax 150 --whatToShow 'heatmap and colorbar'

type="METTL3"
type1=${type}_fwd
type2=${type}_rev
computeMatrix reference-point -p 20 --scoreFileName $Rloop_fwd -R ${type}.fwd.bed --referencePoint center -b 2000 -a 2000 --binSize 10 --skipZeros -o ${type1}.point.plot_matrix_mRNA.gz --outFileNameMatrix ${type1}.point.loadR_matrix.tab --averageTypeBins mean
computeMatrix reference-point -p 20 --scoreFileName $Rloop_rev -R ${type}.rev.bed --referencePoint center -b 2000 -a 2000 --binSize 10 --skipZeros -o ${type2}.point.plot_matrix_mRNA.gz --outFileNameMatrix ${type2}.point.loadR_matrix.tab --averageTypeBins mean
computeMatrixOperations rbind -m ${type1}.point.plot_matrix_mRNA.gz ${type2}.point.plot_matrix_mRNA.gz -o ${type}.merged.point.plot_matrix_mRNA.gz
plotProfile -m ${type}.merged.point.plot_matrix_mRNA.gz -out mean.${type}.point.plot_matrix_mRNA_profile.pdf --plotHeight 7 --plotWidth 7 --averageType mean    #--perGroup 
plotProfile -m ${type}.merged.point.plot_matrix_mRNA.gz -out median.${type}.point.plot_matrix_mRNA_profile.pdf --plotHeight 7 --plotWidth 7 --averageType median    #--perGroup 
plotHeatmap -m ${type}.merged.point.plot_matrix_mRNA.gz -out Spectral.${type}.point.merged_heatmap.pdf --outFileSortedRegions sorted.${type}.bed --sortRegions descend --sortUsing mean --heatmapWidth 1.3 --heatmapHeight 5 --colorMap Spectral --sortUsing mean --zMax 150 --whatToShow 'heatmap and colorbar'


####Figure2H
cat $shCTRL_DDX21 $shMETTL3_DDX21 | sort -k1,1 -k2,2n | bedtools merge -i - | intersectBed -a - -b $R_loop -wa -wb | cut -f1-3 | awk '!a[$0]++' > DDX21.bed
cat $shCTRL_METTL3 $shDDX21_METTL3 | sort -k1,1 -k2,2n | bedtools merge -i - | intersectBed -a - -b $R_loop -wa -wb | cut -f1-3 | awk '!a[$0]++' > METTL3.bed

multiBamSummary BED-file -p 12 --BED DDX21.bed -o DDX21.npz --outRawCounts DDX21.tab --bamfiles \
./shCTRL_DDX21/no_MT.bam \
./shMETTL3_DDX21/no_MT.bam 

multiBamSummary BED-file -p 12 --BED METTL3.bed -o METTL3.npz --outRawCounts METTL3.tab --bamfiles \
./shCTRL_METTL3/no_MT.bam \
./shDDX21_METTL3/no_MT.bam 

sed '1d' DDX21.tab | sed 's/\.0//g' | sed 's/^/chr/' | awk -v OFS="\t" '{print $1,$2,$3,$4*1e9/("'$shCTRL_DDX21'"*($3-$2)),$5*1e9/("'$shMETTL3_DDX21'"*($3-$2))}' | awk '$4>0 && $5>0' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$5/$4}' > DDX21.FC
sed '1d' METTL3.tab | sed 's/\.0//g' | sed 's/^/chr/' | awk -v OFS="\t" '{print $1,$2,$3,$4*1e9/("'$shCTRL_METTL3'"*($3-$2)),$5*1e9/("'$shDDX21_METTL3'"*($3-$2))}' | awk '$4>0 && $5>0' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$5/$4}' > METTL3.FC

options(stringsAsFactors=F)
DDX21 <- read.table("DDX21.FC",header=F)
METTL3 <- read.table("METTL3.FC",header=F)
cols <- c("SteelBlue","DarkOrange")
library(ggplot2)
DDX21.dat <- data.frame(coverage=c(DDX21$V4,DDX21$V5),type=rep(c("DDX21_shCTRL","DDX21_shMETTL3"),each=nrow(DDX21)))
METTL3.dat <- data.frame(coverage=c(METTL3$V4,METTL3$V5),type=rep(c("METTL3_shCTRL","METTL3_shDDX21"),each=nrow(METTL3)))
data <- rbind(DDX21.dat,METTL3.dat)
data$type <- factor(data$type,level=c("DDX21_shCTRL","DDX21_shMETTL3","METTL3_shCTRL","METTL3_shDDX21"))
cols <- c("#999999", "#E69F00", "#999999", "#E69F00")	
library(ggplot2)
wilcox.test(log2(subset(data,type=="DDX21_shCTRL")$coverage),log2(subset(data,type=="DDX21_shMETTL3")$coverage),alternative="greater")
#p=0.154
wilcox.test(log2(subset(data,type=="METTL3_shCTRL")$coverage),log2(subset(data,type=="METTL3_shDDX21")$coverage),alternative="greater")
#p<2.2e-16
pdf("DDX21_METTL3_coverage.pdf",width=3.1,height=3.8)
ggplot(data, aes(x=type,y=log2(coverage)))+
  geom_boxplot(aes(fill = type), outlier.shape=NA,width=0.42)+
  #geom_boxplot(width = 0.5,outlier.shape = NA,trim = FALSE,aes(fill = sample))+
  ylim(c(0,20)) +
  labs(x="",y=expression(paste(log[2],"(coverage)")),title="") +
  #ylim(c(0,8)) +
  scale_fill_manual(values=cols)+
  theme_classic() + theme(
  legend.title=element_blank(),
  axis.title=element_text(size=16),
  axis.text.x=element_text(size=14,angle=30,hjust=1,vjust=1,color="Black"),
  axis.text.y=element_text(size=14,color="Black"),
  legend.position="none",axis.ticks.x=element_blank()) 
dev.off()



