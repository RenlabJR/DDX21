#####FigureS5ABC
computeMatrix scale-regions -p 20 --scoreFileName $shCTRL_fwd $shMETTL3_fwd $shDDX21_fwd --regionsFileName ../fwd.${m6A_type}.mRNA.bed --outFileName ${type1}.${m6A_type}.region.plot_matrix_mRNA.gz --outFileNameMatrix ${type1}.${m6A_type}.region.loadR_matrix.tab --regionBodyLength 5000 --startLabel TSS --endLabel TES --beforeRegionStartLength 2000 --afterRegionStartLength 10000 --binSize $bin_size --averageTypeBins mean --samplesLabel shCTRL shMETTL3 shDDX21 --skipZeros
computeMatrix scale-regions -p 20 --scoreFileName $shCTRL_rev $shMETTL3_rev $shDDX21_rev --regionsFileName ../rev.${m6A_type}.mRNA.bed --outFileName ${type2}.${m6A_type}.region.plot_matrix_mRNA.gz --outFileNameMatrix ${type2}.${m6A_type}.region.loadR_matrix.tab --regionBodyLength 5000 --startLabel TSS --endLabel TES --beforeRegionStartLength 2000 --afterRegionStartLength 10000 --binSize $bin_size --averageTypeBins mean --samplesLabel shCTRL shMETTL3 shDDX21 --skipZeros
computeMatrixOperations rbind -m ${type1}.${m6A_type}.region.plot_matrix_mRNA.gz ${type2}.${m6A_type}.region.plot_matrix_mRNA.gz -o ${m6A_type}.merged.region.m6A_mRNA.mat.gz
zcat total.merged.region.m6A_mRNA.mat.gz | sed '1d' > total.merged.region.m6A_mRNA.mat.txt

genelist <- m6A #or nom6A
type="y6A.m" #or n6A.m
shCTRL <- input.matrix[[1]][genelist, ]
shMETTL3 <- input.matrix[[2]][genelist, ]
shDDX21 <- input.matrix[[3]][genelist, ]
#"black", "yellow"
#"white", "SlateBlue"
shCTRL <- (shCTRL/shCTRL.region_rowMeans[rownames(shCTRL)])
shMETTL3 <- (shMETTL3/shMETTL3.region_rowMeans[rownames(shMETTL3)])
shDDX21 <- (shDDX21/shDDX21.region_rowMeans[rownames(shDDX21)])
shCTRL <- shCTRL[order(rowMeans(shCTRL[21:170]), decreasing=T), ]
shMETTL3 <- shMETTL3[rownames(shCTRL), ]
shDDX21 <- shDDX21[rownames(shCTRL), ]
colnames(shCTRL) <- colnames(shMETTL3) <- colnames(shDDX21) <- c(rep("", 19), "TSS", rep("", 49), "TES", rep("", 100))

########################################shCtrl
pdf(paste(type, "shCTRL_heatmap_scale_genebody.SlateBlue.pdf"), width=3.1, height=3.2)
col_fun <- colorRamp2(
  c(0,  6), 
  c("white", "SlateBlue")
)
Heatmap(shCTRL, col = col_fun, cluster_rows=F, cluster_columns=F, show_row_names=F, show_column_names=T, name = paste(type, "shCTRL_sg"), border = T) #na_col = "black", name = "expression", 
dev.off()
}
########################################shMETTL3
pdf(paste(type, "shMETTL3_heatmap_scale_genebody.SlateBlue.pdf"), width=3.1, height=3.2)
col_fun <- colorRamp2(
  c(0,  6), 
  c("white", "SlateBlue")
)
Heatmap(shMETTL3, col = col_fun, cluster_rows=F, cluster_columns=F, show_row_names=F, show_column_names=T, name = paste(type, "shMETTL3_sg"), border = T) #na_col = "black", name = "expression", 
dev.off()
}
########################################shDDX21
pdf(paste(type, "shDDX21_heatmap_scale_genebody.SlateBlue.pdf"), width=3.1, height=3.2)
col_fun <- colorRamp2(
  c(0,  6), 
  c("white", "SlateBlue")
)
Heatmap(shDDX21, col = col_fun, cluster_rows=F, cluster_columns=F, show_row_names=F, show_column_names=T, name = paste(type, "shDDX21_sg"), border = T) #na_col = "black", name = "expression", 
dev.off()



type <- "m6A"
options(stringsAsFactors=F)
library(RColorBrewer)   
cols <- brewer.pal(9,"Set1")[c(2,4)]
data <- read.table("termination_readthrough_ref_genebody.net_total_gene_readthrough.txt", header=F, sep="\t")

data$FC_shD21 <- (data$V8/data$V7)
data$FC_shMETTL3 <- (data$V9/data$V7)
data$log2FC_shD21 <- log2((data$V8)/(data$V7))
data$log2FC_shMETTL3 <- log2((data$V9)/(data$V7))

m6A.data <- subset(data, V4 %in% m6A)
no_m6A.data <- subset(data, !(V4 %in% nom6A))
m6A_FC <- m6A.data$log2FC_shMETTL3
nom6A_FC <- no_m6A.data$log2FC_shMETTL3
wilcox.test(m6A_FC, nom6A_FC, alternative="greater")
summary(m6A_FC)
summary(nom6A_FC)

m6A_FC <- m6A.data$log2FC_shD21
nom6A_FC <- no_m6A.data$log2FC_shD21
wilcox.test(m6A_FC, nom6A_FC, alternative="greater")
summary(m6A_FC)
summary(nom6A_FC)
#shMETTL3
pdf("termination_6kb_readthrough_ref_genebody_m6A.pdf",width=4.5,height=4.5)
m6A_FC <- m6A.data$log2FC_shMETTL3
nom6A_FC <- no_m6A.data$log2FC_shMETTL3
tmp.data <- m6A.data[c("V7","V9","FC_shMETTL3")]
names(tmp.data) <- c("shCtrl","shtreat","FC")
tmp.data$shCtrl <- log2(tmp.data$shCtrl+1)
tmp.data$shtreat <- log2(tmp.data$shtreat+1)

boxplot(nom6A_FC,m6A_FC,ylim=c(-1,2),outline=F,xaxt="n",ylab="log2(shMETTL3/shCTRL)",main="readthrough",col=cols)
axis(1,at=c(1,2),labels=paste(type, c("-", "+"), sep=""));
#shDDX21
m6A_FC <- m6A.data$log2FC_shD21
nom6A_FC <- no_m6A.data$log2FC_shD21
tmp.data <- m6A.data[c("V7","V8","FC_shD21")]
names(tmp.data) <- c("shCtrl","shtreat","FC")
tmp.data$shCtrl <- log2(tmp.data$shCtrl+1)
tmp.data$shtreat <- log2(tmp.data$shtreat+1)

boxplot(nom6A_FC,m6A_FC,ylim=c(-2,3),outline=F,xaxt="n",ylab="log2(shDDX21/shCTRL)",main="readthrough",col=cols)
axis(1,at=c(1,2),labels=paste(type, c("-", "+"), sep=""));
dev.off()

#####FigureS5DE
setwd("/xtdisk/yangyg_group/liumx/On_DDX21/Figurelayout/20240223/FigureS5DE")
library(RColorBrewer)	##display.brewer.all(type="div")
cols <- brewer.pal(9,"Set1")
m6A <- read.table("m6A_anno_mRNA.bed", header=F)
library(plyr)
reshape.m6A <- ddply(m6A, .(V12, V11), m6A_level=sum(V5), summarise)
#Rloop <- read.table("Rloop.gene", header=F)[[1]]
gene_rpkm <- read.table("RPKM_gene_count.txt", header=F)
gene_rpkm$shCTRL <- rowMeans(gene_rpkm[7:8])

quantiles <- c(2, 5)
m6A.1 <- subset(reshape.m6A, m6A_level <= quantiles[1] )
m6A.2 <- subset(reshape.m6A, m6A_level > quantiles[1] & m6A_level <= quantiles[2])
m6A.3 <- subset(reshape.m6A, m6A_level > quantiles[2])
summary(m6A.1$m6A_level)
summary(m6A.2$m6A_level)
summary(m6A.3$m6A_level)
m6A.1 <- unique(m6A.1$V11)
m6A.2 <- unique(m6A.2$V11)
m6A.3 <- unique(m6A.3$V11)
length(m6A.1)
length(m6A.2)
length(m6A.3)

#gene readthrough fold changes calculation
data$FC_shD21 <- (data$V8/data$V7)
data$FC_shMETTL3 <- (data$V9/data$V7)
data$log2FC_shD21 <- log2((data$V8)/(data$V7))
data$log2FC_shMETTL3 <- log2((data$V9)/(data$V7))
m6A.1.FC <- subset(data, !(V4 %in% c(m6A.2, m6A.3)))
#m6A.1.FC <- subset(data, V4 %in% m6A.1)
m6A.2.FC <- subset(data, V4 %in% m6A.2)
m6A.3.FC <- subset(data, V4 %in% m6A.3)

tmp.data <- rbind(m6A.1.FC, m6A.2.FC, m6A.3.FC)
tmp.data$type <- c(rep("low/-", nrow(m6A.1.FC)), rep("medium", nrow(m6A.2.FC)), rep("high", nrow(m6A.3.FC)))
tmp.data <- merge(tmp.data, reshape.m6A, by.x="V4", by.y="V11", all.x=T)
tmp.data[is.na(tmp.data)] <- "-"
write.table(tmp.data, "mNET_m6A_level_low_medium_high.txt", row.names=F, col.names=T, sep="\t", quote=F)

m6A.1.rpkm <- subset(gene_rpkm, (V4 %in% m6A.1.FC$V4))$shCTRL
m6A.2.rpkm <- subset(gene_rpkm, (V4 %in% m6A.2.FC$V4))$shCTRL
m6A.3.rpkm <- subset(gene_rpkm, (V4 %in% m6A.3.FC$V4))$shCTRL
length(m6A.1.rpkm)
length(m6A.2.rpkm)
length(m6A.3.rpkm)
summary(m6A.1.rpkm)
summary(m6A.2.rpkm)
summary(m6A.3.rpkm)
wilcox.test(m6A.1.rpkm, m6A.2.rpkm)
wilcox.test(m6A.1.rpkm, m6A.3.rpkm)
wilcox.test(m6A.2.rpkm, m6A.3.rpkm)
boxplot(m6A.1.rpkm, m6A.2.rpkm, m6A.3.rpkm, outline=F,xaxt="n",ylab="rpkm",main="",col=cols)
axis(1,at=c(1,2,3),labels=paste("m6A", c("low/-", "medium", "high"), sep=""));

pdf("m6A_level_readthrough_boxplot.pdf")
summary(m6A.1.FC$log2FC_shMETTL3)
summary(m6A.2.FC$log2FC_shMETTL3)
summary(m6A.3.FC$log2FC_shMETTL3)
length(m6A.1.FC$log2FC_shMETTL3)
length(m6A.2.FC$log2FC_shMETTL3)
length(m6A.3.FC$log2FC_shMETTL3)
wilcox.test(m6A.1.FC$log2FC_shMETTL3, m6A.2.FC$log2FC_shMETTL3, alternative="less") #p-value = 2.044e-07
wilcox.test(m6A.1.FC$log2FC_shMETTL3, m6A.3.FC$log2FC_shMETTL3, alternative="less") #p-value < 2.2e-16
wilcox.test(m6A.2.FC$log2FC_shMETTL3, m6A.3.FC$log2FC_shMETTL3, alternative="less") #p-value = 7.844e-05
boxplot(m6A.1.FC$log2FC_shMETTL3, m6A.2.FC$log2FC_shMETTL3, m6A.3.FC$log2FC_shMETTL3, outline=F,xaxt="n",ylab="log2(shMETTL3/shCTRL)",main="mNET readthrough",col=cols)
axis(1,at=c(1,2,3),labels=paste("m6A", c("low/-", "medium", "high"), sep=""));

summary(m6A.1.FC$log2FC_shD21)
summary(m6A.2.FC$log2FC_shD21)
summary(m6A.3.FC$log2FC_shD21)
length(m6A.1.FC$log2FC_shD21)
length(m6A.2.FC$log2FC_shD21)
length(m6A.3.FC$log2FC_shD21)
wilcox.test(m6A.1.FC$log2FC_shD21, m6A.2.FC$log2FC_shD21, alternative="less")
wilcox.test(m6A.1.FC$log2FC_shD21, m6A.3.FC$log2FC_shD21, alternative="less")
wilcox.test(m6A.2.FC$log2FC_shD21, m6A.3.FC$log2FC_shD21, alternative="less")
boxplot(m6A.1.FC$log2FC_shD21, m6A.2.FC$log2FC_shD21, m6A.3.FC$log2FC_shD21, outline=F,xaxt="n",ylab="log2(shDDX21/shCTRL)",main="mNET readthrough",col=cols)
axis(1,at=c(1,2,3),labels=paste("m6A", c("low/-", "medium", "high"), sep=""));
dev.off()




#####FigureS5FG
multiBamSummary BED-file -p 12 --BED DDX21.bed -o DDX21.npz --outRawCounts DDX21.tab --bamfiles \
./DDX21/no_MT.bam \
./METTL3/no_MT.bam

multiBamSummary BED-file -p 12 --BED METTL3.bed -o METTL3.npz --outRawCounts METTL3.tab --bamfiles \
./DDX21/no_MT.bam \
./METTL3/no_MT.bam
sed '1d' DDX21.tab | sed 's/^/chr/' | sed 's/\.0//g' | intersectBed -a DDX21.bed -b - -wa -wb -f 1 -F 1 | cut -f1-6,10-13 > DDX21.txt
sed '1d' METTL3.tab | sed 's/^/chr/' | sed 's/\.0//g' | intersectBed -a METTL3.bed -b - -wa -wb -f 1 -F 1 | cut -f1-6,10-13 > METTL3.txt

#
awk -v OFS="\t" '{D21="'${DDX21}'"; M3="'${METTL3}'"; print $1, $2, $3, $4, $5, $6, 1e6*$7/D21, 1e6*$8/M3}' DDX21.txt > RPKM_norm.DDX21.txt
awk -v OFS="\t" '{D21="'${DDX21}'"; M3="'${METTL3}'"; print $1, $2, $3, $4, $5, $6, 1e6*$7/D21, 1e6*$8/M3}' METTL3.txt > RPKM_norm.METTL3.txt
intersectBed -a RPKM_norm.DDX21.txt -b mNET_tophat.termination_readthrough_ref_genebody.net_total_gene_readthrough.txt -wa -wb -s > DDX21.mNET_tophat.termination_readthrough_ref_genebody.net_total_gene_readthrough.txt
intersectBed -a RPKM_norm.METTL3.txt -b mNET_tophat.termination_readthrough_ref_genebody.net_total_gene_readthrough.txt -wa -wb -s > METTL3.mNET_tophat.termination_readthrough_ref_genebody.net_total_gene_readthrough.txt


options(stringAsFactors=F)
data <- read.table("mNET_tophat.termination_readthrough_ref_genebody.net_total_gene_readthrough.txt", header=F)
data$FC_shD21 <- (data$V8/data$V7)
data$FC_shMETTL3 <- (data$V9/data$V7)
data$log2FC_shD21 <- log2((data$V8)/(data$V7))
data$log2FC_shMETTL3 <- log2((data$V9)/(data$V7))

library(plyr)
DDX21.mNET <- read.table("DDX21.mNET_tophat.termination_readthrough_ref_genebody.net_total_gene_readthrough.txt", header=F)
METTL3.mNET <- read.table("METTL3.mNET_tophat.termination_readthrough_ref_genebody.net_total_gene_readthrough.txt", header=F)
DDX21.mNET$DDX21_enrichment <- DDX21.mNET$V7
METTL3.mNET$METTL3_enrichment <- METTL3.mNET$V8
DDX21.mNET <- ddply(DDX21.mNET, .(V12), DDX21_sum_enrichment=sum(DDX21_enrichment), summarise)
METTL3.mNET <- ddply(METTL3.mNET, .(V12), METTL3_sum_enrichment=sum(METTL3_enrichment), summarise)
names(DDX21.mNET)[1] <- names(METTL3.mNET)[1] <- "V4"
data <- merge(data, DDX21.mNET, by="V4", all.x=T)
data <- merge(data, METTL3.mNET, by="V4", all.x=T)
m6A <- read.table("m6A.gene", header=F)[[1]]
library(RColorBrewer)   
cols <- rev(brewer.pal(9,"Set1"))

pdf("total_gene_readthrough_correlation_with_DDX21_METTL3_binding strength.pdf")
data[is.na(data)] <- 0
data <- subset(data, METTL3_sum_enrichment>0)
data <- subset(data, DDX21_sum_enrichment>0)
data$METTL3_binding_type[data$METTL3_sum_enrichment>0] <- "METTL3+ weak" #quantile50%
data$METTL3_binding_type[data$METTL3_sum_enrichment>10] <- "METTL3+ moderate" #quantile 50%-70%
data$METTL3_binding_type[data$METTL3_sum_enrichment>30] <- "METTL3+ strong"
table(data$METTL3_binding_type)

wilcox.test(subset(data, METTL3_binding_type=="METTL3+ moderate")$FC_shMETTL3, subset(data, METTL3_binding_type=="METTL3+ weak")$FC_shMETTL3, alternative="greater")
wilcox.test(subset(data, METTL3_binding_type=="METTL3+ strong")$FC_shMETTL3, subset(data, METTL3_binding_type=="METTL3+ moderate")$FC_shMETTL3, alternative="greater")
wilcox.test(subset(data, METTL3_binding_type=="METTL3+ moderate")$FC_shD21, subset(data, METTL3_binding_type=="METTL3+ weak")$FC_shD21, alternative="greater")
wilcox.test(subset(data, METTL3_binding_type=="METTL3+ strong")$FC_shD21, subset(data, METTL3_binding_type=="METTL3+ moderate")$FC_shD21, alternative="greater")

data$METTL3_binding_type <- factor(data$METTL3_binding_type, level=c("METTL3+ weak", "METTL3+ moderate", "METTL3+ strong"))
table(data$METTL3_binding_type)
boxplot(log2FC_shMETTL3 ~ METTL3_binding_type, data,ylim=c(-1,2),outline=F,ylab="log2(shMETTL3/shCTRL)",xlab="",main="readthrough",col=cols)
boxplot(log2FC_shD21 ~ METTL3_binding_type, data,ylim=c(-2,3),outline=F,ylab="log2(shDDX21/shCTRL)",xlab="",main="readthrough",col=cols)


data[is.na(data)] <- 0
data$DDX21_binding_type <- rep("DDX21-", nrow(data))
data$DDX21_binding_type[data$DDX21_sum_enrichment>0] <- "DDX21+ weak" #quantile50%
data$DDX21_binding_type[data$DDX21_sum_enrichment>20] <- "DDX21+ moderate" #quantile 50%-70%
data$DDX21_binding_type[data$DDX21_sum_enrichment>55] <- "DDX21+ strong"
data$DDX21_binding_type <- factor(data$DDX21_binding_type, level=c("DDX21+ weak", "DDX21+ moderate", "DDX21+ strong"))
table(data$DDX21_binding_type)
wilcox.test(subset(data, DDX21_binding_type=="DDX21+ moderate")$FC_shMETTL3, subset(data, DDX21_binding_type=="DDX21+ weak")$FC_shMETTL3, alternative="greater")
wilcox.test(subset(data, DDX21_binding_type=="DDX21+ strong")$FC_shMETTL3, subset(data, DDX21_binding_type=="DDX21+ moderate")$FC_shMETTL3, alternative="greater")
wilcox.test(subset(data, DDX21_binding_type=="DDX21+ moderate")$FC_shD21, subset(data, DDX21_binding_type=="DDX21+ weak")$FC_shD21, alternative="greater")
wilcox.test(subset(data, DDX21_binding_type=="DDX21+ strong")$FC_shD21, subset(data, DDX21_binding_type=="DDX21+ moderate")$FC_shD21, alternative="greater")

boxplot(log2FC_shMETTL3 ~ DDX21_binding_type, data,ylim=c(-1,2),outline=F,ylab="log2(shMETTL3/shCTRL)",xlab="",main="readthrough",col=cols)
boxplot(log2FC_shD21 ~ DDX21_binding_type, data,ylim=c(-2,3),outline=F,ylab="log2(shDDX21/shCTRL)",xlab="",main="readthrough",col=cols)
dev.off()


#####FigureS5H
Rloop.NET <- subset(NET, (V4 %in% Rloop))
Rloop_m6A.NET <- subset(NET, (V4 %in% intersect(Rloop, m6A)))
#log2FC METTL3
wilcox.test(Rloop.NET$logFC_DDX21,Rloop_m6A.NET$logFC_DDX21, alternative="less")
#log2FC METTL3
wilcox.test(Rloop.NET$logFC_METTL3,Rloop_m6A.NET$logFC_METTL3, alternative="less")
tmp.dat <- rbind(Rloop.NET, Rloop_m6A.NET)
tmp.dat$type <- c(rep("Rloop+", nrow(Rloop.NET)), rep("Rloop+m6A+", nrow(Rloop_m6A.NET)))
tmp.dat$type <- factor(tmp.dat$type, level=c("Rloop+", "Rloop+m6A+"))
boxplot(logFC_DDX21~type, tmp.dat,ylim=c(-2,3),outline=F,ylab="log2(shDDX21/shCTRL)",main="readthrough",col=cols)
axis(1,at=c(1,2));
boxplot(logFC_METTL3~type, tmp.dat,ylim=c(-2,2),outline=F,ylab="log2(shMETTL3/shCTRL)",main="readthrough",col=cols)

#####FigureS5I
options(stringAsFactors=F)
NET <- read.table("termination_readthrough_ref_genebody.net_total_gene_readthrough.txt", header=F, sep="\t")
NET$FC_DDX21 <- (NET$V8/NET$V7)
NET$FC_METTL3 <- (NET$V9/NET$V7)
NET$logFC_DDX21 <- log2((NET$V8)/(NET$V7))
NET$logFC_METTL3 <- log2((NET$V9)/(NET$V7))
m6A <- read.table("m6A.gene", header=F)[[1]]
m6A.NET <- subset(NET, V4 %in% m6A)

library(ggsci)
library(ggplot2)
cols <- pal_nejm()(9)
new.m6A.NET <- data.frame()
sorted.idx <- sort(m6A.NET$logFC_DDX21, decreasing=F, index.return=T)$ix
for(i in 1:29){
  new.m6A.NET <- rbind(new.m6A.NET, colMeans(m6A.NET[sorted.idx[((i-1)*94+1):(i*94)],c(10:13)]))
}
i <- 30
new.m6A.NET <- rbind(new.m6A.NET, colMeans(m6A.NET[sorted.idx[((i-1)*94+1):(nrow(m6A.NET))],c(10:13)]))
colnames(new.m6A.NET) <- names(m6A.NET[c(10:13)])
cor(m6A.NET$logFC_DDX21, m6A.NET$logFC_METTL3, method="spearman")

cols <- pal_nejm()(9)
min(new.m6A.NET$logFC_DDX21); max(new.m6A.NET$logFC_DDX21); 
min(new.m6A.NET$logFC_METTL3); max(new.m6A.NET$logFC_METTL3); 
ggplot(new.m6A.NET, aes(x = logFC_DDX21, y = logFC_METTL3)) + 
  geom_point(size = 1.5, color="black") + 
  #ylim(0,1e4)+
  scale_x_continuous(breaks = seq(-0.5,1.5,0.5))+
  scale_y_continuous(breaks = seq(-0.5,1,0.5))+
  #xlim(-1,2)+
  #ylim(-0.4,1.2)+
  geom_smooth(method = lm, se = TRUE, fill = cols[2]) +  #se=TRUE,添加95%置信区间
  #scale_color_brewer(palette = "Set1") + 
  #scale_x_continuous(breaks = seq(12,18,1)) + 
  labs(x = "log2FC(shD21/shCtrl)", y = "log2FC(shM3/shCtrl)") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5, size = 16), 
        plot.caption = element_text(size = 12), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 15)) 
dev.off()



