#####FigureS5A
type <- "m6A"
options(stringsAsFactors=F)
library(RColorBrewer)   
cols <- brewer.pal(9,"Set1")[c(2,4)]
data <- read.table("termination_readthrough_ref_genebody.net_total_gene_readthrough.txt", header=F, sep="\t")
m6A <- read.table(paste(type,".gene", sep=""),header=F)[[1]]

data$FC_shD21 <- (data$V8/data$V7)
data$FC_shMETTL3 <- (data$V9/data$V7)
data$log2FC_shD21 <- log2((data$V8)/(data$V7))
data$log2FC_shMETTL3 <- log2((data$V9)/(data$V7))

m6A.data <- subset(data, V4 %in% m6A)
no_m6A.data <- subset(data, !(V4 %in% m6A))
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


#####FigureS5B
type <- "m6A"
options(stringsAsFactors=F)
library(RColorBrewer)   
cols <- brewer.pal(9,"Set1")[c(2,4)]
data <- read.table("upstream_readthrough_ref_genebody.net_total_gene_readthrough.txt", header=F, sep="\t")
m6A <- read.table(paste(type,".gene", sep=""),header=F)[[1]]

data$FC_shD21 <- (data$V8/data$V7)
data$FC_shMETTL3 <- (data$V9/data$V7)
data$log2FC_shD21 <- log2((data$V8)/(data$V7))
data$log2FC_shMETTL3 <- log2((data$V9)/(data$V7))

m6A.data <- subset(data, V4 %in% m6A)
no_m6A.data <- subset(data, !(V4 %in% m6A))
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
pdf("upstream_2kb_readthrough_ref_genebody_m6A.pdf",width=4.5,height=4.5)
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
multiBamSummary BED-file -p 12 --BED DDX21.bed -o DDX21.npz --outRawCounts DDX21.tab --bamfiles \
./DDX21/no_MT.bam \
./METTL3/no_MT.bam

multiBamSummary BED-file -p 12 --BED METTL3.bed -o METTL3.npz --outRawCounts METTL3.tab --bamfiles \
./DDX21/no_MT.bam \
./METTL3/no_MT.bam
sed '1d' DDX21.tab | sed 's/^/chr/' | sed 's/\.0//g' | intersectBed -a DDX21.bed -b - -wa -wb -f 1 -F 1 | cut -f1-6,10-13 > DDX21.txt
sed '1d' METTL3.tab | sed 's/^/chr/' | sed 's/\.0//g' | intersectBed -a METTL3.bed -b - -wa -wb -f 1 -F 1 | cut -f1-6,10-13 > METTL3.txt

#基因上的Peaks计算enrichment的和
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
data$METTL3_binding_type[data$METTL3_sum_enrichment>0] <- "METTL3+ weak"
data$METTL3_binding_type[data$METTL3_sum_enrichment>10] <- "METTL3+ moderate"
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
data$DDX21_binding_type[data$DDX21_sum_enrichment>0] <- "DDX21+ weak"
data$DDX21_binding_type[data$DDX21_sum_enrichment>10] <- "DDX21+ moderate"
data$DDX21_binding_type[data$DDX21_sum_enrichment>30] <- "DDX21+ strong"
data$DDX21_binding_type <- factor(data$DDX21_binding_type, level=c("DDX21+ weak", "DDX21+ moderate", "DDX21+ strong"))
table(data$DDX21_binding_type)
wilcox.test(subset(data, DDX21_binding_type=="DDX21+ moderate")$FC_shMETTL3, subset(data, DDX21_binding_type=="DDX21+ weak")$FC_shMETTL3, alternative="greater")
wilcox.test(subset(data, DDX21_binding_type=="DDX21+ strong")$FC_shMETTL3, subset(data, DDX21_binding_type=="DDX21+ moderate")$FC_shMETTL3, alternative="greater")
wilcox.test(subset(data, DDX21_binding_type=="DDX21+ moderate")$FC_shD21, subset(data, DDX21_binding_type=="DDX21+ weak")$FC_shD21, alternative="greater")
wilcox.test(subset(data, DDX21_binding_type=="DDX21+ strong")$FC_shD21, subset(data, DDX21_binding_type=="DDX21+ moderate")$FC_shD21, alternative="greater")

boxplot(log2FC_shMETTL3 ~ DDX21_binding_type, data,ylim=c(-1,2),outline=F,ylab="log2(shMETTL3/shCTRL)",xlab="",main="readthrough",col=cols)
boxplot(log2FC_shD21 ~ DDX21_binding_type, data,ylim=c(-2,3),outline=F,ylab="log2(shDDX21/shCTRL)",xlab="",main="readthrough",col=cols)
dev.off()