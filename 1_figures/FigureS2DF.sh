####FigureS2D
computeMatrix reference-point -p 12 -R $DDX21 -S $METTL3_bw -b 1000 -a 1000 --referencePoint center --binSize 5 --skipZeros -o METTL3_2_DDX21_peak.gz --outFileSortedRegions METTL3_2_DDX21_peak.bed --outFileNameMatrix METTL3_2_DDX21_peak.tab
plotHeatmap -m METTL3_2_DDX21_peak.gz -out METTL3_2_DDX21_peak.pdf --outFileSortedRegions sorted.METTL3_2_DDX21_peak.bed --sortRegions descend --sortUsing mean  --heatmapWidth 1.5 --heatmapHeight 4  --colorList 'black, yellow' --whatToShow 'heatmap and colorbar' --samplesLabel METTL3_2_DDX21_peak

####FigureS2F
options(stringAsFactors=F)
Rloop <- read.table("mRNA.Rloop.peaks", header=F)
DDX21 <- read.table("DDX21_annoRloop.freq", header=F)[c("V8","V9")]
METTL3 <- read.table("METTL3_annoRloop.freq", header=F)[c("V8","V9")]
#DDX21
Rloop.DDX21 <- merge(Rloop, DDX21, by="V8", all.x=T)
Rloop.DDX21$V9[is.na(Rloop.DDX21$V9)] <- 0 
Rloop.DDX21 <- Rloop.DDX21[order(Rloop.DDX21$V7),]

quantiles <- quantile(Rloop.DDX21$V7, probs = seq(0, 1, by = 0.2))
Rloop.DDX21$group <- 1
Rloop.DDX21$group[Rloop.DDX21$V7>quantiles[2]] <- 2
Rloop.DDX21$group[Rloop.DDX21$V7>quantiles[3]] <- 3
Rloop.DDX21$group[Rloop.DDX21$V7>quantiles[4]] <- 4
Rloop.DDX21$group[Rloop.DDX21$V7>quantiles[5]] <- 5

sum(subset(Rloop.DDX21, group==1)$V9)
sum(subset(Rloop.DDX21, group==2)$V9)
sum(subset(Rloop.DDX21, group==3)$V9)
sum(subset(Rloop.DDX21, group==4)$V9)
sum(subset(Rloop.DDX21, group==5)$V9)

Rloop.DDX21$group <- as.character(Rloop.DDX21$group)
Rloop.DDX21$group <- factor(Rloop.DDX21$group, level=1:5)
library(RColorBrewer)	
library(ggsci)
cols = pal_npg("nrc", alpha = 0.7)(10)
library(ggplot2)
library(plyr)
library(RColorBrewer)	
display.brewer.all()
cols <- brewer.pal(6,"YlGnBu")

tmp.dat <- ddply(Rloop.DDX21, .(group), sum_count=sum(V9), count=sum(V9>0), summarise)
pdf("Rloop_peaks_5_groups_DDX21_binding_count.pdf", width=5.6, height=4.5)
ggplot(tmp.dat, aes(x=group, y=sum_count, fill=group)) + 
            geom_bar(stat = "identity", width=0.8, colour = "black",size=0.2)+ 
            scale_fill_manual(values=cols)+
            theme_classic()+
            ylab("DDX21 bingding peaks count")+
            xlab("R-loop peak group orderd by enrichment level")+
            theme(axis.text=element_text(size=12,colour = "black"),
                  axis.title.x=element_text(size=12), axis.title.y=element_text(size=12),legend.position="none")
dev.off()



