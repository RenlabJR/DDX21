####Figure3C
options(stringsAsFactors=F)
data <- read.table(m6A_file,header=F)
colnames(data) <- c("chr","start","end","strand","shCTRL","shDDX21","FC")
cols <- c("grey","#ef8a62","#67a9cf")	#stable,up,down

shCTRL.x <- knots(ecdf(data[,5]+seq(1e-5,by=1e-12,length.out=length(data[,5]))));
shCTRL.y <-(0:(length(data[,5])-1))/(length(data[,5])-1);

shDDX21.x <- knots(ecdf(data[,6]+seq(1e-5,by=1e-12,length.out=length(data[,6]))));
shDDX21.y <-(0:(length(data[,6])-1))/(length(data[,6])-1);

pdf("m6A_shDDX21_differential_cumulative_curve.pdf",width=4.5,height=4.5)
plot(shCTRL.x, shCTRL.y, type="l",col=cols[1], main="", xlab="m6A level",ylab="Cumulative frequency", bty="n",xaxt="n",yaxt="n",lwd=1,xlim=c(0,6));	#
lines(shDDX21.x, shDDX21.y,col=cols[2],lwd=1);
axis(1,pos=0,at=seq(0,6,3),lwd=2);
axis(2,pos=0,at=c(0,0.5,1),lwd=2);
legend(-1,1,c("shCTRL","shDDX21"),col=cols,bty="n",lty=rep(1,9),lwd=2,merge=TRUE);
dev.off()


####Figure3D
#bash code
bedtools slop -i $mRNA -g hg19.sizes -b 3000 | intersectBed -a - -b combine.middle_cluster.bed -wa -wb -F 0.5 -s| awk -v OFS="\t" 'NR==FNR{x[$4]=$2"\t"$3}NR>FNR{if(x[$4]) print $1,x[$4],$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}' $mRNA - | awk -v OFS="\t" '{if ($6=="+" && $13<$2) print $0,3000-($2-$13); 
else if ($6=="+" && $13>=$2 && $14<=$3) print $0,3000+5000*($13-$2)/($3-$2); 
else if ($6=="+" && $14>$3) print $0,3000+5000+($14-$3); 
else if ($6=="-" && $14>$3) print $0,3000-($14-$3); 
else if ($6=="-" && $13>=$2 && $14<=$3) print $0,3000+5000*($3-$14)/($3-$2); 
else if ($6=="-" && $13<$2) print $0,3000+5000+($2-$13)}' - > combine.enrichment.txt

####R code
combine <- read.table("combine.enrichment.txt",header=F,stringsAsFactors=F)
library(plyr)
combine.dat <- data.frame(pos=round(combine$V15/11000*110),FoldChange=combine$V11/combine$V10)

down <- subset(combine.dat,FoldChange<1/1.3)
down.dat <- data.frame(pos=as.numeric(names(table(down$pos))),count=as.numeric(table(down$pos)))
up <- subset(combine.dat,FoldChange>1.3)
up.dat <- data.frame(pos=as.numeric(names(table(up$pos))),count=as.numeric(table(up$pos)))
normal <- subset(combine.dat,FoldChange>=1/1.3 & FoldChange<=1.3)
normal.dat <- data.frame(pos=as.numeric(names(table(normal$pos))),count=as.numeric(table(normal$pos)))

down.dat$density <- down.dat$count/sum(down.dat$count)
up.dat$density <- up.dat$count/sum(up.dat$count)
normal.dat$density <- normal.dat$count/sum(normal.dat$count)

down.x <- smooth.spline(down.dat$pos,down.dat$density,df=30)$x
down.y <- smooth.spline(down.dat$pos,down.dat$density,df=30)$y
up.x <- smooth.spline(up.dat$pos,up.dat$density,df=30)$x
up.y <- smooth.spline(up.dat$pos,up.dat$density,df=30)$y
normal.x <- smooth.spline(normal.dat$pos,normal.dat$density,df=30)$x
normal.y <- smooth.spline(normal.dat$pos,normal.dat$density,df=30)$y

down.label <- paste0("A Hypo-peaks ",nrow(down))
up.label <- paste0("A Hyper-peaks ",nrow(up))
normal.label <- paste0("A stable-peaks ",nrow(normal))

cols <- c("SteelBlue","FireBrick","Black")
pdf("mRNA_gene_peak_density.pdf",width=4.2,height=3.8)
plot(down.x,down.y,type='l',xaxt="n",yaxt="n",xlab="",ylab="",bty="n",col=cols[1],ylim=c(0,0.1),lwd=2)
lines(up.x,up.y,col=cols[2],lwd=2)
lines(normal.x,normal.y,col=cols[3],lwd=2)
legend("topleft",c(expression(paste(m^6, "A Hypo-peaks 2126")),expression(paste(m^6, "A Hyper-peaks 1096")) ,expression(paste(m^6, "A stable-peaks 4954"))),col=cols,lty=1,lwd=2,cex=1,box.lty=0)
axis(1,at=c(0,30,80,110),label=c("-3kb","TSS","TES","3kb"),lwd=2,cex=1.3)
axis(2,at=c(0,0.05,0.1),label=c(0,0.05,0.1),lwd=2,las=1)	#las=1,
mtext("Density",side=2,las=0,line=2.2,cex=1.3)
dev.off()


####Figure3E
library(plyr)
options(stringsAsFactors=F)
CoRloop <- read.table("Rloop.m6A.bed",header=F,stringsAsFactors=F)
CoRloop <- ddply(CoRloop,.(V1, V2, V3, V4, V5, V6, V7, V8),R_loop=sum(V9), summarise)

CoRloop$type <- rep("low",nrow(CoRloop))
CoRloop$type[CoRloop$R_loop>2] <- "high"
table(CoRloop$type)
p1 <- wilcox.test(subset(CoRloop,type=="high")$V7, subset(CoRloop,type=="low")$V7,alternative="less")$p.value
#data$type <- factor(data$type,level=c("None","low","high"))
cols<-c("#459943","#bdc3d2","#e8c559","#ea9c9d")
pdf(pdf_name,width=5.1,height=5.5)
boxplot(log2(V7) ~ type,data=CoRloop,col=cols,outline=F,boxwex = 0.45,ylab="m6A log2(shDDX21/shCTRL)",xlab="")
text(1,-1,paste0("p=",p1))
dev.off()






