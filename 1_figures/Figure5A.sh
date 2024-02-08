options(stringsAsFactors=F)
m6A <- read.table("m6A.gene",header=F)[[1]]
NET <- read.table("net_total_gene_readthrough.txt",header=F)

NET$FC_DDX21 <- NET$V8/NET$V7
NET$FC_METTL3 <- NET$V9/NET$V7
NET$logFC_DDX21 <- log2(NET$V8/NET$V7)
NET$logFC_METTL3 <- log2(NET$V9/NET$V7)
m6A.NET <- subset(NET,V4 %in% m6A)
non_m6A.NET <- subset(NET,!(V4 %in% m6A))

library(ggplot2)
data <- data.frame(count=c(1856,72,1427,30),sample=rep(c("shDDX21","shMETTL3"),each=2),type=rep(c("up","down"),2))
data$sample <- factor(data$sample,level=c("shMETTL3","shDDX21"))
data$type <- factor(data$type,level=c("up","down"))
cols <- c("SandyBrown","#444c59")
#cols <- c("#e9a3c9","#a1d76a")
pdf("m6A_gene_readthrough_result.pdf",height=3.8,width=4.2)
ggplot(data)+
  geom_bar(aes(x=sample,y=count,fill=type),color="black",width=0.6,stat="identity",position=position_dodge())+
  ylim(0,2000) + 
  #scale_y_continuous(labels=scales::percent) +
  scale_fill_manual(values=cols) +		#c("#BEBEBE","#00BFC4","#F8766D")
  theme(title=element_text(size=14),axis.title.y=element_text(size=15))+
  theme(axis.text.x = element_text(size = 14, color = "black"),legend.text=element_text(size=12,color="black"))+ 
  theme(axis.text.y = element_text(size = 14, color = "black"),legend.title=element_blank())+
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),panel.border = element_rect(color="Black",fill=NA))+
  labs(x="",y="Number of genes") 
dev.off()


