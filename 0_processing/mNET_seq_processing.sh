####clean data
cutadapt -a AGATCGGAAG  -o 1-trim_1.fq.gz  $raw_1
fastqc $raw_1 -t 6 -o fastqc/
cutadapt -a AAAAAAAAAA -e 0.1 -o 2-trim_1.fq.gz 1-trim_1.fq.gz
fastqc 1-trim_1.fq.gz -t 6 -o fastqc/
rm 1-trim_1.fq.gz
cutadapt -u 4 -o 3-trim_1.fq.gz 2-trim_1.fq.gz
fastqc 3-trim_1.fq.gz -t 6 -o fastqc/
rm 2-trim_1.fq.gz
java -Xmx4g -jar trimmomatic-0.33.jar SE -phred33 ./3-trim_1.fq.gz ./4-trim_1.fq.gz LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:10
fastqc 4-trim_1.fq.gz -t 6 -o fastqc/
rm 3-trim_1.fq.gz


####map to genome
tophat -a 10 -p 6 -G ${refGene} -o ./ ${ref_genome} ../clean_data/4-trim_1.fq.gz

samtools view -bS -1 -q 20 -h accepted_hits.bam | samtools sort -@ 10 -m 3G -l 9 -o sort.bam - 
samtools index sort.bam
mkdir bams
cd bams
for chrom in `seq 1 22` X
do
samtools view -@ 6 -bh ../sort.bam chr${chrom} > chr${chrom}.bam
samtools index chr${chrom}.bam
done
samtools merge no_MT.bam chr1.bam chr2.bam chr3.bam chr4.bam chr5.bam chr6.bam chr7.bam chr8.bam chr9.bam chr10.bam chr11.bam chr12.bam chr13.bam chr14.bam chr15.bam chr16.bam chr17.bam chr18.bam chr19.bam chr20.bam chr21.bam chr22.bam chrX.bam

####separate strand
samtools view -b -F 16 no_MT.bam > fwd.bam
samtools view -b -f 16 no_MT.bam > rev.bam
samtools index fwd.bam
samtools index rev.bam

####readthrough calculation
multiBamSummary BED-file -p 12 --BED fwd_terminator.bed -o fwd_terminator.npz --outRawCounts fwd_terminator.tab --bamfiles \
./shCTRL-NET-seq-1/tophat2/fwd.bam \
./shCTRL-NET-seq-2/tophat2/fwd.bam \
./shDDX21-NET-seq-1/tophat2/fwd.bam \
./shDDX21-NET-seq-2/tophat2/fwd.bam \
./shMETTL3-NET-seq-1/tophat2/fwd.bam \
./shMETTL3-NET-seq-2/tophat2/fwd.bam
sed '1d' fwd_terminator.tab | sed 's/\.0//g' | sort -k1,1 -k2,2 -k3,3n | paste fwd_terminator.bed - | awk -v OFS="\t" '$1==$7 && $2==$8 && $3==$9' - | cut -f1-6,10-15 > fwd_terminator_count.txt
####
multiBamSummary BED-file -p 12 --BED rev_terminator.bed -o rev_terminator.npz --outRawCounts rev_terminator.tab --bamfiles \
./shCTRL-NET-seq-1/tophat2/rev.bam \
./shCTRL-NET-seq-2/tophat2/rev.bam \
./shDDX21-NET-seq-1/tophat2/rev.bam \
./shDDX21-NET-seq-2/tophat2/rev.bam \
./shMETTL3-NET-seq-1/tophat2/rev.bam \
./shMETTL3-NET-seq-2/tophat2/rev.bam 
sed '1d' rev_terminator.tab|sed 's/\.0//g' | sort -k1,1 -k2,2 -k3,3n | paste rev_terminator.bed - | awk -v OFS="\t" '$1==$7 && $2==$8 && $3==$9' - | cut -f1-6,10-15 > rev_terminator_count.txt

cat fwd_terminator_count.txt rev_terminator_count.txt > terminator_count.txt
rm fwd_terminator_count.txt rev_terminator_count.txt


####readthrough level calculation
awk -v OFS="\t" '{if($6=="+") print $1,$3-1,$3,$4,$5,$6}' $ref|bedtools slop -r 10000 -l -1 -i - -g hg19.sizes > fwd_downstream.bed
awk -v OFS="\t" '{if($6=="+") print $1,$2,$3,$4,$5,$6}' $ref|bedtools slop -r 10000 -l 0 -i - -g hg19.sizes > fwd_extension.bed
intersectBed -a fwd_downstream.bed -b fwd_extension.bed -s -wa -wb|awk -v OFS="\t" '$4!=$10' -|cut -f1-6|sort|uniq > fwd_overlapped.bed

awk -v OFS="\t" '{if($6=="-") print $1,$2,$2+1,$4,$5,$6}' $ref|bedtools slop -r -1 -l 10000 -i - -g hg19.sizes > rev_downstream.bed
awk -v OFS="\t" '{if($6=="-") print $1,$2,$3,$4,$5,$6}' $ref|bedtools slop -r 0 -l 10000 -i - -g hg19.sizes > rev_extension.bed
intersectBed -a rev_downstream.bed -b rev_extension.bed -s -wa -wb|awk -v OFS="\t" '$4!=$10' -|cut -f1-6|sort|uniq > rev_overlapped.bed

cat fwd_overlapped.bed rev_overlapped.bed > overlapped_10000.bed
rm fwd_downstream.bed fwd_extension.bed rev_downstream.bed rev_extension.bed

awk -v OFS="\t" 'NR==FNR{a[$4]=$0}NR>FNR{if(!a[$4]) print $0}' overlapped_10000.bed $ref > mRNA_10000.bed
awk '$6=="+"' mRNA_10000.bed > fwd_mRNA_10000.bed
awk '$6=="-"' mRNA_10000.bed > rev_mRNA_10000.bed

awk -v OFS="\t" '$3-$2 >= 5000' fwd_mRNA_10000.bed > fwd_long_mRNA_10000.bed
awk -v OFS="\t" '$3-$2 >= 5000' rev_mRNA_10000.bed > rev_long_mRNA_10000.bed


cat $fwd_genebody | awk -v OFS="\t" '{print $1,$3-1,$3,$4,$5,$6}' - | bedtools slop -r 6000 -l 0 -i - -g $chromsize | bedtools slop -r 0 -l 0 -i - -g $chromsize | sort -k1,1 -k2,2 -k3,3n > fwd_terminator.bed
cat $rev_genebody | awk -v OFS="\t" '{print $1,$2,$2+1,$4,$5,$6}' - | bedtools slop -r 0 -l 6000 -i - -g $chromsize | bedtools slop -r 0 -l 0 -i - -g $chromsize | sort -k1,1 -k2,2 -k3,3n > rev_terminator.bed

cat $fwd_genebody | bedtools slop -r -500 -l -500 -i - -g $chromsize | sort -k1,1 -k2,2 -k3,3n > fwd_gene.bed
cat $rev_genebody | bedtools slop -r -500 -l -500 -i - -g $chromsize | sort -k1,1 -k2,2 -k3,3n > rev_gene.bed

#genebody count
###########################################################################################################################
multiBamSummary BED-file -p 12 --BED fwd_gene.bed -o fwd_gene.npz --outRawCounts fwd_gene.tab --bamfiles \
./shCTRL-NET-seq-1/fwd.bam \
./shCTRL-NET-seq-2/tophat2/fwd.bam \
.shDDX21-NET-seq-1/tophat2/fwd.bam \
./shDDX21-NET-seq-2/tophat2/fwd.bam \
./shMETTL3-NET-seq-1/tophat2/fwd.bam \
./shMETTL3-NET-seq-2/tophat2/fwd.bam
sed '1d' fwd_gene.tab | sed 's/\.0//g' | sort -k1,1 -k2,2 -k3,3n | paste fwd_gene.bed - | awk -v OFS="\t" '$1==$7 && $2==$8 && $3==$9' - | cut -f1-6,10-15 > fwd_gene_count.txt
####
multiBamSummary BED-file -p 12 --BED rev_gene.bed -o rev_gene.npz --outRawCounts rev_gene.tab --bamfiles \
./shCTRL-NET-seq-1/tophat2/rev.bam \
./shCTRL-NET-seq-2/tophat2/rev.bam \
./shDDX21-NET-seq-1/tophat2/rev.bam \
./shDDX21-NET-seq-2/tophat2/rev.bam \
./shMETTL3-NET-seq-1/tophat2/rev.bam \
./shMETTL3-NET-seq-2/tophat2/rev.bam 
sed '1d' rev_gene.tab|sed 's/\.0//g' | sort -k1,1 -k2,2 -k3,3n | paste rev_gene.bed - | awk -v OFS="\t" '$1==$7 && $2==$8 && $3==$9' - | cut -f1-6,10-15 > rev_gene_count.txt

cat fwd_gene_count.txt rev_gene_count.txt > gene_count.txt
rm fwd_gene_count.txt rev_gene_count.txt


awk -v OFS="\t" 'NR==FNR{ct1=$1;ct2=$2;shd1=$3;shd3=$4;shm1=$5;shm2=$6}NR>FNR{print $1,$2,$3,$4,$5,$6,$7*1e9/(ct1*($3-$2)),$8*1e9/(ct2*($3-$2)),$9*1e9/(shd1*($3-$2)),$10*1e9/(shd3*($3-$2)),$11*1e9/(shm1*($3-$2)),$12*1e9/(shm2*($3-$2))}' total_reads.txt terminator_count.txt | sort -k4,4 > RPKM_terminator_count.txt

awk -v OFS="\t" 'NR==FNR{ct1=$1;ct2=$2;shd1=$3;shd3=$4;shm1=$5;shm2=$6}NR>FNR{print $1,$2,$3,$4,$5,$6,$7*1e9/(ct1*($3-$2)),$8*1e9/(ct2*($3-$2)),$9*1e9/(shd1*($3-$2)),$10*1e9/(shd3*($3-$2)),$11*1e9/(shm1*($3-$2)),$12*1e9/(shm2*($3-$2))}' total_reads.txt gene_count.txt | sort -k4,4 > RPKM_gene_count.txt


i=0.2 
paste RPKM_gene_count.txt RPKM_terminator_count.txt | awk '$4==$16' | cut -f1-12,19-24 | awk '$7>"'$i'"&&$8>"'$i'"&&$9>"'$i'"&&$10>"'$i'"&&$11>"'$i'"&&$12>"'$i'"&&$13>"'$i'"&&$14>"'$i'"&&$15>"'$i'"&&$16>"'$i'"&&$17>"'$i'"&&$18>"'$i'"' >  net_total_gene.txt
awk -v OFS="\t" '{shCTRL=($13+$14)/($7+$8);shDDX21=($15+$16)/($9+$10);shMETTL3=($17+$18)/($11+$12);print $1,$2,$3,$4,$5,$6,shCTRL,shDDX21,shMETTL3}' net_total_gene.txt > termination_readthrough_ref_genebody.net_total_gene_readthrough.txt


