####clean data
read1=$rawdata/*1.fastq.gz
read2=$rawdata/*2.fastq.gz
mkdir ${sample}
cd ${sample}

mkdir clean_data
cd clean_data
fastqc $read1 -t 4 -o ./
fastqc $read2 -t 4 -o ./
#############cutadapter
cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG -o cutadapt_1.fq.gz -p cutadapt_2.fq.gz $read1 $read2
fastqc cutadapt_1.fq.gz -t 4 -o ./
fastqc cutadapt_2.fq.gz -t 4 -o ./
#############trim low quality reads
java -Xmx4g -jar trimmomatic-0.33.jar PE -phred33 ./cutadapt_1.fq.gz ./cutadapt_2.fq.gz ./trim_1.fq.gz ./unpaired_1.fq.gz ./trim_2.fq.gz ./unpaired_2.fq.gz LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:18;
fastqc trim_1.fq.gz -t 4 -o ./
fastqc trim_2.fq.gz -t 4 -o ./
rm cutadapt_1.fq.gz cutadapt_2.fq.gz

####map to genome
hisat2 -p 6 --dta -x $ref -1 ../clean_data/trim_2.fq.gz -2 ../clean_data/trim_2.fq.gz -S hisat2.sam
samtools view -@ 4 -Sh hisat2.sam -q 20 > hisat2_20.sam
grep -E "^@|NH:i:1$|NH:i:1[^0-9]" hisat2_20.sam > uniqmap.sam  
samtools view -S uniqmap.sam -b > uniqmap.bam
samtools sort -@ 4 uniqmap.bam > sort.bam	
samtools index sort.bam 
rm hisat2.sam uniqmap.sam uniqmap.bam hisat2_20.sam 

samtools view -h sort.bam | awk '$1~/^@/||$3!="MT"' | samtools view -Sb - | samtools sort -@ 10 -o no_MT.bam -
samtools index no_MT.bam
samtools flagstat no_MT.bam > readsCount
rm sort.bam

####Separate strand
mkdir forward
cd forward
samtools view -b -f 128 -F 16 ../../sort_rm_duplicates.bam > fwd1.bam
samtools index fwd1.bam
samtools view -b -f 80 ../../sort_rm_duplicates.bam > fwd2.bam
samtools index fwd2.bam
samtools merge -f fwd.bam fwd1.bam fwd2.bam
samtools index fwd.bam
rm fwd1.bam fwd2.bam fwd1.bam.bai fwd2.bam.bai


cd ../
mkdir reverse
cd reverse
samtools view -b -f 144 ../../sort_rm_duplicates.bam > rev1.bam
samtools index rev1.bam
samtools view -b -f 64 -F 16 ../../sort_rm_duplicates.bam > rev2.bam
samtools index rev2.bam
samtools merge -f rev.bam rev1.bam rev2.bam
samtools index rev.bam
rm rev1.bam rev2.bam rev1.bam.bai rev2.bam.bai

####call peaks
treatment=$dir/fwd.bam
control=$dir/fwd.bam
macs2 callpeak -t $treatment -c $control -f BAM --nomodel -g hs --keep-dup all -q 0.01 -n fwd --extsize 200

treatment=$dir/rev.bam
control=$dir/rev.bam
macs2 callpeak -t $treatment -c $control -f BAM --nomodel -g hs --keep-dup all -q 0.01 -n rev --extsize 200

####filter peaks with miCLIP sites
intersectBed -a $miCLIP -b combined.enrichment.txt -wa -wb| awk '$4==$8' | cut -f1-4,13,14 | awk '!a[$0]++' | awk -v OFS="\t" '{print $0,$6/$5}' > site.filtered.combined.enrichment.bed

####calculate m6A enrichment
echo "GeneID Chr Start End Strand" | sed 's/ /\t/g' > rev.saf
sort -k1,1 -k2,2n rev_combined.bed | sed 's/^chr//' | awk -v OFS="\t" '{print "peak_"NR,$0}'  >> rev.saf
featureCounts -T 6 -F 'SAF' -a rev.saf --minOverlap 30 --fraction -O -o rev_combined.enrichment.bed \
./shCTRL-IP/rev.bam \
./shCTRL-Input/rev.bam \
./shDDX21-IP/rev.bam \
./shDDX21-Input/rev.bam

echo "GeneID Chr Start End Strand" | sed 's/ /\t/g' > fwd.saf
sort -k1,1 -k2,2n rev_combined.bed | sed 's/^chr//' | awk -v OFS="\t" '{print "peak_"NR,$0}'  >> fwd.saf
featureCounts -T 6 -F 'SAF' -a fwd.saf --minOverlap 30 --fraction -O -o rev_combined.enrichment.bed \
./shCTRL-IP/fwd.bam \
./shCTRL-Input/fwd.bam \
./shDDX21-IP/fwd.bam \
./shDDX21-Input/fwd.bam







