####clean data
cutadapt -a AGATCGGAAGAGCACACGTCT -A AGATCGGAAGAGCGTCGTGTAG -o 1-trim_1.fastq -p 1-trim_2.fastq raw1.fastq raw2.fastq
cutadapt -u 10 -u -10 -U 10 -U -10 -o 2-qua_1.fastq -p 2-qua_2.fastq 1-trim_1.fastq 1-trim_2.fastq
java -jar trimmomatic-0.33.jar PE -phred33 ./2-qua_1.fastq ./2-qua_2.fastq ./3-qua_1.fastq ./unpaired_1.fastq ./3-qua_2.fastq ./unpaired_2.fastq LEADING:1 TRAILING:1 SLIDINGWINDOW:4:1 MINLEN:50
fastqc raw1.fastq -t 12 -o fastqc/
fastqc raw2.fastq -t 12 -o fastqc/
fastqc 1-trim_1.fastq -t 12 -o fastqc/
fastqc 1-trim_2.fastq -t 12 -o fastqc/
fastqc 2-qua_1.fastq -t 12 -o fastqc/
fastqc 2-qua_2.fastq -t 12 -o fastqc/
fastqc 3-qua_1.fastq -t 12 -o fastqc/
fastqc 3-qua_2.fastq -t 12 -o fastqc/

####map to genome
mkdir bowtie2_dup
cd bowtie2_dup
bowtie2 --phred33 -p 12 -x $hs -1 ../clean_data/3-qua_1.fastq -2 ../clean_data/3-qua_2.fastq -S tmp.sam 2> map.log
samtools view -bS -1 -q 20 -h tmp.sam | samtools sort -@ 10 -m 3G -l 9 -o sort.bam -
#保留tail需要加--local参数 -1, use fast BAM compression 
#-m, Set maximum memory per thread; -l Set compression level
rm tmp.sam
samtools index sort.bam
java -jar picard.jar MarkDuplicates REMOVE_DUPLICATES=true METRICS_FILE=sort.matrix INPUT=./sort.bam OUTPUT=sorted.nodup.bam ASSUME_SORTED=true
samtools index sorted.nodup.bam
rm sort.bam
samtools view -h sorted.nodup.bam | awk '$1~/^@/||$3!="MT"' | samtools view -Sb - | samtools sort -@ 10 -o no_MT.bam -
rm sorted.nodup.bam

####separate strand
mkdir forward
cd forward
samtools view -b -f 128 -F 16 ../no_MT.bam > fwd1.bam
samtools index fwd1.bam
samtools view -b -f 80 ../no_MT.bam > fwd2.bam
samtools index fwd2.bam
samtools merge -f fwd.bam fwd1.bam fwd2.bam
samtools index fwd.bam
rm fwd1.bam fwd2.bam


cd ../
mkdir reverse
cd reverse
samtools view -b -f 144 ../no_MT.bam > rev1.bam
samtools index rev1.bam
samtools view -b -f 64 -F 16 ../no_MT.bam > rev2.bam
samtools index rev2.bam
samtools merge -f rev.bam rev1.bam rev2.bam
samtools index rev.bam
rm rev1.bam rev2.bam

bamCoverage --normalizeUsing RPKM -b fwd.bam --binSize 20 -o fwd.bw
bamCoverage --normalizeUsing RPKM -b rev.bam --binSize 20 -o rev.bw


####call peak
fwd=fwd.bam
rev=rev.bam
samtools flagstat $fwd > fwd_bam.flagstat
samtools flagstat $rev > rev_bam.flagstat

macs2 callpeak -t $fwd -g hs -f BAM --nomodel --keep-dup all -n fwd --broad --broad-cutoff 0.01
macs2 callpeak -t $rev -g hs -f BAM --nomodel --keep-dup all -n rev --broad --broad-cutoff 0.01
awk -v OFS="\t" '{print "chr"$1,$2,$3,"*","*","+"}' fwd_peaks.broadPeak > R_loop.bed
awk -v OFS="\t" '{print "chr"$1,$2,$3,"*","*","-"}' rev_peaks.broadPeak >> R_loop.bed

echo "GeneID Chr Start End Strand" | sed 's/ /\t/g' > fwd.saf
awk '$6=="+"' R_loop.bed | cut -f1-3,6 | sort -k1,1 -k2,2n - | awk -v OFS="\t" '{print "peak_"NR,$0}'  >> fwd.saf
featureCounts -T 6 -F 'SAF' -a fwd.saf --minOverlap 1 --fraction -O -o fwd.enrichment.bed $fwd

echo "GeneID Chr Start End Strand" | sed 's/ /\t/g' > rev.saf
awk '$6=="-"' R_loop.bed | cut -f1-3,6 | sort -k1,1 -k2,2n - | awk -v OFS="\t" '{print "peak_"NR,$0}'  >> rev.saf
featureCounts -T 6 -F 'SAF' -a rev.saf --minOverlap 1 --fraction -O -o rev.enrichment.bed $rev
sed '1,2d' fwd.enrichment.bed > R_loop.enrichment.bed
sed '1,2d' rev.enrichment.bed >> R_loop.enrichment.bed
rm R_loop.bed fwd.saf rev.saf fwd.enrichment.bed rev.enrichment.bed
awk -v OFS="\t" '{print $2,$3,$4,$6,$7,$5}' R_loop.enrichment.bed | awk -v OFS="\t" '{print $0,($5*1e9)/($4*"'$total_reads'")}' > R_loop.RPKM.bed







