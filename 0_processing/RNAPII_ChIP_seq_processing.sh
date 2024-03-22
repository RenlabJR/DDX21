####clean data
raw_1=/$data_path/shCtrl/*1.fq.gz
raw_2=/$data_path/shCtrl/*2.fq.gz
/software/biosoft/software/fastp/fastp -w $thread -i $raw_1 -o ./trim_1.fq.gz -I $raw_2 -O ./trim_2.fq.gz --unpaired1 ./unpaired_1.fq.gz --unpaired2 ./unpaired_2.fq.gz --trim_poly_g --adapter_fasta /xtdisk/yangyg_group/liumx/reference/adapter_sum.fa -l 50 -W 4 -M 20
fastqc -t $thread -o ./fastqc $raw_1 $raw_2 trim_1.fq.gz trim_2.fq.gz
####map to genome hg19-mm10 mix
cd ../
mkdir bowtie2_dup
cd bowtie2_dup
bowtie2 --phred33 -p $thread -x $hs_mm -1 ../clean_data/trim_1.fq.gz -2 ../clean_data/trim_2.fq.gz -S tmp.sam 2> map.log
samtools view -bS -1 -q 20 -h tmp.sam | samtools sort -@ $thread -l 9 -o sort.bam -
rm tmp.sam
samtools index sort.bam
echo "q>20|rmPCR|noMT" > readsNumber
samtools view sort.bam | cut -f1 | awk '!a[$0]++' | wc -l >> readsNumber 
sambamba-0.8.1 markdup -r -t $thread sort.bam sorted.nodup.bam
samtools index sorted.nodup.bam
samtools view sorted.nodup.bam | cut -f1 | awk '!a[$0]++' | wc -l >> readsNumber 
#rm sort.bam
samtools flagstat sorted.nodup.bam > total.reads
bamCoverage -p 6 -bs 20 -b sorted.nodup.bam --normalizeUsing RPKM -o bamCoverage.RPKM.sorted.nodup.bw


#separate bam into human and mouse
echo "sort.nodup.bam->Separate human and mouse:" > sort.nodup.hs_mm.readsCount
samtools view -f 2 -bh $map_bam > paired.mapped.bam
samtools view paired.mapped.bam | cut -f1 | awk '!a[$0]++' | wc -l >> sort.nodup.hs_mm.readsCount

samtools view -H paired.mapped.bam > sort.nodup.hs.sam
samtools view paired.mapped.bam | awk -v OFS="\t" '$3 ~ /^H/' >> sort.nodup.hs.sam
samtools view -Sbh sort.nodup.hs.sam -o sort.nodup.hs.bam
samtools index sort.nodup.hs.bam

samtools view -H paired.mapped.bam > sort.nodup.mm.sam
samtools view paired.mapped.bam | awk -v OFS="\t" '$3 ~ /^M/' >> sort.nodup.mm.sam
samtools view -Sbh sort.nodup.mm.sam -o sort.nodup.mm.bam
samtools index sort.nodup.mm.bam
rm *sam *bedgraph

samtools view sort.nodup.hs.bam| cut -f1 | awk '!a[$0]++' | wc -l >> sort.nodup.hs_mm.readsCount
samtools view sort.nodup.mm.bam| cut -f1 | awk '!a[$0]++' | wc -l >> sort.nodup.hs_mm.readsCount

mouse_nodup_reads=`awk 'NR==4{print $1}' sort.nodup.hs_mm.readsCount`
scalefactor=`echo $mouse_nodup_reads | awk '{print 1e6/"'$mouse_nodup_reads'"}'`
bamCoverage --bam sort.nodup.hs.bam --outFileName RRPM_norm.bamcoverage.nodup.hs.bedgraph --binSize 20 --scaleFactor $scalefactor  --numberOfProcessors 12 --normalizeUsing None -of bedgraph
grep "^H" RRPM_norm.bamcoverage.nodup.hs.bedgraph | sed 's/^H//' | grep -v "chrMT" > RRPM_norm_excludeH.bamcoverage.nodup.hs.bedgraph
bedGraphToBigWig RRPM_norm_excludeH.bamcoverage.nodup.hs.bedgraph hg19.size RRPM_norm_excludeH.bamcoverage.nodup.hs.bw
rm *bedgraph


