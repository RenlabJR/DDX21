####clean data
cutadapt -a AGATCGGAAG -o cutadapt_1.fastq.gz $raw_1
cutadapt -a AAAAAAAAAA -o cutadapt_2.fastq.gz  cutadapt_1.fastq.gz 
cutadapt -u 4 -o cutadapt_3.fastq.gz cutadapt_2.fastq.gz

java -Xmx4g -jar trimmomatic-0.33.jar SE -phred33 cutadapt_3.fastq.gz trim_1.fastq.gz LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:18
####
fastqc -t 6 $raw_1 -o ./fastqc 
fastqc -t 6 cutadapt_1.fastq.gz  -o ./fastqc 
fastqc -t 6 cutadapt_2.fastq.gz -o ./fastqc 
fastqc -t 6 cutadapt_3.fastq.gz  -o ./fastqc 
fastqc -t 6 trim_1.fastq.gz  -o ./fastqc 
rm cutadapt_*fastq.gz 

####map to genome
gunzip -c ../clean_data/trim_1.fastq.gz > trim.fq
bowtie $hs -v 2 -m 10 -p 12 --best --strata trim.fq ./mapping.sam
rm trim.fq

awk -v FS="\t" '$3!="MT"' ./mapping.sam| awk 'BEGIN{FS="\t";OFS="\t"}length($6) < 100 {$3="chr"$3;print $0}' > mapping_noMT.sam

echo 'mapped_reads:' > ct_conversion.info
wc -l mapping_noMT.sam >> ct_conversion.info
echo ''  >> ct_conversion.info

echo 'mutated reads:' >> ct_conversion.info
cat mapping_noMT.sam |cut -f8|awk '$1'|wc -l >> ct_conversion.info
echo ''  >> ct_conversion.info

echo 'mutated nucleotide:' >> ct_conversion.info
cat mapping_noMT.sam |cut -f8|awk '$1'|awk '{split($1,a,",");for(i=1;i<=length(a);i++){print a[i]}}'|awk -v FS=":" '{print $2}'|sort|uniq -c|sed -e 's/^ *//'|sed -e 's/ /\t/'  >> ct_conversion.info
echo ''  >> ct_conversion.info

echo 'conversion ratio:' >> ct_conversion.info
cat mapping_noMT.sam |cut -f8|awk '$1'|awk '{split($1,a,",");for(i=1;i<=length(a);i++){print a[i]}}'|awk -v FS=":" '{print $2}'|sort|uniq -c|sed -e 's/^ *//'|sed -e 's/ /\t/'|awk -v OFS="\t" '{num+=$1;tmp[$2]=$1}END{for(i in tmp){print i,tmp[i],num,tmp[i]/num}}'|sort  >> ct_conversion.info

####PARalyzer call cluster
sed -i "s?sam_path?$(pwd)?" PAR.ini
sed -i "s?output_path?$(pwd)?" PAR.ini
sh PARalyzer 6G PAR.ini

awk -v OFS="\t" '{print $1, $2, $3, "peak_"NR, "peak_"NR, $4}' shCTRL_DDX21.bed | intersectBed -a - -b $genebody -wa -wb -s | grep mRNA | cut -f1-6 | awk '!a[$0]++' > anno_mRNA.shCTRL_DDX21.bed
awk -v OFS="\t" '{print $1, $2, $3, "peak_"NR, "peak_"NR, $4}' shCTRL_METTL3.bed | intersectBed -a - -b $genebody -wa -wb -s | grep mRNA | cut -f1-6 | awk '!a[$0]++' > anno_mRNA.shCTRL_METTL3.bed



