####clean data
gunzip -c  $dir*_1.fq.gz > 1.fq
gunzip -c  $dir*_2.fq.gz > 2.fq
fastqc 1.fq -t 6 -o ./
fastqc 2.fq -t 6 -o ./

#####Trimming of 3' linker sequences
fastx_clipper -a AGATCGGAAGAGCACACG -l 24 -n -i 1.fq -Q 33 | fastq_quality_trimmer -t 5 -l 24 -Q 33 -o trim_1.fastq
fastx_clipper -a AGATCGGAAGAGCGTCGT -l 24 -n -i 2.fq -Q 33 | fastq_quality_trimmer -t 5 -l 24 -Q 33 -o trim_2.fastq
rm 1.fq 2.fq
fastqc trim_1.fastq -t 6 -o ./
fastqc trim_2.fastq -t 6 -o ./

awk '{if(NR%4==1){flag=0;tmp=$1"#1"} if(NR%4==2 && length($0)<150){print(tmp);flag=1} if(flag==1){print($0)}}' trim_1.fastq > trim-f_1.fastq
awk '{if(NR%4==1){flag=0;tmp=$1"#2"} if(NR%4==2 && length($0)<150){print(tmp);flag=1} if(flag==1){print($0)}}' trim_2.fastq > trim-f_2.fastq
rm trim_1.fastq trim_2.fastq
fastqc trim-f_1.fastq -t 6 -o ./
fastqc trim-f_2.fastq -t 6 -o ./

fastx_reverse_complement -Q 33 -i trim-f_2.fastq -o trim-f_2r.fastq 
cat trim-f_1.fastq trim-f_2r.fastq > trim-all.fastq
rm trim-f_2.fastq trim-f_2r.fastq trim-f_1.fastq
fastqc trim-all.fastq -t 6 -o ./

perl $ctk/fastq2collapse.pl trim-all.fastq deduplicate_total.fastq 
rm trim-all.fastq
fastqc deduplicate_total.fastq -t 6 -o ./

cutadapt -a "A{10}" -e 0.1 -o rmPOLYA_total.fastq deduplicate_total.fastq
awk '{if(NR%4==1){split($1,a,"#");print a[1]"#"a[2]}else{print $0}}' rmPOLYA_total.fastq > temp1.fastq  ##去除末尾的polyA尾
rm deduplicate_total.fastq
fastqc rmPOLYA_total.fastq -t 6 -o ./ 

perl $ctk/stripBarcode.pl -format fastq -len 3 temp1.fastq debarcode_total.fastq
rm rmPOLYA_total.fastq temp1.fastq
fastqc debarcode_total.fastq -t 6 -o ./

java -Xmx4g -jar trimmomatic-0.33.jar SE -phred33 debarcode_total.fastq qf_trim_total.fastq LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:18
rm debarcode_total.fastq
fastqc qf_trim_total.fastq -t 6 -o ./


####map to genome
bwa aln -t 6 -n 0.06 -q 20 $hs ../clean_data/qf_trim_total.fastq > final.sai
bwa samse $hs final.sai ../clean_data/qf_trim_total.fastq > final.sam
samtools view -S final.sam -b -o final.bam
/software/biosoft/software/bamtools/tools/bin/bamtools stats -in final.bam > bwa_mapping_report
bedtools bamtobed -i final.bam -split > final.bed
gzip ../clean_data/qf_trim_total.fastq

#####CTK Processing
mkdir CTK_Procedure
cd CTK_Procedure
perl $ctk/parseAlignment.pl -v --map-qual 1 --min-len 18 --mutation-file final.mutation.txt ../final.sam final.tag.bed 
rm ../final.sam ../final.sai

perl $ctk/tag2collapse.pl -v -big --random-barcode -EM 30 --seq-error-model alignment -weight --weight-in-name --keep-max-score --keep-tag-name final.tag.bed final.tag.uniq.bed 
perl $ctk/selectRow.pl -q 3 -f 3 final.mutation.txt final.tag.uniq.bed > final.tag.uniq.mutation.txt 
rm final.tag.bed final.mutation.txt

perl $ctk/bed2rgb.pl -v -col "128,0,0" final.tag.uniq.bed final.tag.uniq.rgb.bed 

awk '{print $3-$2}' final.tag.uniq.bed | sort -n | uniq -c | awk '{print $2"\t"$1}' > final.uniq.len.dist.txt
awk '{if($9=="+") {print $0}}' final.tag.uniq.mutation.txt | cut -f 1-6 > final.tag.uniq.ins.bed
awk '{if($9==">") {print $0}}' final.tag.uniq.mutation.txt | cut -f 1-6 > final.tag.uniq.sub.bed
awk '{if($9=="-") {print $0}}' final.tag.uniq.mutation.txt | cut -f 1-6 > final.tag.uniq.del.bed



mkdir CIMS
cd CIMS
perl $ctk/CIMS.pl -n 20 -p -v --keep-cache -c cache_mut ../final.tag.uniq.bed ../final.tag.uniq.sub.bed final.tag.uniq.mut.CIMS.txt
cut -f 1-6 final.tag.uniq.mut.CIMS.txt > final.tag.uniq.mut.CIMS.bed
sed '1d' final.tag.uniq.mut.CIMS.txt|awk -v OFS="\t" '$8>=1 && $8/$7>=0.01 && $8/$7<=0.5{if($6=="+"){print $1,$2-4,$3+4,$4,$5,$6}else{print $1,$2-4,$3+4,$4,$5,$6}}'|cut -f1-6| awk '{if($2 >= 0){print $0}}' |fastaFromBed -fi $hs -bed - -s -fo final.tag.uniq.mut.fasta
awk -v OFS="\t" 'NR%2==1{tmp=$1}NR%2==0{if($1~/[G,A][G,A]AC[A,T,C]/){a=$1; gsub("[G,A][G,A]AC[A,T,C]", "\t", a); split(a,b,"\t"); split(b[1],c,""); l=length(c); print tmp, l, $1}}' final.tag.uniq.mut.fasta|sed 's/>//'|sed 's/:/\t/'|sed 's/-/\t/'|sed 's/(/\t/'|sed 's/)//'|awk -v OFS="\t" '$4~/+/{print $1,$2+$5, $2+$5+5, "mut_"$1"_"$2+4"_"$3-4"_"$4, "mut_"$6, $4}$4~/-/{print $1, $3-$5-5, $3-$5, "mut_"$1"_"$2+4"_"$3-4"_"$4, "mut_"$6, $4}' > final.mut.bed
cat final.mut.bed|awk -vOFS="\t" '{print $1,$2,$3,"*","*",$6}'|fastaFromBed -fi $hs -bed - -s -fo tmp.fasta
cat tmp.fasta|sed s'/T/U/g' > motif.fa
rm tmp.fasta


####CIMS_CT
cd ../
mkdir CIMS_CT
cd CIMS_CT
awk '($8=="C" && $10=="T" && $6=="+") || ($8=="G" && $10=="A" && $6=="-")' ../final.tag.uniq.mutation.txt |cut -f 1-6 > final.tag.uniq.sub-CT.bed


#####Mutation Mode
perl $ctk/CIMS.pl -n 10 -p -v --keep-cache -c cache_mut-CT ../final.tag.uniq.bed final.tag.uniq.sub-CT.bed final.tag.uniq.mut-CT.CIMS.txt
sed '1d' final.tag.uniq.mut-CT.CIMS.txt|awk '$8>=1 && $8/$7>=0.01 && $8/$7<=0.5'|cut -f1-6 > filter_final.tag.uniq.mut-CT.CIMS.txt
awk -v OFS="\t" '{if($6=="+"){print $1,$2-3,$3+1,$4,$5,$6}else{print $1,$2-3,$3+1,$4,$5,$6}}' filter_final.tag.uniq.mut-CT.CIMS.txt | awk '$2>=0' | fastaFromBed -fi $hs -bed - -s -fo ./filter_final.tag.uniq.mut-CT.fasta
awk '{if(NR%2==1){tmp=$0}else{split($1,a,"");if(a[3]=="A"){print tmp; print $0}}}' filter_final.tag.uniq.mut-CT.fasta > final.tag.uniq.mut-CT.fasta
sed 's/>//' final.tag.uniq.mut-CT.fasta|sed 's/:/\t/'|sed 's/-/\t/'|sed 's/(/\t/' |sed 's/)//'|awk -v OFS="\t" '{if(NR%2==1){chr=$1; chrStart=$2; chrEnd=$3; strand=$4;}else{print chr, chrStart, chrEnd, $0, strand}}' > final.tag.uniq.mut-CT.bed
