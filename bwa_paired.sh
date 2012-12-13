#!/bin/bash
#USAGE ./bwa.sh reference reads1 reads2 result
threads=4
hashsz=10 
if [ $# -ge 5 ]; then
    threads=$5
fi
if [ $# = 6 ]; then
    hashsz=$6
fi
ref_base=`echo $1 | perl -pe "s/^.*?([^\/]+?)(\.[^\/]+)?$/\1/; s/\n//" `
reads1_base=`echo $2 | perl -pe "s/^.*?([^\/]+?)(\.[^\/]+)?$/\1/; s/\n//" `
reads2_base=`echo $3 | perl -pe "s/^.*?([^\/]+?)(\.[^\/]+)?$/\1/; s/\n//" `
out=$4
(echo "treads number = $threads; hashsize = $hashsz" | tee $out/log 
echo "build index started" | tee -a $out/log
date | tee -a $out/log
bwa index -a bwtsw $1 && \
echo "align started" | tee -a $out/log && \
date | tee -a $out/log && \
bwa aln -t $threads $1 $2 > $out/$reads1_base.sai && \
bwa aln -t $threads $1 $3 > $out/$reads2_base.sai && \
echo "making sam started" | tee -a $out/log && \
date | tee -a $out/log && \
bwa sampe $1 $out/$reads1_base.sai $out/$reads2_base.sai $2 $3 > $out/$ref_base.sam && \
echo "convert sam to bam started" | tee -a $out/log && \
date | tee -a $out/log && \
samtools view -bt $1 $out/$ref_base.sam > $out/$ref_base.bam && \
echo "sort and index bam started" | tee -a $out/log && \
date | tee -a $out/log && \
samtools sort $out/$ref_base.bam $out/$ref_base.sorted && \
samtools index $out/$reads_base.sorted.bam && \
echo "building vcf started" | tee -a $out/log && \
date | tee -a $out/log && \
samtools faidx $1 && \
samtools index $out/$ref_base.sorted.bam && \
samtools mpileup -uf $1 $out/$ref_base.sorted.bam | bcftools view -vcg - > $out/$ref_base.vcf 
echo "done" | tee -a $out/log
date | tee -a $out/log
samtools flagstat $out/$reads_base.sorted.bam | tee -a $out/log) | tee $out/log.full
