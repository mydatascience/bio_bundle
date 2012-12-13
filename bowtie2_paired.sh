#!/bin/bash
#USAGE ./bowtie2_paired.sh reference reads1 reads2 result threads hashsz 
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
echo "build started" | tee -a $out/log 
date | tee -a $out/log
bowtie2-build $1 $out/$ref_base.ind && \
echo "mapping started" | tee -a $out/log && \
date | tee -a $out/log && \
bowtie2 -t -q -p $threads -a --fr --minins 0 --maxins 500 $out/$ref_base.ind -1 $2 -2 $3 > $out/$ref_base.sam && \
echo "convert sam to bam started" | tee -a $out/log && \
date | tee -a $out/log && \
samtools view -bt $1 $out/$ref_base.sam > $out/$ref_base.bam && \
echo "sort and index bam started" | tee -a $out/log && \
date | tee -a $out/log && \
samtools sort $out/$ref_base.bam $out/$ref_base.sorted && \
samtools index $out/$ref_base.sorted.bam && \
echo "building vcf started" | tee -a $out/log && \
date | tee -a $out/log && \
samtools faidx $1 && \
samtools mpileup -uf $1 $out/$ref_base.sorted.bam | bcftools view -cvg - > $out/$ref_base.vcf 
echo "done" | tee -a $out/log
date | tee -a $out/log
samtools flagstat $out/$ref_base.sorted.bam | tee -a $out/log) | tee $out/log.full
