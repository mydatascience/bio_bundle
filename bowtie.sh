#!/bin/bash
#USAGE ./bowtie.sh reference reads result threads hashsz 
threads=4
hashsz=10 
if [ $# -ge 4 ]; then
    threads=$4
fi
if [ $# = 5 ]; then
    hashsz=$5
fi
ref_base=`echo $1 | perl -pe "s/^.*?([^\/]+?)(\.[^\/]+)?$/\1/; s/\n//" `
reads_base=`echo $2 | perl -pe "s/^.*?([^\/]+?)(\.[^\/]+)?$/\1/; s/\n//" `
out=$3
(echo "treads number = $threads; hashsize = $hashsz" | tee $out/log 
echo "build started" | tee -a $out/log 
date | tee -a $out/log
bowtie-build $1 $out/$ref_base.ind && \
echo "mapping started" | tee -a $out/log && \
date | tee -a $out/log && \
bowtie -t -q -p $threads --best --sam $out/$ref_base.ind $2 > $out/$reads_base.sam && \
echo "convert sam to bam started" | tee -a $out/log && \
date | tee -a $out/log && \
samtools view -bt $1 $out/$reads_base.sam > $out/$reads_base.bam && \
echo "sort and index bam started" | tee -a $out/log && \
date | tee -a $out/log && \
samtools sort $out/$reads_base.bam $out/$reads_base.sorted && \
samtools index $out/$reads_base.sorted.bam && \
echo "building vcf started" | tee -a $out/log && \
date | tee -a $out/log && \
samtools faidx $1 && \
samtools mpileup -uf $1 $out/$reads_base.sorted.bam | bcftools view -cvg - > $out/$reads_base.vcf 
echo "done" | tee -a $out/log
date | tee -a $out/log
samtools flagstat $out/$reads_base.sorted.bam | tee -a $out/log) | tee $out/log.full
