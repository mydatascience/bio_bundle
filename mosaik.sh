#!/bin/bash
#USAGE ./mosaik.sh reference reads result threads hashsz
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
MosaikBuild -fr $1 -oa $out/$ref_base.dat && \
MosaikBuild -q $2 -out $out/$reads_base.dat -st Illumina && \
echo "align started" | tee -a $out/log && \
date | tee -a $out/log && \
MosaikAligner -p $threads -hs $hashsz -in $out/$reads_base.dat -out $out/$reads_base -ia $out/$ref_base.dat -annpe ~/pe.ann -annse ~/se.ann && \
echo "sort and index bam started" | tee -a $out/log && \
date | tee -a $out/log && \
samtools sort $out/$reads_base.bam $out/$reads_base.sorted && \
samtools index $out/$reads_base.sorted.bam && \
echo "building vcf started" | tee -a $out/log && \
date | tee -a $out/log && \
samtools faidx $1 && \
samtools mpileup -uf $1 $out/$reads_base.sorted.bam | bcftools view -vcg - > $out/$reads_base.vcf && \
echo "done" | tee -a $out/log
date | tee -a $out/log
samtools flagstat $out/$reads_base.sorted.bam | tee -a $out/log) | tee $out/log.full
