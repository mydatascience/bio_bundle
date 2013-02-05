#!/bin/bash
#USAGE ./mosaik.sh reference reads1 read2 result threads hashsz
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
MosaikBuild -fr $1 -oa $out/$ref_base.dat && \
MosaikBuild -q $2 -out $out/$reads1_base.dat -st Illumina && \
MosaikBuild -q $3 -out $out/$reads2_base.dat -st Illumina && \
echo "align started" | tee -a $out/log && \
date | tee -a $out/log && \
MosaikAligner -p $threads -hs $hashsz -in $out/$reads1_base.dat -in $out/$reads2_base.dat -out $out/$ref_base -ia $out/$ref_base.dat -annpe ~/pe.ann -annse ~/se.ann && \
echo "sort and index bam started" | tee -a $out/log && \
date | tee -a $out/log && \
samtools sort $out/$ref_base.bam $out/$ref_base.sorted
samtools index $out/$ref_base.sorted.bam
echo "building vcf started" | tee -a $out/log && \
date | tee -a $out/log && \
samtools faidx $1 
samtools mpileup -uf $1 $out/$ref_base.sorted.bam | bcftools view -vcg - > $out/$ref_base.vcf
echo "done" | tee -a $out/log
date | tee -a $out/log
samtools flagstat $out/$ref_base.sorted.bam | tee -a $out/log) | tee $out/log.full
