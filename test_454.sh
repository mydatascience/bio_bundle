#!/bin/bash
#USAGE ./test.sh reference reads result

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

rm $1.*

#bwa
if [ ! -d $out/bwa ]; then
    mkdir $out/bwa
else 
    rm $out/bwa/*
fi

samtools faidx $1

echo "bwa started"
(echo "treads number = $threads; hashsize = $hashsz" | tee $out/bwa/log 
echo "build index started" | tee -a $out/bwa/log
date | tee -a $out/bwa/log
bwa index -a bwtsw $1 && \
echo "align started" | tee -a $out/bwa/log && \
date | tee -a $out/bwa/log && \
bwa bwasw -t $threads $1 $2 > $out/bwa/$reads_base.sam && \
echo "convert sam to bam started" | tee -a $out/bwa/log && \
date | tee -a $out/bwa/log && \
samtools view -bt $1 $out/bwa/$reads_base.sam > $out/bwa/$reads_base.bam && \
echo "sort and index bam started" | tee -a $out/bwa/log && \
date | tee -a $out/bwa/log && \
samtools sort $out/bwa/$reads_base.bam $out/bwa/$reads_base.sorted && \
samtools index $out/bwa/$reads_base.sorted.bam && \
echo "building vcf started" | tee -a $out/bwa/log && \
date | tee -a $out/bwa/log && \
samtools mpileup -uf $1 $out/bwa/$reads_base.sorted.bam | bcftools view -vcg - > $out/bwa/$reads_base.vcf 
echo "done" | tee -a $out/log
date | tee -a $out/log
samtools flagstat $out/bwa/$reads_base.sorted.bam | tee -a $out/bwa/log) | tee $out/bwa/log.full
echo "bwa finished"

#mosaik
if [ ! -d "$out/mosaik" ]; then
    mkdir $out/mosaik
else 
    rm $out/mosaik/*
fi
echo "mosaik started"
(echo "treads number = $threads; hashsize = $hashsz" | tee $out/mosaik/log 
echo "build started" | tee -a $out/mosaik/log
date | tee -a $out/mosaik/log
MosaikBuild -fr $1 -oa $out/mosaik/$ref_base.dat && \
MosaikBuild -q $2 -out $out/mosaik/$reads_base.dat -st 454 && \
echo "align started" | tee -a $out/mosaik/log && \
date | tee -a $out/mosaik/log && \
MosaikAligner -p $threads -hs $hashsz -in $out/mosaik/$reads_base.dat -out $out/mosaik/$reads_base -ia $out/mosaik/$ref_base.dat -annpe ~/pe.ann -annse ~/se.ann && \
echo "sort and index bam started" | tee -a $out/mosaik/log && \
date | tee -a $out/mosaik/log && \
samtools sort $out/mosaik/$reads_base.bam $out/$reads_base.sorted && \
samtools index $out/mosaik/$reads_base.sorted.bam && \
echo "building vcf started" | tee -a $out/mosaik/log && \
date | tee -a $out/mosaik/log && \
samtools faidx $1 && \
samtools mpileup -uf $1 $out/mosaik/$reads_base.sorted.bam | bcftools view -vcg - > $out/mosaik/$reads_base.vcf && \
echo "done" | tee -a $out/log
date | tee -a $out/mosaik/log
samtools flagstat $out/mosaik/$reads_base.sorted.bam | tee -a $out/mosaik/log) | tee $out/mosaik/log.full
echo "mosaik finished"
