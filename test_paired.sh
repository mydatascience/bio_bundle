#!/bin/bash
#USAGE ./test.sh reference reads1 reads2 result

rm $1.*

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

samtools faidx $1

#bwa
if [ ! -d $4/bwa ]; then
    mkdir $4/bwa
else 
    rm $4/bwa/*
fi

echo "bwa started"
(echo "treads number = $threads; hashsize = $hashsz" | tee $out/bwa/log 
echo "build index started" | tee -a $out/bwa/log
date | tee -a $out/bwa/log
bwa index -a bwtsw $1 && \
echo "align started" | tee -a $out/bwa/log && \
date | tee -a $out/bwa/log && \
bwa aln -t $threads $1 $2 > $out/bwa/$reads1_base.sai && \
bwa aln -t $threads $1 $3 > $out/bwa/$reads2_base.sai && \
echo "making sam started" | tee -a $out/bwa/log && \
date | tee -a $out/bwa/log && \
bwa sampe $1 $out/bwa/$reads1_base.sai $out/bwa/$reads2_base.sai $2 $3 > $out/bwa/$ref_base.sam && \
echo "convert sam to bam started" | tee -a $out/bwa/log && \
date | tee -a $out/bwa/log && \
samtools view -bt $1 $out/bwa/$ref_base.sam > $out/bwa/$ref_base.bam && \
echo "sort and index bam started" | tee -a $out/bwa/log && \
date | tee -a $out/bwa/log && \
samtools sort $out/bwa/$ref_base.bam $out/bwa/$ref_base.sorted && \
samtools index $out/bwa/$reads_base.sorted.bam && \
echo "building vcf started" | tee -a $out/bwa/log && \
date | tee -a $out/bwa/log && \
samtools faidx $1 && \
samtools index $out/bwa/$ref_base.sorted.bam && \
samtools mpileup -uf $1 $out/bwa/$ref_base.sorted.bam | bcftools view -vcg - > $out/bwa/$ref_base.vcf 
echo "done" | tee -a $out/bwa/log
date | tee -a $out/bwa/log
samtools flagstat $out/bwa/$reads_base.sorted.bam | tee -a $out/bwa/log) | tee $out/bwa/log.full
echo "bwa finished"

#bowtie
if [ ! -d "$4/bowtie" ]; then
    mkdir $4/bowtie
else 
    rm $4/bowtie/*
fi
echo "bowtie started"
(echo "treads number = $threads; hashsize = $hashsz" | tee $out/bowtie/log 
echo "build started" | tee -a $out/bowtie/log 
date | tee -a $out/bowtie/log
bowtie-build $1 $out/bowtie/$ref_base.ind && \
echo "mapping started" | tee -a $out/bowtie/log && \
date | tee -a $out/bowtie/log && \
bowtie -t -q -p $threads -a --fr --minins 0 --maxins 500 --sam $out/bowtie/$ref_base.ind -1 $2 -2 $3 > $out/bowtie/$ref_base.sam && \
echo "convert sam to bam started" | tee -a $out/bowtie/log && \
date | tee -a $out/bowtie/log && \
samtools view -bt $1 $out/bowtie/$ref_base.sam > $out/bowtie/$ref_base.bam && \
echo "sort and index bam started" | tee -a $out/bowtie/log && \
date | tee -a $out/bowtie/log && \
samtools sort $out/bowtie/$ref_base.bam $out/bowtie/$ref_base.sorted && \
samtools index $out/bowtie/$ref_base.sorted.bam && \
echo "building vcf started" | tee -a $out/bowtie/log && \
date | tee -a $out/bowtie/log && \
samtools faidx $1 && \
samtools mpileup -uf $1 $out/bowtie/$ref_base.sorted.bam | bcftools view -cvg - > $out/bowtie/$ref_base.vcf 
echo "done" | tee -a $out/bowtie/log
date | tee -a $out/bowtie/log
samtools flagstat $out/bowtie/$ref_base.sorted.bam | tee -a $out/bowtie/log) | tee $out/bowtie/log.full
echo "bowtie finished"

#bowtie2
if [ ! -d "$4/bowtie2" ]; then
    mkdir $4/bowtie2
else 
    rm $4/bowtie2/*
fi
echo "bowtie2 started"
(echo "treads number = $threads; hashsize = $hashsz" | tee $out/bowtie2/log 
echo "build started" | tee -a $out/bowtie2/log 
date | tee -a $out/bowtie2/log
bowtie2-build $1 $out/bowtie2/$ref_base.ind && \
echo "mapping started" | tee -a $out/bowtie2/log && \
date | tee -a $out/bowtie2/log && \
bowtie2 -t -q -p $threads -a --fr --minins 0 --maxins 500 $out/bowtie2/$ref_base.ind -1 $2 -2 $3 > $out/bowtie2/$ref_base.sam && \
echo "convert sam to bam started" | tee -a $out/bowtie2/log && \
date | tee -a $out/bowtie2/log && \
samtools view -bt $1 $out/bowtie2/$ref_base.sam > $out/bowtie2/$ref_base.bam && \
echo "sort and index bam started" | tee -a $out/bowtie2/log && \
date | tee -a $out/bowtie2/log && \
samtools sort $out/bowtie2/$ref_base.bam $out/bowtie2/$ref_base.sorted && \
samtools index $out/bowtie2/$ref_base.sorted.bam && \
echo "building vcf started" | tee -a $out/bowtie2/log && \
date | tee -a $out/bowtie2/log && \
samtools faidx $1 && \
samtools mpileup -uf $1 $out/bowtie2/$ref_base.sorted.bam | bcftools view -cvg - > $out/bowtie2/$ref_base.vcf 
echo "done" | tee -a $out/bowtie2/log
date | tee -a $out/bowtie2/log
samtools flagstat $out/bowtie2/$ref_base.sorted.bam | tee -a $out/bowtie2/log)
echo "bowtie2 finished"

#mosaik
if [ ! -d "$4/mosaik" ]; then
    mkdir $4/mosaik
else 
    rm $4/mosaik/*
fi
echo "mosaik started"
(echo "treads number = $threads; hashsize = $hashsz" | tee $out/mosaik/log 
echo "build started" | tee -a $out/mosaik/log
date | tee -a $out/mosaik/log
MosaikBuild -fr $1 -oa $out/mosaik/$ref_base.dat && \
MosaikBuild -q $2 -out $out/mosaik/$reads1_base.dat -st Illumina && \
MosaikBuild -q $3 -out $out/mosaik/$reads2_base.dat -st Illumina && \
echo "align started" | tee -a $out/mosaik/log && \
date | tee -a $out/mosaik/log && \
MosaikAligner -p $threads -hs $hashsz -in $out/mosaik/$reads1_base.dat -in $out/mosaik/$reads2_base.dat -out $out/mosaik/$ref_base -ia $out/$ref_base.dat -annpe ~/pe.ann -annse ~/se.ann && \
echo "sort and index bam started" | tee -a $out/mosaik/log && \
date | tee -a $out/mosaik/log && \
samtools sort $out/mosaik/$ref_base.bam $out/mosaik/$ref_base.sorted
samtools index $out/mosaik/$ref_base.sorted.bam
echo "building vcf started" | tee -a $out/mosaik/log && \
date | tee -a $out/mosaik/log && \
samtools faidx $1 
samtools mpileup -uf $1 $out/mosaik/$ref_base.sorted.bam | bcftools view -vcg - > $out/mosaik/$ref_base.vcf
echo "done" | tee -a $out.log
date | tee -a $out/mosaik/log
samtools flagstat $out/mosaik/$ref_base.sorted.bam | tee -a $out/log) | tee $out/log.full
echo "mosaik finished"
