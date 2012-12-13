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
ref_base=`echo $1 | perl -pe "s/^.*?([^\/]+?)(\.[^\/]+)?$/\1/; s/\n//" `    # получение базового имени от пути к файлу референса
reads_base=`echo $2 | perl -pe "s/^.*?([^\/]+?)(\.[^\/]+)?$/\1/; s/\n//" `  # получение базового имени от пути к файлу с ридами
out=$3

rm $1.*

samtools faidx $1 # индексирование референса

#bwa
if [ ! -d $out/bwa ]; then
    mkdir $out/bwa
else 
    rm $out/bwa/*
fi

echo "bwa started"
(echo "treads number = $threads; hashsize = $hashsz" | tee $out/bwa/log 
echo "build index started" | tee -a $out/bwa/log
date | tee -a $out/bwa/log
bwa index -a bwtsw $1 && \ # индексирование генома по алгоритму burrows-wheeler
date | tee -a $out/bwa/log && \
bwa aln -t $threads $1 $2 > $out/bwa/$reads_base.sai && \ # выравнивание ридов на референс
echo "making sam started" | tee -a $out/bwa/log && \
date | tee -a $out/bwa/log && \
bwa samse $1 $out/bwa/$reads_base.sai $2 > $out/bwa/$reads_base.sam && \ # преобразование выравнивания из внутреннего формата bwa в sam
echo "convert sam to bam started" | tee -a $out/bwa/log && \
date | tee -a $out/bwa/log && \
samtools view -bt $1 $out/bwa/$reads_base.sam > $out/bwa/$reads_base.bam && \ | # преобразование sam в bam
echo "sort and index bam started" | tee -a $out/bwa/log && \
date | tee -a $out/bwa/log && \
samtools sort $out/bwa/$reads_base.bam $out/bwa/$reads_base.sorted && \ # сортировка bam файла
samtools index $out/bwa/$reads_base.sorted.bam && \ # индексирование получившегося bam файла
echo "building vcf started" | tee -a $out/bwa/log && \
date | tee -a $out/bwa/log && \
samtools mpileup -uf $1 $out/bwa/$reads_base.sorted.bam | bcftools view -vcg - > $out/bwa/$reads_base.vcf # получение vcf файла из bam файла
echo "done" | tee -a $out/bwa/log
date | tee -a $out/bwa/log
samtools flagstat $out/bwa/$reads_base.sorted.bam | tee -a $out/bwa/log) | tee $out/bwa/log.full
echo "bwa finished"

#bowtie
if [ ! -d "$out/bowtie" ]; then
    mkdir $out/bowtie
else 
    rm $out/bowtie/*
fi
echo "bowtie started"
(echo "treads number = $threads; hashsize = $hashsz" | tee $out/bowtie/log 
echo "build started" | tee -a $out/bowtie/log 
date | tee -a $out/bowtie/log
bowtie-build $1 $out/bowtie/$ref_base.ind && \ # построение индекса для референса
echo "mapping started" | tee -a $out/bowtie/log && \
date | tee -a $out/bowtie/log && \
bowtie -t -q -p $threads --best --sam $out/bowtie/$ref_base.ind $2 > $out/bowtie/$reads_base.sam && \ # выравнивание ридов на референс
echo "convert sam to bam started" | tee -a $out/bowtie/log && \
date | tee -a $out/bowtie/log && \
samtools view -bt $1 $out/bowtie/$reads_base.sam > $out/bowtie/$reads_base.bam && \ # преобразование sam в bam
echo "sort and index bam started" | tee -a $out/bowtie/log && \
date | tee -a $out/bowtie/log && \
samtools sort $out/bowtie/$reads_base.bam $out/bowtie/$reads_base.sorted && \ # сортировка bam файла
samtools index $out/bowtie/$reads_base.sorted.bam && \ # индексирование получившегося bam файла
echo "building vcf started" | tee -a $out/bowtie/log && \
date | tee -a $out/bowtie/log && \
samtools mpileup -uf $1 $out/bowtie/$reads_base.sorted.bam | bcftools view -cvg - > $out/bowtie/$reads_base.vcf # получение vcf файла из bam файла
echo "done" | tee -a $out/bowtie/log
date | tee -a $out/bowtie/log
samtools flagstat $out/bowtie/$reads_base.sorted.bam | tee -a $out/bowtie/log) | tee $out/bowtie/log.full
echo "bowtie finished"

#bowtie2
if [ ! -d "$out/bowtie2" ]; then
    mkdir $out/bowtie2
else 
    rm $out/bowtie2/*
fi
echo "bowtie2 started"
(echo "treads number = $threads; hashsize = $hashsz" | tee $out/bowtie2/log 
echo "build started" | tee -a $out/bowtie2/log 
date | tee -a $out/bowtie2/log
bowtie2-build $1 $out/bowtie2/$ref_base.ind && \ # построение индекса для референса
echo "mapping started" | tee -a $out/bowtie2/log && \
date | tee -a $out/bowtie2/log && \
bowtie2 -t -q -p $threads --sensitive --sam-rg -x $out/bowtie2/$ref_base.ind -U $2 > $out/bowtie2/$reads_base.sam && \ # выравнивание ридов на рееренс
echo "convert sam to bam started" | tee -a $out/bowtie2/log && \
date | tee -a $out/bowtie2/log && \
samtools view -bt $1 $out/bowtie2/$reads_base.sam > $out/bowtie2/$reads_base.bam && \ # преобразование sam в bam
echo "sort and index bam started" | tee -a $out/bowtie2/log && \
date | tee -a $out/bowtie2/log && \
samtools sort $out/bowtie2/$reads_base.bam $out/bowtie2/$reads_base.sorted && \ # сортировка bam файла
samtools index $out/bowtie2/$reads_base.sorted.bam && \ # индексирование получившегося bam файла
echo "building vcf started" | tee -a $out/bowtie2/log && \
date | tee -a $out/bowtie2/log && \
samtools mpileup -uf $1 $out/bowtie2/$reads_base.sorted.bam | bcftools view -cvg - > $out/bowtie2/$reads_base.vcf # получение vcf файла из bam файла
echo "done" | tee -a $out/bowtie2/log
date | tee -a $out/bowtie2/log
samtools flagstat $out/bowtie2/$reads_base.sorted.bam | tee -a $out/bowtie2/log) | tee $out/bowtie2/log.full
echo "bowtie2 finished"

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
MosaikBuild -fr $1 -oa $out/mosaik/$ref_base.dat && \ # преобразование файла референса во внутренний формат mosaik
MosaikBuild -q $2 -out $out/mosaik/$reads_base.dat -st Illumina && \ # преобразование файла ридов во внутренний формат mosaik
echo "align started" | tee -a $out/mosaik/log && \
date | tee -a $out/mosaik/log && \
MosaikAligner -p $threads -hs $hashsz -in $out/mosaik/$reads_base.dat -out $out/mosaik/$reads_base -ia $out/mosaik/$ref_base.dat -annpe ~/pe.ann -annse ~/se.ann && \ # выравнивание ридов
echo "sort and index bam started" | tee -a $out/mosaik/log && \
date | tee -a $out/mosaik/log && \
samtools sort $out/mosaik/$reads_base.bam $out/mosaik/$reads_base.sorted && \ # сортировка bam файла
samtools index $out/mosaik/$reads_base.sorted.bam && \ # индексирование получившегося bam файла
echo "building vcf started" | tee -a $out/mosaik/log && \
date | tee -a $out/mosaik/log && \
samtools mpileup -uf $1 $out/mosaik/$reads_base.sorted.bam | bcftools view -vcg - > $out/mosaik/$reads_base.vcf && \ # получение vcf файла из bam файла
echo "done" | tee -a $out/mosaik/log
date | tee -a $out/mosaik/log
samtools flagstat $out/mosaik/$reads_base.sorted.bam | tee -a $out/mosaik/log) | tee $out/mosaik/log.full
echo "mosaik finished"
