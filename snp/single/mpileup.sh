#!/bin/bash
#usage ./mpileup.sh ref.fa aln.bam res_name out_dir bed

ref=$1
bam=$2
name=$3
out=$4
bed=$5

samtools faidx $ref
date
echo -e "0\tStart\t"`date +%s` >> $out/snp_time.log 
if [ -n "$bed" ]; then
    samtools mpileup -q 0 -Q 0 -d 1000 -uf $ref -l $bed $bam | bcftools view -cvg - > $out/$name.vcf
else
    samtools mpileup -q 0 -Q 0 -d 1000 -uf $ref $bam | bcftools view -cvg - > $out/$name.vcf
fi
echo -e "1\tEnd\t"`date +%s` >> $out/snp_time.log 
date
