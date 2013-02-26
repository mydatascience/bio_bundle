#!/bin/bash
#usage ./mpileup.sh ref.fa aln.bam res_name out_dir additional

ref=$1
bam=$2
name=$3
out=$4
additional=$5

samtools faidx $ref
date
echo -e "0\tStart\t"`date +%s` >> $out/snp_time.log 
samtools mpileup -uf $ref $bam | bcftools view -cvg - > $out/$name.vcf
echo -e "1\tEnd\t"`date +%s` >> $out/snp_time.log 
date
