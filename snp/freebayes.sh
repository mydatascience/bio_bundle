#!/bin/bash
#usage ./freebayes.sh ref.fa aln.bam res_name out_dir additional

ref=$1
bam=$2
name=$3
out=$4
additional=$5

date
echo -e "0\tStart\t"`date +%s` >> $out/snp_time.log 
freebayes $additional -b $bam -f $ref -v $out/$name.vcf 
echo -e "1\tEnd\t"`date +%s` >> $out/snp_time.log 
date
