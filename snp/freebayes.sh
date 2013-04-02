#!/bin/bash
#usage ./freebayes.sh ref.fa aln.bam res_name out_dir bed

ref=$1
bam=$2
name=$3
out=$4
bed=$5

date
echo -e "0\tStart\t"`date +%s` >> $out/snp_time.log 
if [ -n "$bed" ]; then
    freebayes -b $bam -f $ref -t $bed -C 10 --min-coverage 10 -m 30 -q 20 -v $out/$name.vcf 
else
    freebayes -b $bam -f $ref -C 10 --min-coverage 10 -m 30 -q 20 -v $out/$name.vcf 
fi
echo -e "1\tEnd\t"`date +%s` >> $out/snp_time.log 
date
