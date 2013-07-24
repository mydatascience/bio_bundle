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
    freebayes -b $bam -f $ref -t $bed -C 1 -R 0 -S 0 -F 0.05 --min-coverage 1 -m 0 -q 0 -v $out/$name.vcf 
else
    freebayes -b $bam -f $ref -C 1 -R 0 -S 0 -F 0.05 --min-coverage 1 -m 0 -q 0 -v $out/$name.vcf 
fi
echo -e "1\tEnd\t"`date +%s` >> $out/snp_time.log 
date
