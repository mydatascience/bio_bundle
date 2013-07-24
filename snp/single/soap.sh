#!/bin/bash
#usage ./soap.sh ref.fa aln.bam res_name out_dir bed

ref=$1
bam=$2
name=$3
out=$4
bed=$5

date
echo -e "0\tStart\t"`date +%s` >> $out/snp_time.log 
if [ -n "$bed" ]; then
    soapsnp -B $bam -q -m -b 0.1 -P 4 -d $ref -o $out/$name.soap -L 300 -T $bed
else
    soapsnp -B $bam -q -m -b 0.1 -P 4 -d $ref -o $out/$name.soap -L 300 
fi
echo -e "1\tEnd\t"`date +%s` >> $out/snp_time.log 
date

dir=$(dirname $(readlink -e $0))
gunzip $out/$name.soap.gz
python $dir/../../misc/soapsnpToVCF.py $out/$name.soap $bam > $out/$name.vcf
date
