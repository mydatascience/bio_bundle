#!/bin/bash
#usage freebayes.sh ref.fa aln.bam res_name out_dir bed

ref=$1
bam=$2
name=$3
out=$4
bed=$5

date
echo -e "0\tStart\t"`date +%s` >> $out/snp_time.log 
args=""
while read line
do
    basename=$(basename $line)
    args="$args -V $out/$basename.vcf"
    if [ -n "$bed" ]; then
        java -jar $TOOLS_PATH/SNVerIndividual.jar -b 0.05 -mq 0 -bq 0 -a 1 -i $line -r $ref -l $bed -o $out/$basename
    else
        java -jar $TOOLS_PATH/SNVerIndividual.jar -b 0.05 -mq 0 -bq 0 -a 1 -i $line -r $ref -o $out/$basename
    vcf-concat $out/$basename.filter.vcf $out/$basename.indel.filter.vcf | vcf-sort > $out/$basename.vcf
    fi
    echo $args
done <$bam
echo $args
GenomeAnalysisTK.jar -T CombineVariants $args -o $out/$name.vcf -genotypeMergeOptions UNIQUIFY -R $ref
echo -e "1\tEnd\t"`date +%s` >> $out/snp_time.log && \
date && \
