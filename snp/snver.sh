#!/bin/bash
#usage freebayes.sh ref.fa aln.bam res_name out_dir additional

ref=$1
bam=$2
name=$3
out=$4
additional=$5

date
java -jar $SNVER/SNVerIndividual.jar -i $bam -r $ref -o $out/$name
date
