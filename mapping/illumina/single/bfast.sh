#!/bin/bash
#USAGE ./bowtie.sh reference ref_base_name reads reads_base_name result threads hashsz additional

ref=$1
ref_base=$2
reads=$3
reads_base=$4
out=$5
threads=$6
hashsz=$7
additional=$8

echo "build index started" 
date
echo -e "0\tBuilding index\t"`date +%s` >> $out/mapping_time.log
bfast fasta2brg -f $ref && \
bfast index -n $threads -t -f $ref -m 1111111111111111111111 -w $hashsz -i 1 && \
bfast index -n $threads -t -f $ref -m 1111101110111010100101011011111 -w $hashsz -i 2 && \
bfast index -n $threads -t -f $ref -m 1011110101101001011000011010001111111 -w $hashsz -i 3 && \
bfast index -n $threads -t -f $ref -m 10111001101001100100111101010001011111 -w $hashsz -i 4 && \
bfast index -n $threads -t -f $ref -m 11111011011101111011111111 -w $hashsz -i 5 && \
bfast index -n $threads -t -f $ref -m 111111100101001000101111101110111 -w $hashsz -i 6 && \
bfast index -n $threads -t -f $ref -m 11110101110010100010101101010111111 -w $hashsz -i 7 && \
bfast index -n $threads -t -f $ref -m 111101101011011001100000101101001011101 -w $hashsz -i 8 && \
bfast index -n $threads -t -f $ref -m 1111011010001000110101100101100110100111 -w $hashsz -i 9 && \
bfast index -n $threads -t -f $ref -m 1111010010110110101110010110111011 -w $hashsz -i 10
echo "match started" 
date 
echo -e "1\tMaking alignment\t"`date +%s` >> $out/mapping_time.log && \
bfast match -n $threads -f $ref -r $reads > $out/$reads_base.matches.s && \
echo "localalign started" && \
date && \
bfast localalign -n $threads -f $ref -m $out/$reads_base.matches.s > $out/$reads_base.aligned.s && \
echo "postprocess started" && \
date && \
bfast postprocess -n $threads -f $ref -O 1 -i $out/$reads_base.aligned.s > $out/$reads_base.sam && \
echo "postprocess done" && \
date && \
echo -e "2\tAlignment done\t"`date +%s` >> $out/mapping_time.log
