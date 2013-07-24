#!/usr/bin/python
import sys
import os.path

#alphabet = ['A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'N': 'N',
#    'R': 'A,G', 'Y': 'C,T', 'K': 'G,T', 'M': 'A,C', 'S': 'C,G', 'W': 'A,T',
#    'B': 'C,G,T', 'D': 'A,G,T', 'H': 'A,C,T', 'V': 'A,C,G']
usage = 'soapsnpToVCF.py in.soap <aln.bam>'
header = r'''##fileformat=VCFv4.1
##samtoolsVersion=0.1.18 (r982:295)
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
##INFO=<ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">
##INFO=<ID=MQ,Number=1,Type=Integer,Description="Root-mean-square mapping quality of covering reads">
##INFO=<ID=FQ,Number=1,Type=Float,Description="Phred probability of all samples being the same">
##INFO=<ID=AF1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele frequency (assuming HWE)">
##INFO=<ID=AC1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele count (no HWE assumption)">
##INFO=<ID=G3,Number=3,Type=Float,Description="ML estimate of genotype frequencies">
##INFO=<ID=HWE,Number=1,Type=Float,Description="Chi^2 based HWE test P-value based on G3">
##INFO=<ID=CLR,Number=1,Type=Integer,Description="Log ratio of genotype likelihoods with and without the constraint">
##INFO=<ID=UGT,Number=1,Type=String,Description="The most probable unconstrained genotype configuration in the trio">
##INFO=<ID=CGT,Number=1,Type=String,Description="The most probable constrained genotype configuration in the trio">
##INFO=<ID=PV4,Number=4,Type=Float,Description="P-values for strand bias, baseQ bias, mapQ bias and tail distance bias">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
##INFO=<ID=PC2,Number=2,Type=Integer,Description="Phred probability of the nonRef allele frequency in group1 samples being larger (,smaller) than in group2.">
##INFO=<ID=PCHI2,Number=1,Type=Float,Description="Posterior weighted chi^2 P-value for testing the association between group1 and group2 samples.">
##INFO=<ID=QCHI2,Number=1,Type=Integer,Description="Phred scaled PCHI2.">
##INFO=<ID=PR,Number=1,Type=Integer,Description="# permutations yielding a smaller PCHI2.">
##INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias">            
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">                       
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">              
##FORMAT=<ID=GL,Number=3,Type=Float,Description="Likelihoods for RR,RA,AA genotypes (R=ref,A=alt)">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="# high-quality bases">          
##FORMAT=<ID=SP,Number=1,Type=Integer,Description="Phred-scaled strand bias P-value">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
''' + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"

if __name__ == '__main__':
    if len(sys.argv) < 3 or not os.path.exists(sys.argv[1]):
        sys.stderr.write(usage + '\n')
        sys.exit(-1)
    f = open(sys.argv[1])

    print header + '\t' + sys.argv[2]
    for l in f:
        fields = l.strip().split('\t')
        if fields[2] == fields[3]:
            continue
        chr_name = fields[0]
        pos = fields[1]
        ref = fields[2]
        qual = fields[4]
        dp = fields[13]
        if fields[5] == fields[2]:
            alt = fields[9]
            dp4 = fields[8] + ',0,' + fields[12] + ',0'
        else:
            alt = fields[5]
            dp4 = fields[12] + ',0,' + fields[8] + ',0'
        if alt == fields[3]:
            gt = '1/1'
        else:
            gt = '1/0'
        print '\t'.join([chr_name, pos, '.', ref, alt, qual, '.', 'DP=' + dp + ';DP4=' + dp4, 'GT', gt]) 

