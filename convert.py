#!/usr/bin/python
# convert snver and varscan vcf to the format supported by IGV
# ./convert.py from.vcf res.vcf

import sys

if (len(sys.argv) == 3):
    target = open(sys.argv[1])
    res = open(sys.argv[2], "w")
else:
    target = sys.stdin
    res = sys.stdout

for l in target:
    if (l[0] == "#"):
        res.write(l)
    else:
        fields = l.strip().split('\t')
        if (fields[4][0] == "+"):
            fields[4] = fields[3] + fields[4][1:]
            fields[7] = "INDEL;" + fields[7]
        elif (fields[4][0] == "-"):
            fields[3] = fields[3] + fields[4][1:]
            fields[4] = fields[3][:-(len(fields[4]) - 1)]
            fields[7] = "INDEL;" + fields[7]
        for f in fields[:-1]:
            res.write(f + '\t')
        res.write(fields[-1] + '\n')
