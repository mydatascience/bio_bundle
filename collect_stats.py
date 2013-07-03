#!/usr/bin/python
# -*- coding: utf-8 -*-

from argparse import ArgumentParser
from lib import statslib
from lib import stats_writer
import os.path
import os
import subprocess
import shutil
import sys

import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    parser = ArgumentParser(description="Collect statistics from test.py")
    parser.add_argument("-i", action="store", dest="in_dir", help="directory with test.py script results", required=True)
    parser.add_argument("--vcf", action="store", dest="vcf", help="golden file with variants")
    parser.add_argument("-n", "--name", action="store", dest="name", help="base name for vcf files", required=True)

    args = parser.parse_args()

    print args

    golden_vcf = statslib.VarCaller("golden")

    if (os.path.exists(args.in_dir + "/test_stats")):
        shutil.rmtree(args.in_dir + "/test_stats")
    os.makedirs(args.in_dir + "/test_stats")
    os.makedirs(args.in_dir + "/test_stats/plots")

    if args.vcf != None and not os.path.exists(args.vcf):
        print "Golden file not exists"
        sys.exit(-1)
    elif args.vcf != None:
        golden_vcf.process_vcf("golden", "golden", args.in_dir + "/test_stats", args.vcf, None, os.path.dirname(sys.argv[0]))

    aligners = [] # [Aligner1 Aligner2 ...]
    var_callers = {} # {aligner_name->[var_caller1 ...] ... }
    var_caller_names = []
    aln_dirs = os.listdir(args.in_dir)

    for aln_name in aln_dirs:
        curr_aln_dir = args.in_dir + "/" + aln_name
        if (not os.path.isdir(curr_aln_dir)):
            continue

        aln_files = os.listdir(curr_aln_dir)

        curr_aligner = statslib.Aligner(aln_name)
        curr_aligner.fill_metrics(curr_aln_dir)

        if (not curr_aligner.empty()):
            aligners.append(curr_aligner)

        for var_caller_name in aln_files:
            curr_var_dir = curr_aln_dir + "/" + var_caller_name
            if (not os.path.isdir(curr_var_dir)):
                continue

            curr_var_caller = statslib.VarCaller(var_caller_name)
            curr_var_caller.fill_stats(curr_var_dir, aln_name, var_caller_name, args.in_dir, args.name, args.vcf, os.path.dirname(sys.argv[0]))

            if (not curr_var_caller.empty()):
                if curr_var_caller.name not in var_caller_names:
                    var_caller_names.append(curr_var_caller.name)
                if (aln_name not in var_callers):
                    var_callers[aln_name] = [curr_var_caller]
                else:
                    var_callers[aln_name].append(curr_var_caller)

    var_caller_names = sorted(var_caller_names)

    aligners = sorted(aligners, key=lambda aln: aln.name)
    stats_writer.remove_empty_metrics(aligners)

#for aln in aligners:
#    print str(aln)
#
#for aln in var_callers.keys():
#    print aln + ":"
#    for var_c in var_callers[aln]:
#        print "\t" + str(var_c)

# Output statistics

    tex_f = open(args.in_dir + "/test_stats/stats.tex", "w")
    stats_writer.begin_document(tex_f)


# Aligners stats
    if len(aligners):
        stats_writer.new_section(tex_f, 'Aligner results')
        stats_writer.output_aligners_time(aligners, tex_f, args.in_dir)
        stats_writer.output_aligners_stats(aligners, tex_f, args.in_dir)

# variant callers stats
    stats_writer.new_section(tex_f, 'Variant callers results')
    stats_writer.vc_single_stats(var_callers, tex_f, args.in_dir, var_caller_names)

    if args.vcf != None:
        stats_writer.new_section(tex_f, 'Variant caller shared statstics')
        stats_writer.vc_shared_stats(var_callers, tex_f, args.in_dir, var_caller_names, golden_vcf)

    stats_writer.end_document(tex_f)
