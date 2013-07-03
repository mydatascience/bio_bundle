#!/usr/bin/python

from argparse import ArgumentParser
from datetime import datetime
import os
import re
import sys
import shutil
from lib import testlib

if __name__ == "__main__":
    dirname = os.path.dirname(sys.argv[0])
    parser = ArgumentParser(description='Test mappers and snp callers')
    parser.add_argument("-i", action="store", dest="input", help="fast[a,q] file with reads or sorted bam file for snp calling")
    parser.add_argument("-1", action="store", dest="reads1", help="fast[a,q] file with reads")
    parser.add_argument("-2", action="store", dest="reads2", help="fast[a,q] file with reads")
    parser.add_argument("-r", "--ref", action="store", dest="ref", required=True, help="reference file")
    parser.add_argument("-t", "--threads", action="store", dest="threads", default=4, help="threads number, used in multithreaded programms")
    parser.add_argument("--hash", action="store", dest="hashsz", default=10, help="hash size, used where this argument exists")
    parser.add_argument("--tech", action="store", dest="technology", default="illumina", choices=['illumina', '454', 'ionTorrent'], help="NGS technology used to get reads")
    parser.add_argument("--skip_test", action="store_true", dest="skip_test", help="skip test of sam files")
    parser.add_argument("-v", action="store_true", dest="make_vcf", help="make snp calling for mapped files")
    parser.add_argument("-m", "--mappers", action="store", dest="mappers", help="list of programms used for mapping")
    parser.add_argument("-s", "--snp", action="store", dest="snp_callers", help="list of programms used for snp calling")
    parser.add_argument("-o", "--out", action="store", required=True, dest="storage", help="directory to store results")
    parser.add_argument("--bed", action="store", dest="bed", help="region file in bed format")

    args = parser.parse_args()
    tmp_dir = args.storage + '/tmp'

    if (not args.input and not args.reads1):
        sys.stderr.write("Please, specify reads or alignment file\n")
        sys.exit(-1)

    print args

    if (dirname == ""):
        dirname = "."
    (need_map, need_snp) = (False, False)

    if (args.input and args.input.endswith('.bam')):
        need_snp = True
    else:
        need_map = True
        if (args.make_vcf or args.snp_callers):
            need_snp = True

    mappers = []
    snp_callers = []
    is_paired = True if args.reads2 else False

    if (need_map):
        if (args.mappers):
            mappers = args.mappers.split(" ")
        mappers = testlib.get_mappers(args.technology, is_paired, mappers, dirname)

    if (need_snp):
        if (args.snp_callers):
            snp_callers = args.snp_callers.split(" ")
        snp_callers = testlib.get_snp_callers(snp_callers, dirname)

    if not os.path.exists(args.storage):
        os.makedirs(args.storage)

    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    if (not need_map):
        for caller in snp_callers:
            out_dir_snp = testlib.make_snp_calling(caller, args.ref, args.input, args.bed, tmp_dir)
            testlib.flushres(out_dir_snp, args.storage)
        sys.exit(0)
    else:
        reads1 = args.input if args.input else args.reads1
        reads2 = args.reads2 if is_paired else ""

        for mapper in mappers:
            (out_dir, bam) = testlib.make_mapping(mapper, is_paired, reads1, reads2, args.ref, args.threads, args.hashsz, tmp_dir)
            if (os.path.exists(bam) and os.path.getsize(bam) > 0):
                print "    Succeeded"
                if (not args.skip_test):
                    testlib.test_bam(bam, args.ref, out_dir, is_paired)

                for caller in snp_callers:
                    testlib.make_snp_calling(caller, args.ref, bam, args.bed, out_dir)
            else:
                print "    Failed"
            testlib.flushres(out_dir, args.storage)

    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)
