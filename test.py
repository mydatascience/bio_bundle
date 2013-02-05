#!/usr/bin/python

from argparse import ArgumentParser
from datetime import datetime
import os
import re
import sys
import shutil
import subprocess

os.environ['MOSAIK_ANN'] = "/home/kate"
os.environ['STAMPY'] = "/home/kate/bioinf/soft/stampy-1.0.21/stampy.py"
os.environ['SNVER'] = "/home/kate/bioinf/soft/"
os.environ['SHORE_SCORE'] = "/home/kate/bioinf/soft/shore-Linux-x86_64-0.7.1beta/scoring_matrices"
dirname = os.path.dirname(sys.argv[0])
base_name_pattern = re.compile("(\w+)(\.\w+)$")

#TODO doesn't work
def check_space(directory, storage):
    stat = os.statvfs(out_dir)
    free_sp = stat.f_bsize * stat.f_bavail / 1024.0 / 1024 / 1024
    if (free_sp <= 10):
        sys.stderr.write("Only " + str(free_sp) + "Gb of space left")
        if (args.storage):
            stat = os.statvfs(args.storage)
            free_store = stat.f_bsize * stat.f_bavail / 1024.0 / 1024 / 1024
            if (free_store < 5):
                sys.stderr("Storage has only " + str(free_store) + "Gb left")
            else:
                sys.stderr("Moving all sorted bam files, it's indexes and vcf files to external storage, removing " + out_dir)
                while True:
                    return
                        

def test_bam(bam, ref, out_dir, is_paired):
    log = open(out_dir + "/metrics_log", "a")
    subprocess.call("picard-tools CollectAlignmentSummaryMetrics I=\"" + bam 
            + "\" O=\"" + out_dir + "/metrics_summ\"" 
            + "R=\"" + ref + "\"",
            stdout=log, stderr=log, shell=True)
    subprocess.call("picard-tools CollectMultipleMetrics I=\"" + bam 
            + "\" O=\"" + out_dir + "/metrics_mult\"" 
            + "R=\"" + ref + "\"",
            stdout=log, stderr=log, shell=True)
    subprocess.call("picard-tools CalculateHsMetrics I=\"" + bam 
            + "\" O=\"" + out_dir + "/metrics_hs\"" 
            + "R=\"" + ref + "\"",
            stdout=log, stderr=log, shell=True)
    subprocess.call("picard-tools CollectGcBiasMetrics I=\"" + bam 
            + "\" O=\"" + out_dir + "/metrics_gc\"" 
            + "R=\"" + ref + "\"",
            stdout=log, stderr=log, shell=True)
    if (is_paired):
        subprocess.call("picard-tools CollectInsertSizeMetrics I=\"" + bam 
                + "\" O=\"" + out_dir + "/metrics_ins\"" 
                + "R=\"" + ref + "\" H=\"" + out_dir + "/metrics_hist",
                stdout=log, stderr=log, shell=True)

def make_bam(directory, ref):
    files = os.listdir(directory)
    print files
    print directory
    aln = ""
    for f in files:
        if (len(f) <= 4):
            continue
        if (f[-4:] == ".bam"):
            aln = f
            break
        if (f[-4:] == ".sam"):
            aln = f
    if (aln == ""):
        return -1

    log = open(out_dir + "/mapping.log", "a")

    if (aln[-4:] == ".sam"):
        subprocess.call("samtools view -bt \"" + ref 
                + "\" \"" + directory + "/" + aln + "\" "
                + " > \"" + directory + "/" + aln[0:-3] + "bam\"" ,
                stderr=log, stdout=log, shell=True)
        aln = aln[0:-3] + "bam"
    rs = subprocess.call("samtools sort \"" + directory + "/" + aln 
            + "\" \"" + directory + "/" + aln[0:-3] + "sorted\"",
            stderr=log, stdout=log, shell=True)
    print rs
    rs = subprocess.call("samtools index \"" + directory + "/" + aln[0:-3] + "sorted.bam\"",
            stderr=log, stdout=log, shell=True)
    print rs
    return 0

def mapping_paired(mapper, reads1, reads2, ref, out_dir, threads, hashsz, additional):
    if (os.path.exists(out_dir)):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir)

    log = open(out_dir + "/mapping.log", "a")
    subprocess.call(mapper 
            + " \"" + ref + "\" " + base_name_pattern.search(ref).group(1) 
            + " \"" + reads1 + "\" " + base_name_pattern.search(reads1).group(1) 
            + " \"" + reads2 + "\" " + base_name_pattern.search(reads2).group(1) 
            + " \"" + out_dir + "\" " + str(threads) + " " + str(hashsz) + " " + additional,
            stderr=log, stdout=log, shell=True)


def mapping_single(mapper, reads, ref, out_dir, threads, hashsz, additional):
    if (os.path.exists(out_dir)):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir)

    log = open(out_dir + "/mapping.log", "a")
    subprocess.call(mapper 
            + " \"" + ref + "\" " + base_name_pattern.search(ref).group(1) 
            + " \"" + reads + "\" " + base_name_pattern.search(reads).group(1) 
            + " \"" + out_dir + "\" " + str(threads) + " " + str(hashsz) + " " + additional,
            stderr=log, stdout=log, shell=True)
    
def snp_calling(caller, bam, ref, out_dir, additional):
    if (os.path.exists(out_dir)):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir)

    log = open(out_dir + "/snp.log", "a")
    subprocess.call(caller 
            + " \"" + ref + "\""
            + " \"" + bam + "\" " + base_name_pattern.search(bam).group(1) 
            + " \"" + out_dir + "\" " + additional,
            stderr=log, stdout=log, shell=True)

def get_mappers(technology, is_paired, mappers):
    directory = dirname + "/mapping/" + technology + "/" + ("paired" if is_paired else "single") + "/"
    if (not len(mappers)):
        mappers = os.listdir(directory)
    for i in range(len(mappers)):
        mappers[i] = directory + mappers[i] + (".sh" if (mappers[i][-3:] != ".sh") else "")
    return mappers

def get_snp_callers(snp_callers):
    directory = dirname + "/snp/"
    if (not len(snp_callers)):
        snp_callers = os.listdir(directory)
    for i in range(len(snp_callers)):
        snp_callers[i] = directory + snp_callers[i] + (".sh" if (snp_callers[i][-3:] != ".sh") else "")
    return snp_callers

parser = ArgumentParser(description='Test mappers and snp callers')
parser.add_argument("-i", action="store", dest="input", help="fast[a,q] file with reads or sorted bam file for snp calling")
parser.add_argument("-1", action="store", dest="reads1", help="fast[a,q] file with reads")
parser.add_argument("-2", action="store", dest="reads2", help="fast[a,q] file with reads")
parser.add_argument("-r", "--ref", action="store", dest="ref", required=True, help="reference file")
parser.add_argument("-t", "--threads", action="store", dest="threads", default=4, help="threads number, used in multithreaded programms")
parser.add_argument("--hash", action="store", dest="hashsz", default=10, help="hash size, used where this argument exists")
parser.add_argument("-o", "--out", action="store", dest="out", required=True, help="directory to store results")
parser.add_argument("--tech", action="store", dest="technology", default="illumina", choices=['illumina', '454', 'ionTorrent'], help="NGS technology used to get reads")
parser.add_argument("--skip_test", action="store_true", dest="skip_test", help="skip test of sam files")
parser.add_argument("-v", action="store_true", dest="make_vcf", help="make snp calling for mapped files")
parser.add_argument("-m", "--mappers", action="store", dest="mappers", help="list of programms used for mapping")
parser.add_argument("-s", "--snp", action="store", dest="snp_callers", help="list of programms used for snp calling")
parser.add_argument("-c", "--storage", action="store", dest="storage", help="directory to store results for a long lime")

args = parser.parse_args()

if (not args.input and not args.reads1):
    sys.stderr.write("Please, specify reads or alignment file\n")
    sys.exit(-1)

print args

if (dirname == ""):
    dirname = "."
need_map = False
need_snp = False

if (args.input and re.match(".+\.bam", args.input)):
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
    mappers = get_mappers(args.technology, is_paired, mappers)

if (need_snp):
    if (args.snp_callers):
        snp_callers = args.snp_callers.split(" ")
    snp_callers = get_snp_callers(snp_callers)

#print mappers
#print snp_callers

if (not os.path.exists(args.out)):
    os.makedirs(args.out)

if (not need_map):
    for caller in snp_callers:
        snp_calling(caller, args.input, args.ref, args.out, "")
        sys.exit(0)
else:
    reads1 = args.input if args.input else args.reads1
    reads2 = args.reads2 if is_paired else ""

    for mapper in mappers:

        prog_name = mapper[mapper.rfind('/') + 1:-3]
        out_dir = args.out + "/" + prog_name

        if (is_paired):
            mapping_paired(mapper, reads1, reads2, args.ref, out_dir, args.threads, args.hashsz, "")
        else:
            mapping_single(mapper, reads1, args.ref, out_dir, args.threads, args.hashsz, "")
        log = open(out_dir + "/mapping.log", "a")
        log.write("making indexed sorted bam file started: " + str(datetime.now()) + "\n")
        make_bam(out_dir, args.ref)
        log.write("making indexed sorted bam file done: " + str(datetime.now()))
        log.close()
        bam = out_dir + "/" + base_name_pattern.search(reads1).group(1) + ".sorted.bam"
        print bam

        subprocess.call("ls " + out_dir, shell=True)
        print str(os.path.exists(bam))

        if (os.path.exists(bam)):
            if (not args.skip_test):
                test_bam(bam, args.ref, out_dir, is_paired)

            for caller in snp_callers:
                prog_name = caller[caller.rfind('/') + 1:-3]
                out_dir_snp = out_dir + "/" + prog_name
                snp_calling(caller, bam, args.ref, out_dir_snp, "")
