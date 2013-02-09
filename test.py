#!/usr/bin/python

from argparse import ArgumentParser
from datetime import datetime
import os
import re
import sys
import shutil
import subprocess

dirname = os.path.dirname(sys.argv[0])
base_name_pattern = re.compile("([^/.]+)([^/]+)$")

def cache(directory, storage):
    if (not os.path.exists(storage)):
        os.makedirs(storage)
    stat = os.statvfs(storage)
    free_space = stat.f_bsize * stat.f_bavail / 1024.0 / 1024 / 1024
    if (free_space < 5):
        sys.stderr("Storage has only " + str(free_space) + "Gb left. Not enought space to continue.")
        sys.exit(-1)
    else:
        sys.stderr.write("Moving all sorted bam files, it's indexes and vcf files to external storage, removing " + directory + "\n")
        base_path = storage + "/" + os.path.basename(directory)
        for root, dirs, files in os.walk(directory):
            store = []
            for f in files:
                if (f.find("metrics") != -1 or
                        (len(f) > 4 and f[-4:] == ".log") or
                        (len(f) > 4 and f[-4:] == ".vcf") or
                        (len(f) > 11 and f[-11:] == ".sorted.bam")):
                    store.append(f)
            if (len(store) > 0):
                print root + " " + directory
                rel_path = os.path.relpath(root, directory)
                new_dir = base_path + "/" + rel_path if (root != directory) else base_path
                print new_dir
                if (not os.path.exists(new_dir)):
                    os.makedirs(new_dir)
                for f in store:
                    print f + " " + new_dir + "/" + f
                    shutil.copyfile(root + "/" + f, new_dir + "/" + f)
        shutil.rmtree(directory)

def test_bam(bam, ref, out_dir, is_paired):
    log = open(out_dir + "/metrics_log", "a")
    subprocess.call("java -jar $TOOLS_PATH/picard-tools-1.84/CollectAlignmentSummaryMetrics.jar I=\"" + bam 
            + "\" O=\"" + out_dir + "/metrics_summ\"" 
            + " R=\"" + ref + "\"",
            stdout=log, stderr=log, shell=True)
    subprocess.call("java -jar $TOOLS_PATH/picard-tools-1.84/CollectMultipleMetrics.jar I=\"" + bam 
            + "\" O=\"" + out_dir + "/metrics_mult\"" 
            + " R=\"" + ref + "\"",
            stdout=log, stderr=log, shell=True)
    subprocess.call("java -jar $TOOLS_PATH/picard-tools-1.84/CollectGcBiasMetrics.jar I=\"" + bam 
            + "\" O=\"" + out_dir + "/metrics_gc\"" 
            + " R=\"" + ref 
            + "\" CHART=\"" + out_dir + "/metrics_gc.pdf\"",
            stdout=log, stderr=log, shell=True)
    if (is_paired):
        subprocess.call("java -jar $TOOLS_PATH/picard-tools-1.84/CollectInsertSizeMetrics.jar I=\"" + bam 
                + "\" O=\"" + out_dir + "/metrics_ins\"" 
                + " R=\"" + ref 
                + "\" H=\"" + out_dir + "/metrics_ins.pdf\"",
                stdout=log, stderr=log, shell=True)

def make_bam(directory, ref):
    files = os.listdir(directory)
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
    subprocess.call("samtools sort \"" + directory + "/" + aln 
            + "\" \"" + directory + "/" + aln[0:-3] + "sorted\"",
            stderr=log, stdout=log, shell=True)
    subprocess.call("samtools index \"" + directory + "/" + aln[0:-3] + "sorted.bam\"",
            stderr=log, stdout=log, shell=True)
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
        prog_name = caller[caller.rfind('/') + 1:-3]
        print "Variant calling with " + prog_name + " started"
        out_dir_snp = args.out + "/" + prog_name
        snp_calling(caller, bam, args.ref, out_dir_snp, "")
        print "Variant calling with " + caller + " done"
        sys.exit(0)
else:
    reads1 = args.input if args.input else args.reads1
    reads2 = args.reads2 if is_paired else ""

    for mapper in mappers:
        prog_name = mapper[mapper.rfind('/') + 1:-3]
        print "Mapping with " + prog_name + " started"
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

        print "Mapping with " + prog_name + " done: "

        if (os.path.exists(bam)):
            print "    Succeeded"
            if (not args.skip_test):
                test_bam(bam, args.ref, out_dir, is_paired)

            for caller in snp_callers:
                prog_name = caller[caller.rfind('/') + 1:-3]
                print "Variant calling with " + prog_name + " started"
                out_dir_snp = out_dir + "/" + prog_name
                snp_calling(caller, bam, args.ref, out_dir_snp, "")
                print "Variant calling with " + prog_name + " done"
        else:
            print "    Failed"
        if (args.storage):
            cache(out_dir, args.storage)
