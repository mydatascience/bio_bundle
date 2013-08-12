import os
import shutil
import sys
import subprocess
import time
from datetime import datetime

def flushres(directory, storage):
    if (not os.path.exists(storage)):
        os.makedirs(storage)
    stat = os.statvfs(storage)
    free_space = stat.f_bsize * stat.f_bavail / 1024.0 / 1024 / 1024
    if (free_space < 1):
        sys.write.stderr("Storage has only " + str(free_space) + "Gb left. Not enought space to continue.")
        sys.exit(-1)
    else:
        sys.stderr.write("Moving all bam files, it's indexes and vcf files to external storage, removing " + directory + "\n")
        base_path = storage + "/" + os.path.basename(directory)
        for root, dirs, files in os.walk(directory):
            store = []
            for f in files:
                if (f.find("metrics") != -1 or f.endswith(".log") or 
                        f.endswith(".vcf") or f.endswith(".bam")):
                    store.append(f)
            if (len(store) > 0):
                rel_path = os.path.relpath(root, directory)
                new_dir = base_path + "/" + rel_path if (root != directory) else base_path
                if (not os.path.exists(new_dir)):
                    os.makedirs(new_dir)
                for f in store:
                    shutil.copyfile(root + "/" + f, new_dir + "/" + f)
        shutil.rmtree(directory)

def test_bam(bam, ref, out_dir, is_paired):
    log = open(out_dir + "/metrics_log", "a")
    subprocess.call("java -jar $TOOLS_PATH/picard-tools-1.84/CollectMultipleMetrics.jar I=\"" + bam 
            + "\" O=\"" + out_dir + "/metrics\" ASSUME_SORTED=true" 
            + " R=\"" + ref + "\" VALIDATION_STRINGENCY=SILENT",
            stdout=log, stderr=log, shell=True)

def create_repeats_reg(ref):
    reg = open(ref + '.rep.bed', 'w')
    subprocess.call("$TOOLS_PATH/RepeatMasker/RepeatMasker -lib $TOOLS_PATH/RepeatMasker/Libraries/RepBase17.07.fa -pa 4 " + ref, shell=True)
    subprocess.call('$TOOLS_PATH/RepeatMasker/Libraries/rmOutToBED.pl ' + ref + '.out 0', stdout = reg, shell=True)

def make_bam(directory, ref):
    files = os.listdir(directory)
    aln = ""
    for f in files:
        if (len(f) <= 4):
            continue
        if (f.endswith(".bam")):
            aln = f
            break
        if (f.endswith(".sam")):
            aln = f
    if (aln == ""):
        return ''

    log = open(directory + "/mapping.log", "a")

    aln_base = os.path.splitext(aln)[0]
    if (aln.endswith(".sam")):
        subprocess.call("samtools view -bt \"" + ref + '.fai'
                + "\" \"" + directory + "/" + aln + "\" "
                + " > \"" + directory + "/" + aln_base + ".bam\"" ,
                stderr=log, stdout=log, shell=True)
        os.remove(directory + "/" + aln)
        aln = aln_base + ".bam"
    subprocess.call("samtools sort \"" + directory + "/" + aln 
            + "\" \"" + directory + "/" + aln_base + ".sorted\"",
            stderr=log, stdout=log, shell=True)
    os.remove(directory + "/" + aln)
    os.rename(directory + '/' + aln_base + '.sorted.bam', directory + '/' + aln_base + '.bam')
    subprocess.call("samtools index \"" + directory + "/" + aln_base + ".bam\"",
            stderr=log, stdout=log, shell=True)

    time_log = open(directory + "/mapping_time.log", "a")
    time_log.write("3\tSorted bam file done\t" + str(int(time.mktime(datetime.now().timetuple()))))
    return directory + '/' + aln_base + '.bam'

def mapping_paired(mapper, reads1, reads2, ref, out_dir, threads, hashsz):
    if (os.path.exists(out_dir)):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir)

    log = open(out_dir + "/mapping.log", "a")
    subprocess.call(mapper 
            + " \"" + ref + "\" " + os.path.splitext(os.path.basename(ref))[0]
            + " \"" + reads1 + "\" " + os.path.splitext(os.path.basename(reads1))[0] 
            + " \"" + reads2 + "\" " + os.path.splitext(os.path.basename(reads2))[0] 
            + " \"" + out_dir + "\" " + str(threads) + " " + str(hashsz),
            stderr=log, stdout=log, shell=True)


def mapping_single(mapper, reads, ref, out_dir, threads, hashsz):
    if (os.path.exists(out_dir)):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir)

    log = open(out_dir + "/mapping.log", "a")
    subprocess.call(mapper 
            + " \"" + ref + "\" " + os.path.splitext(os.path.basename(ref))[0] 
            + " \"" + reads + "\" " + os.path.splitext(os.path.basename(reads))[0]
            + " \"" + out_dir + "\" " + str(threads) + " " + str(hashsz),
            stderr=log, stdout=log, shell=True)
    
def snp_calling(caller, bam, ref, out_dir, bed):
    if (os.path.exists(out_dir)):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir)

    log = open(out_dir + "/snp.log", "a")
    subprocess.call(caller 
            + " \"" + ref + "\""
            + " \"" + bam + "\" " + os.path.splitext(os.path.basename(bam))[0]
            + " \"" + out_dir + "\" " + (bed if bed else ""),
            stderr=log, stdout=log, shell=True)
    log.close()

def make_seq_dict(ref):
    subprocess.call("java -jar $TOOLS_PATH/picard-tools-1.84/CreateSequenceDictionary.jar R=\"" + ref
            + "\" O=\"" + os.path.splitext(ref)[0] + ".dict\"", shell = True)

def get_mappers(technology, is_paired, mappers, dirname):
    directory = dirname + "/mapping/" + technology + "/" + ("paired" if is_paired else "single") + "/"
    if (not len(mappers)):
        mappers = os.listdir(directory)
    for i in range(len(mappers)):
        mappers[i] = directory + mappers[i] + (".sh" if (mappers[i][-3:] != ".sh") else "")
    return mappers

def get_snp_callers(snp_callers, dirname, is_mult_bam):
    directory = dirname + "/snp/" + ('multiple/' if is_mult_bam else 'single/')
    if (not len(snp_callers)):
        snp_callers = os.listdir(directory)
    for i in range(len(snp_callers)):
        snp_callers[i] = directory + snp_callers[i] + (".sh" if (snp_callers[i][-3:] != ".sh") else "")
    return snp_callers

def filter_vcf(vcf, ref):
    basename = os.path.splitext(vcf)[0]
    dirname = os.path.dirname(vcf) if os.path.dirname(vcf) else '.'
    log = open(dirname + "/snp.log", "a")
#    if not os.path.exists(ref + '.rep.bed'):
#        create_repeats_reg(ref)

    subprocess.call("vcftools --minQ 30 --hwe 0.05 "
            + " --non-ref-ac 2 --minDP 5 --recode --recode-INFO-all --out " + basename + " --vcf " + vcf,
#            + " --exclude-bed " + ref + '.rep.bed', 
            stdout=log, stderr=log, shell=True)
    if not os.path.exists(basename + '.recode.vcf'):
        subprocess.call("vcftools --hwe 0.05 "
                + " --non-ref-ac 2 --minDP 5 --recode --recode-INFO-all --out " + basename + " --vcf " + vcf,
#                + " --exclude-bed " + ref + '.rep.bed', 
                stdout=log, stderr=log, shell=True)
    os.rename(basename + '.vcf', basename + '.raw.vcf')
    if os.path.exists(basename + '.recode.vcf'):
        os.rename(basename + '.recode.vcf', basename + '.vcf')

def make_snp_calling(caller, ref, bam, bed, tmp_dir):
    prog_name = os.path.splitext(os.path.basename(caller))[0]
    print "Variant calling with " + prog_name + " started"
    out_dir_snp = tmp_dir + "/" + prog_name
    snp_calling(caller, bam, ref, out_dir_snp, bed)
    vcf_name = out_dir_snp + '/' + os.path.splitext(os.path.basename(bam))[0] + '.vcf' 
    if os.path.exists(vcf_name):
        print "Filter vcf file"
        filter_vcf(vcf_name, ref)
    print "Variant calling with " + prog_name + " done"
    return out_dir_snp

def make_mapping(mapper, is_paired, reads1, reads2, ref, threads, hashsz, tmp_dir):
    prog_name = os.path.splitext(os.path.basename(mapper))[0]
    print "Mapping with " + prog_name + " started"
    out_dir = tmp_dir + "/" + prog_name

    if (is_paired):
        mapping_paired(mapper, reads1, reads2, ref, out_dir, threads, hashsz)
    else:
        mapping_single(mapper, reads1, ref, out_dir, threads, hashsz)
    log = open(out_dir + "/mapping.log", "a")
    log.write("making indexed sorted bam file started: " + str(datetime.now()) + "\n")
    bam = make_bam(out_dir, ref)
    log.write("making indexed sorted bam file done: " + str(datetime.now()))
    log.close()
    return (out_dir, bam)

    print "Mapping with " + prog_name + " done: "
