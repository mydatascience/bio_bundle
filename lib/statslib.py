#!/usr/bin/python
# -*- coding: utf-8 -*-

from datetime import datetime
from datetime import timedelta
import os.path
import os
import subprocess
import shutil
import sys

import numpy as np
import matplotlib.pyplot as plt

class Aligner:

    def __init__(self, name, 
            build_ind_time = timedelta(0), 
            aln_time = timedelta(0), 
            make_bam_time = timedelta(0), 
            metrics = [], ins_metrics = {}):
        self.name = name
        self.build_ind_time = build_ind_time
        self.aln_time = aln_time
        self.make_bam_time = make_bam_time
        self.summ_time = build_ind_time + aln_time + make_bam_time
        self.metrics = metrics # [{metric -> value}] (because of paired, unpaired, first and second of pair metrics) 
        self.ins_metrics = ins_metrics

    def __str__(self):
        return ("Aligner: name = " + self.name 
            + "\nbuilding index time = " + str(self.build_ind_time)
            + ", alignment time = " + str(self.aln_time)
            + ", making bam time = " + str(self.make_bam_time)
            + ", summary time = " + str(self.summ_time)
            + "\n summary metrics = " + str(self.metrics)
            + "\n insert size metrics = " + str(self.ins_metrics))

    def empty(self):
        return (self.build_ind_time == timedelta(0) 
            and self.aln_time == timedelta(0)
            and self.make_bam_time == timedelta(0)
            and self.metrics == []
            and self.ins_metrics == {})

    def parse_aln_summ_metrics(self, filename):
        f = open(filename)
        keys = []
        self.metrics = []
        for l in f:
            if (l[0] == '#' or l.strip() == ""):
                continue
            values = l.strip().split('\t')
            if keys != []:
                self.metrics.append(dict(zip(keys, l.strip().split('\t')))) 
            else:
                keys = values

    def parse_aln_ins_metrics(self, filename):
        f = open(filename)
        keys = []
        for l in f:
            if (l[0] == '#' or l.strip() == ""):
                continue
            values = l.strip().split('\t')
            if keys != []:
                self.ins_metrics = dict(zip(keys, l.strip().split('\t')))
                return 
            else:
                keys = values

    def parse_mapping_time_logs(self, filename):
        f = open(filename)
        prev = 0
        for l in f:
            tmp = l.strip().split('\t')
            curr = datetime.fromtimestamp(int(tmp[-1]))
            if (prev != 0):
                time = curr - prev
                if tmp[0] == '1':
                    self.build_ind_time = time
                elif tmp[0] == '2':
                    self.aln_time = time
                elif tmp[0] == '3':
                    self.make_bam_time = time
            prev = curr
        self.summ_time = self.build_ind_time + self.aln_time + self.make_bam_time

    def fill_metrics(self, curr_aln_dir):
        aln_files = os.listdir(curr_aln_dir)
        if ("metrics.alignment_summary_metrics" in aln_files):
            self.parse_aln_summ_metrics(curr_aln_dir + "/metrics.alignment_summary_metrics")
        if ("metrics.insert_size_metrics" in aln_files):
            self.parse_aln_ins_metrics(curr_aln_dir + "/metrics.insert_size_metrics")
        if ("mapping_time.log" in aln_files):
            self.parse_mapping_time_logs(curr_aln_dir + "/mapping_time.log")

class VCFStats:

    def __init__(self):
        self.ID = []
        self.SN = []
        self.SiS = []
        self.AF = []
        self.empty = True

    def __str__(self):
        return ("ID = " + str(self.ID) 
            + "\nSN = " + str(self.SN)
            + "\nSiS = " + str(self.SiS) 
            + "\nAF = " + str(self.AF))
    def clear(self):
        self.ID = []
        self.SN = []
        self.SiS = []
        self.AF = []
        self.empty = True

    def is_empty(self):
        return (self.ID == [] 
            and self.SN == []
            and self.SiS == []
            and self.AF == []
            and self.empty)

    def parse_vcf_stats(self, filename, vcf_name):
        f = open(filename)
        keys = []
        for l in f:
            if (l[0] == "#"):
                keys = l.strip()[2:].split('\t')[1:]
                for i in range(len(keys)):
                    keys[i] = keys[i][3:]
            elif (l.strip() != ""):
                values = l.strip().split('\t')
                if (values[0] == "ID"):
                    self.ID.append(values[2:])
                    self.SN.append({})
                    self.SiS.append({})
                    self.AF.append([])
                elif (values[0] == "SN"):
                    if (values[2] != "number of samples:"):
                        self.SN[int(values[1])][values[2][:-1]] = values[-1]
                elif (values[0] == "SiS"):
                    self.SiS[int(values[1])] = dict(zip(keys[1:], values[2:]))
                elif (values[0] == "AF"):
                    self.AF[int(values[1])].append(dict(zip(keys[1:], values[2:])))
        for i in range(len(self.ID)):
            if (len(self.ID[i]) == 1):
                if (os.path.basename(self.ID[i][0]) == os.path.basename(vcf_name + ".gz")):
                    self.ID[i] = "FP"
                else:
                    self.ID[i] = "FN"
            else:
                self.ID[i] = "TP"
        if (not self.is_empty() and self.SN[0] == {}):
            self.clear()

class VarCaller:

    def __init__(self, name):
        self.name = name
        self.time = timedelta(0)
        self.aligner = ""
        self.single_stats = VCFStats()
        self.shared_stats = VCFStats()

    def __str__(self):
        return ("Variant caller: name = " + self.name + ", aligner = " + self.aligner
            + "\ntime = " + str(self.time)
            + "\nsingle statistics: \n" + str(self.single_stats)
            + "\nshared statistics: \n" + str(self.shared_stats))

    def empty(self):
        return (self.single_stats.is_empty()
            and self.shared_stats.is_empty())

    def parse_var_call_time_logs(self, filename):
        f = open(filename)
        prev = 0
        for l in f:
            tmp = l.strip().split('\t')
            curr = datetime.fromtimestamp(int(tmp[-1]))
            if (prev != 0):
                self.time = curr - prev
            prev = curr

    def index_vcf(self, filename, err_log, dirname):
        if filename.endswith('.vcf'):
            vcf_bgzip = open(filename + ".gz", "w")
            subprocess.call("vcf-sort " + filename 
                + " | " + dirname + "/lib/convert.py"
                + " | bgzip -c", 
                stdout=vcf_bgzip, stderr=err_log, shell=True)
            vcf_bgzip.close()
            filename += '.gz'
        subprocess.call("tabix -p vcf " + filename, stderr=err_log, shell=True)
        return filename

    def process_vcf(self, aligner_name, var_caller_name, res_dir, filename, golden, dirname):
        print aligner_name + " " + var_caller_name
        res_base = aligner_name + "." + var_caller_name
        vcf_name = os.path.basename(filename)
        err_log = open(res_dir + "/err.log", "a")
        filename = self.index_vcf(filename, err_log, dirname)
        vcf_chk = open(res_dir + "/" + res_base + ".single.vchk", "w", 0)
        if (subprocess.call("vcf-validator " + filename, stderr=err_log, shell=True) == 0):
            subprocess.call("vcf check " + filename, stdout=vcf_chk, stderr=err_log, shell=True)
            vcf_chk.close()
            self.single_stats.parse_vcf_stats(res_dir + "/" + res_base + ".single.vchk", vcf_name)
            if (golden != None):
                vcf_chk = open(res_dir + "/" + res_base + ".shared.vchk", "w")
                subprocess.call("vcf check " + filename + " " + golden, stdout=vcf_chk, stderr=err_log, shell=True)
                self.shared_stats.parse_vcf_stats(res_dir + "/" + res_base + ".shared.vchk", vcf_name)

    def fill_stats(self, curr_var_dir, aln_name, var_caller_name, in_dir, name, vcf, dirname):
        if (os.path.exists(curr_var_dir + "/snp_time.log")):
            self.parse_var_call_time_logs(curr_var_dir + "/snp_time.log")

        if (os.path.exists(curr_var_dir + "/" + name + ".vcf")):
            self.process_vcf(aln_name, var_caller_name, in_dir + "/test_stats", curr_var_dir + "/" + name + ".vcf", vcf, dirname)

