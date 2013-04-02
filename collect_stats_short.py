#!/usr/bin/python
# -*- coding: utf-8 -*-

from argparse import ArgumentParser
from datetime import datetime
from datetime import timedelta
import os.path
import os
import subprocess
import shutil
import sys

aln_metrics_defs = {'CATEGORY':'Категория', 'PF_READS_ALIGNED':'Число выровненных ридов', 'PF_HQ_ALIGNED_READS':'Число ридов, выровненных с высоким качеством', 'READS_ALIGNED_IN_PAIRS':'Число ридов, выровненных в парах', 'PF_MISMATCH_RATE':'Число ошибок заменой', 'PF_INDEL_RATE':'Число ошибок вставкой'}
ins_size_metrics_defs = {'MEAN_INSERT_SIZE':'Среднее расстояние вставки', 'MIN_INSERT_SIZE':'Минимальное расстояние вставки', 'MAX_INSERT_SIZE':'Максимальное расстояние вставки', 'STANDARD_DEVIATION':'Среднеквадратичное отклонение расстояния вставки'}

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


class VCFStats:

    def __init__(self):
        self.ID = []
        self.SN = []
        self.SiS = []
        self.AF = []

    def __str__(self):
        return ("ID = " + str(self.ID) 
            + "\nSN = " + str(self.SN)
            + "\nSiS = " + str(self.SiS) 
            + "\nAF = " + str(self.AF))

    def empty(self):
        return (self.ID == [] 
            and self.SN == []
            and self.SiS == []
            and self.AF == [])

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
        return (self.single_stats.empty()
            and self.shared_stats.empty())

def parse_aln_summ_metrics(filename):
    f = open(filename)
    metrics = []
    keys = []
    for l in f:
        if (l[0] == '#' or l.strip() == ""):
            continue
        values = l.strip().split('\t')
        if keys != []:
            metrics.append(dict(zip(keys, l.strip().split('\t')))) 
        else:
            keys = values
    return metrics

def parse_aln_ins_metrics(filename):
    f = open(filename)
    metrics = {}
    keys = []
    for l in f:
        if (l[0] == '#' or l.strip() == ""):
            continue
        values = l.strip().split('\t')
        if keys != []:
            return dict(zip(keys, l.strip().split('\t')))
        else:
            keys = values

def parse_mapping_time_logs(filename, aligner):
    f = open(filename)
    prev = 0
    for l in f:
        tmp = l.strip().split('\t')
        curr = datetime.fromtimestamp(int(tmp[-1]))
        if (prev != 0):
            time = curr - prev
            if tmp[0] == '1':
                aligner.build_ind_time = time
            elif tmp[0] == '2':
                aligner.aln_time = time
            elif tmp[0] == '3':
                aligner.make_bam_time = time
        prev = curr
    aligner.summ_time = aligner.build_ind_time + aligner.aln_time + aligner.make_bam_time
    return aligner

def parse_var_call_time_logs(filename, var_caller):
    f = open(filename)
    prev = 0
    for l in f:
        tmp = l.strip().split('\t')
        curr = datetime.fromtimestamp(int(tmp[-1]))
        if (prev != 0):
            var_caller.time = curr - prev
        prev = curr
    return var_caller

def parse_vcf_stats(filename, vcf_name):
    f = open(filename)
    vcf_stats = VCFStats()
    keys = []
    for l in f:
        if (l[0] == "#"):
            keys = l.strip()[2:].split('\t')[1:]
            for i in range(len(keys)):
                keys[i] = keys[i][3:]
        elif (l.strip() != ""):
            values = l.strip().split('\t')
            if (values[0] == "ID"):
                vcf_stats.ID.append(values[2:])
                vcf_stats.SN.append({})
                vcf_stats.SiS.append({})
                vcf_stats.AF.append([])
            elif (values[0] == "SN"):
                if (values[2] != "number of samples:"):
                    vcf_stats.SN[int(values[1])][values[2][:-1]] = values[-1]
            elif (values[0] == "SiS"):
                vcf_stats.SiS[int(values[1])] = dict(zip(keys[1:], values[2:]))
            elif (values[0] == "AF"):
                vcf_stats.AF[int(values[1])].append(dict(zip(keys[1:], values[2:])))
    for i in range(len(vcf_stats.ID)):
        if (len(vcf_stats.ID[i]) == 1):
#            print vcf_stats.ID[i][0] + " " + vcf_name
            if (os.path.basename(vcf_stats.ID[i][0]) == os.path.basename(vcf_name + ".gz")):
                vcf_stats.ID[i] = "результат"
            else:
                vcf_stats.ID[i] = "эталон"
        else:
            vcf_stats.ID[i] = "общий"
    if (not vcf_stats.empty() and vcf_stats.SN[0] == {}):
        return VCFStats()
    return vcf_stats

def process_vcf(aligner_name, var_caller_name, res_dir, filename, golden, var_caller):
    print aligner_name + " " + var_caller_name
    res_base = aligner_name + "." + var_caller_name
    vcf_name = os.path.basename(filename)
    err_log = open(res_dir + "/err.log", "a")
    vcf_bgzip = open(filename + ".gz", "w")
    subprocess.call("vcf-sort " + filename 
        + " | " + os.path.dirname(sys.argv[0]) + "/convert.py"
        + " | bgzip -c", 
        stdout=vcf_bgzip, stderr=err_log, shell=True)
    vcf_bgzip.close()
    subprocess.call("tabix -p vcf " + filename + ".gz", stderr=err_log, shell=True)
    vcf_chk = open(res_dir + "/" + res_base + ".single.vchk", "w", 0)
    if (subprocess.call("vcf check " + filename + ".gz | cat ", stdout=vcf_chk, stderr=err_log, shell=True) == 0):
        vcf_chk.close()
        subprocess.call("plot-vcfcheck -p " + res_dir + "/plots/" 
            + res_base + ".single" + " " 
            + res_dir + "/" + res_base + ".single.vchk", stderr=err_log,
            shell=True)
        var_caller.single_stats = parse_vcf_stats(res_dir + "/" + res_base + ".single.vchk", os.path.basename(filename))
        if (golden != None):
            vcf_chk = open(res_dir + "/" + res_base + ".shared.vchk", "w")
            subprocess.call("vcf check " + filename + ".gz " + golden, stdout=vcf_chk, stderr=err_log, shell=True)
            subprocess.call("plot-vcfcheck -p " + res_dir + "/plots/" 
                + res_base + ".shared" + " " 
                + res_dir + "/" + res_base + ".shared.vchk", stderr=err_log,
                shell=True)
            var_caller.shared_stats = parse_vcf_stats(res_dir + "/" + res_base + ".shared.vchk", os.path.basename(filename))
    return var_caller

parser = ArgumentParser(description="Collect statistics from test.py")
parser.add_argument("-i", action="store", dest="in_dir", help="directory with test.py script results", required=True)
parser.add_argument("--vcf", action="store", dest="vcf", help="golden file with variants")
parser.add_argument("-n", "--name", action="store", dest="name", help="base name for vcf files", required=True)

args = parser.parse_args()

print args

if (args.vcf != None and not os.path.exists(args.vcf)):
    print "Golden file not exists"
    sys.exit(-1)

aligners = [] # [Aligner1 Aligner2 ...]
var_callers = {} # {aligner_name->[var_caller1 ...] ... }
aln_dirs = os.listdir(args.in_dir)

if (os.path.exists(args.in_dir + "/test_stats")):
    shutil.rmtree(args.in_dir + "/test_stats")
os.makedirs(args.in_dir + "/test_stats")
os.makedirs(args.in_dir + "/test_stats/plots")

for aln_name in aln_dirs:
    curr_aln_dir = args.in_dir + "/" + aln_name
    if (not os.path.isdir(curr_aln_dir)):
        continue

    aln_files = os.listdir(curr_aln_dir)

    curr_aligner = Aligner(aln_name)

    if ("metrics.alignment_summary_metrics" in aln_files):
        curr_aligner.metrics = parse_aln_summ_metrics(curr_aln_dir + "/metrics.alignment_summary_metrics")
    if ("metrics.insert_size_metrics" in aln_files):
        curr_aligner.ins_metrics = parse_aln_ins_metrics(curr_aln_dir + "/metrics.insert_size_metrics")
    if ("mapping_time.log" in aln_files):
        curr_aligner = parse_mapping_time_logs(curr_aln_dir + "/mapping_time.log", curr_aligner)

    if (not curr_aligner.empty()):
        aligners.append(curr_aligner)

    for var_caller_name in aln_files:
        curr_var_dir = curr_aln_dir + "/" + var_caller_name
        if (not os.path.isdir(curr_var_dir)):
            continue

        curr_var_caller = VarCaller(var_caller_name)
        
        if (os.path.exists(curr_var_dir + "/snp_time.log")):
            curr_var_caller = parse_var_call_time_logs(curr_var_dir + "/snp_time.log", curr_var_caller)

        if (os.path.exists(curr_var_dir + "/" + args.name + ".vcf")):
            curr_var_caller = process_vcf(aln_name, var_caller_name, args.in_dir + "/test_stats", curr_var_dir + "/" + args.name + ".vcf", args.vcf, curr_var_caller)

        if (not curr_var_caller.empty()):
            if (aln_name not in var_callers):
                var_callers[aln_name] = [curr_var_caller]
            else:
                var_callers[aln_name].append(curr_var_caller)

def_metrics = []
for aln in aligners:
    for m in aln.metrics:
        for key in m.keys():
            if (m[key] != '0' and key not in def_metrics):
                def_metrics.append(key)
    for key in aln.ins_metrics.keys():
        if (aln.ins_metrics[key] != '0' and key not in def_metrics):
            def_metrics.append(key)

for key in aln_metrics_defs.keys():
    if key not in def_metrics:
        del aln_metrics_defs[key]
for key in ins_size_metrics_defs.keys():
    if key not in def_metrics:
        del ins_size_metrics_defs[key]

#for aln in aligners:
#    print str(aln)
#
#for aln in var_callers.keys():
#    print aln + ":"
#    for var_c in var_callers[aln]:
#        print "\t" + str(var_c)

# Начало вывода статистики

tex_f = open(args.in_dir + "/test_stats/stats.tex", "w")
tex_f.write(
r'''\documentclass[a4paper]{article}
\usepackage [T2A]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage[russian]{babel}
\usepackage{multirow}
\usepackage{float}
\usepackage{geometry}
\geometry{verbose,a4paper,tmargin=2cm,bmargin=2cm,lmargin=2.5cm,rmargin=1.5cm}
\begin{document}
\hyphenpenalty=9999
''')

# Статистика по aligners

if (len(aligners)):
    tex_f.write(r'''
\section{Результаты работы инструментов для выравнивания}
''')

    tex_f.write(r'''
\begin{table}[H]
\caption{Сравнительная таблица продолжительности работы инструментов для выравнивания}
\begin{center}
\begin{tabular}{|p{2cm}|p{3cm}p{3cm}p{3cm}p{3cm}|}
\hline\hline
''')
    tex_f.write(r'''Название & Время построения индекса & Время выравнивания & Время построения сортированного bam файла & Общее время \\ [0.5ex]
\hline\hline
''')

    for aln in aligners:
        tex_f.write(aln.name + " & " + str(aln.build_ind_time) + " & " + str(aln.aln_time) + " & " + str(aln.make_bam_time) + " & " + str(aln.summ_time) + "\\\\" + "\n" + "\\hline\n")
    tex_f.write(r'''\end{tabular}
\end{center}
\end{table}
''')

    tex_f.write(r'''
\begin{table}[H]
\caption{Сравнительная таблица результатов работы инструментов для выравнивания}
\begin{center}
\begin{tabular}{|p{2cm}|''' 
        + (("p{" + str(12.0 / (len(aln_metrics_defs) + len(ins_size_metrics_defs))) 
        + "cm}") * (len(aln_metrics_defs) + len(ins_size_metrics_defs)))
        + r'''|}
\hline\hline
Название & ''')

    for metric_name in aln_metrics_defs.values():
        tex_f.write(metric_name)
        if (metric_name != aln_metrics_defs.values()[-1] or ins_size_metrics_defs != {}):
            tex_f.write(" & ")
        else:
            tex_f.write(" \\\\ [0.5ex]\n")
    for metric in ins_size_metrics_defs.values():
        tex_f.write(metric)
        if (metric != ins_size_metrics_defs.values()[-1]):
            tex_f.write(" & ")
        else:
            tex_f.write(" \\\\ [0.5ex]\n")
    tex_f.write("\\hline\\hline\n")

    for aln in aligners:
        if (aln.metrics != []):
            tex_f.write("\\multirow{" + str(len(aln.metrics)) + "}{*}{" + aln.name + "}")
            for metric in aln.metrics:
                for metric_name in aln_metrics_defs.keys():
                    if (metric_name in metric.keys()):
                        tex_f.write(" & \\verb|" + metric[metric_name] + "|")
                    else:
                        tex_f.write(" & 0")
                if (metric['CATEGORY'] == "PAIR"):
                    for metric_name in ins_size_metrics_defs.keys():
                        tex_f.write(" & " + aln.ins_metrics[metric_name])
                else:
                    tex_f.write(" & " * len(ins_size_metrics_defs.keys()))
                tex_f.write(" \\\\\n\\hline\n")

    tex_f.write("\\hline\n")
    tex_f.write(r'''\end{tabular}
\end{center}
\end{table}
''')

# Статистика variant callers

tex_f.write(r'''
\section{Результаты работы инструментов для поиска полиморфизмов}
''')

print var_callers.keys()
for aligner in var_callers.keys():
    if (not len(var_callers[aligner])):
        continue
    tex_f.write("\\subsection{" + aligner + r'''}
\begin{table}[H]
\caption{Общая таблица результатов работы инструментов для поиска полиморфизмов}
\begin{center}
\begin{tabular}{|p{1.5cm}|''' 
        + (("p{1.5cm}") * 6) + "|p{1.5cm}|}")

    tex_f.write(r'''
\hline\hline
Название & Число SNP & Число MNP & Число вставок & Число других полиморфизмов & Число multiallelic sites & ts/tv & Время работы \\ [0.5ex]
\hline\hline
''')

    for var_caller in var_callers[aligner]:
        if (var_caller.single_stats.SN != [] and var_caller.single_stats.SN[0] != {}):
            tex_f.write(var_caller.name + " & " 
                + var_caller.single_stats.SN[0]["number of SNPs"] + " & "
                + var_caller.single_stats.SN[0]["number of MNPs"] + " & "
                + var_caller.single_stats.SN[0]["number of indels"] + " & "
                + var_caller.single_stats.SN[0]["number of others"] + " & "
                + var_caller.single_stats.SN[0]["number of multiallelic sites"] + " & " 
                + var_caller.single_stats.SN[0]["ts/tv"] + " & "
                + str(var_caller.time)  + "\\\\\n"
                + "\\hline\n")
    tex_f.write("\\hline\n")
    tex_f.write(r'''
\end{tabular}
\end{center}
\end{table}
''')

if args.vcf != None:

    tex_f.write(r'''
\section{Сравнение полученных результатов с эталонным}
''')

    for aligner in var_callers.keys():
        if (not len(var_callers[aligner])):
            continue
        tex_f.write("\\subsection{" + aligner + r'''}
\begin{table}[H]
\caption{Сравнение общих результатов работы инструментов для поиска полиморфизмов с эталонным результатом}
\begin{center}
\begin{tabular}{|p{1.5cm}|p{1.5cm}|''' 
            + (("p{1.5cm}") * 6) + "|}")

        tex_f.write(r'''
\hline\hline
Название & Категория & Число SNP & Число MNP & Число вставок & Число других полиморфизмов & Число multiallelic sites & ts/tv \\ [0.5ex]
\hline\hline
''')

        for var_caller in var_callers[aligner]:
            if var_caller.shared_stats.empty():
                continue
            tex_f.write("\\multirow{" + str(len(var_caller.shared_stats.SN)) + "}{*}{" + var_caller.name + "}")
            for i in range(len(var_caller.shared_stats.SN)):
                tex_f.write(" & " + var_caller.shared_stats.ID[i] + " & "
                    + var_caller.shared_stats.SN[i]["number of SNPs"] + " & "
                    + var_caller.shared_stats.SN[i]["number of MNPs"] + " & "
                    + var_caller.shared_stats.SN[i]["number of indels"] + " & "
                    + var_caller.shared_stats.SN[i]["number of others"] + " & "
                    + var_caller.shared_stats.SN[i]["number of multiallelic sites"] + " & " 
                    + var_caller.shared_stats.SN[i]["ts/tv"] + "\\\\\n")
            tex_f.write("\\hline\n")
        tex_f.write("\\hline\n")
        tex_f.write(r'''
\end{tabular}
\end{center}
\end{table}
''')

tex_f.write(r'''\end{document}''')
