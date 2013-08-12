import numpy as np
import matplotlib.pyplot as plt
width = 0.3

colors = ['blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black', '#99ff99', '#ccccff', '#ffcc00']

aln_metrics_defs = {'CATEGORY':'Category', 'PF_READS_ALIGNED':'Aligned reads', 'PF_HQ_ALIGNED_READS':'Reads aligned with high qual', 'READS_ALIGNED_IN_PAIRS':'Reads aligned in pairs', 'PF_MISMATCH_RATE':'Mismatch rate', 'PF_INDEL_RATE':'Indel rate'}
ins_size_metrics_defs = {'MEAN_INSERT_SIZE':'Mean inset size', 'MIN_INSERT_SIZE':'Min insert size', 'MAX_INSERT_SIZE':'Max insert size', 'STANDARD_DEVIATION':'Standard deviation'}

def remove_empty_metrics(aligners):
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

def begin_document(tex_f):
    tex_f.write(
    r'''\documentclass[a4paper]{report}
\usepackage [T2A]{fontenc}
\usepackage{graphicx}
\usepackage[utf8]{inputenc}
\usepackage[english]{babel}
\usepackage{multirow}
\usepackage{float}
\usepackage{pdflscape}
\usepackage{geometry}
\usepackage{longtable}
\geometry{verbose,a4paper,tmargin=2cm,bmargin=2cm,lmargin=2.5cm,rmargin=1.5cm}
\begin{document}
\hyphenpenalty=9999
''')

def end_document(tex_f):
    tex_f.write(r'''\end{document}''')

def new_section(tex_f, name):
    tex_f.write(r'''
\section{''' + name +'''}
''')

def add_picture(tex_f, path, scale):
    tex_f.write(r'''\begin{figure}[h]
    \begin{center}
    \includegraphics[scale=''' + str(scale) + "]{" + path + r'''}
    \end{center}
    \end{figure}
    ''')

def aligners_time_plot(aligners, in_dir):
    fig = plt.figure()
    plot = fig.add_subplot(111)

    rects = []
    i = 0 
    for aln in aligners:
        rects.append([plot.bar(i * width, aln.summ_time.seconds, width - 0.02, color=colors[i%len(colors)]), aln.name])
        i += 1

    plot.set_ylabel(u"Time, s")
    plot.set_title(u"Alignment time")
    plot.set_xticklabels([])
    box = plot.get_position()
    plot.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    plot.legend(tuple([tmp[0][0] for tmp in rects]), tuple([tmp[1] for tmp in rects]), loc='center left', bbox_to_anchor=(1, 0.5))
    fig.savefig(in_dir + "/test_stats/plots/aln_summ_time.png")

def output_aligners_time(aligners, tex_f, in_dir):
    tex_f.write(r'''
    \begin{table}[h]
    \begin{center}
    \begin{tabular}{|p{2cm}|p{3cm}p{3cm}p{3cm}p{3cm}|}
    \hline\hline
    ''')
    tex_f.write(r'''Name & Building index time & Alignment time & Bam making time & Summary time\\ [0.5ex]
    \hline\hline
    ''')

    for aln in aligners:
        tex_f.write(aln.name + " & " + str(aln.build_ind_time) + " & " + str(aln.aln_time) + " & " + str(aln.make_bam_time) + " & " + str(aln.summ_time) + "\\\\" + "\n" + "\\hline\n")

    tex_f.write(r'''\end{tabular}
    \end{center}
    \end{table}
    ''')

    aligners_time_plot(aligners, in_dir)

def aligners_stats_plots(aligners, in_dir):
    fig = plt.figure()
    plot = fig.add_subplot(111)

    i = 0
    rects = []
    tmp_data = {'aln': [], 'mismatch': [], 'indel': []}
    ind = np.arange(len(aligners))

    for aln in aligners:
        if (aln.metrics != []):
            for metric in aln.metrics:
                if metric == aln.metrics[-1]:
                    rects.append([
                        plot.bar(i * width, 
                        int(metric['PF_READS_ALIGNED']), width-0.02, color="#E0E0E0"),
                        plot.bar(i * width, int(metric['PF_HQ_ALIGNED_READS']), 
                        width - 0.02, color=colors[i%len(colors)]), 
                        aln.name])
                    tmp_data['aln'].append(aln.name)
                    tmp_data['mismatch'].append(float(metric['PF_MISMATCH_RATE']))
                    tmp_data['indel'].append(float(metric['PF_INDEL_RATE']))
                    i += 1
        else:
            tmp_data['aln'].append(aln.name)
            tmp_data['mismatch'].append(0)
            tmp_data['indel'].append(0)
            i += 1

    plot.set_ylabel(u"Number of aligned reads")
    plot.set_title(u"Aligned reads")
    plot.set_xticklabels([])
    box = plot.get_position()
    plot.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    plot.legend(tuple([tmp[1][0] for tmp in rects] + [rects[0][0][0]]), tuple([tmp[2] for tmp in rects] + [u'LQ']), loc='center left', bbox_to_anchor=(1, 0.5))
    fig.savefig(in_dir + "/test_stats/plots/aln_aligned_reads.png")

    fig.clear()
    plot = fig.add_subplot(111)

    rects = []
    print "height is " + str(tmp_data['mismatch'])
    rects.append(plot.bar(ind, tmp_data['mismatch'], width - 0.02, color=colors[0]))
    rects.append(plot.bar(ind + width, tmp_data['indel'], width - 0.02, color=colors[1]))

    plot.set_ylabel("Rates")
    plot.set_title("Mismatch and indel rates")
    plot.set_xticks(ind + width * 2)
    plot.set_xticklabels(tmp_data['aln'])
    box = plot.get_position()
    plot.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    lgd = plot.legend((rects[0][0], rects[1][0]), ('mismatch rate', 'indel rate'), loc='center left', bbox_to_anchor=(1, 0.5))
    fig.autofmt_xdate()
    fig.savefig(in_dir + "/test_stats/plots/aln_snp_indel_rate.png", bbox_extra_artists=(lgd,), bbox_inches='tight')

def output_aligners_stats(aligners, tex_f, in_dir):
    tex_f.write(r'''
    \begin{landscape}
    \begin{table}[h]
    \begin{center}
    \begin{tabular}{|p{1.5cm}|''' 
        + (("p{" + str(15.0 / (len(aln_metrics_defs) + len(ins_size_metrics_defs))) 
        + "cm}") * (len(aln_metrics_defs) + len(ins_size_metrics_defs)))
        + r'''|}
    \hline\hline
    Name & ''')

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
                        if metric_name == 'CATEGORY':
                            tex_f.write(" & \\verb|" + metric[metric_name][:6] + "|")
                        else:
                            tex_f.write(" & \\verb|" + str(round(float(metric[metric_name]), 4)) + "|")
                    else:
                        tex_f.write(" & 0")
                if (metric['CATEGORY'] == "PAIR"):
                    for metric_name in ins_size_metrics_defs.keys():
                        tex_f.write(" & " + str(round(float(aln.ins_metrics[metric_name]), 4)))
                else:
                    tex_f.write(" & " * len(ins_size_metrics_defs.keys()))
                tex_f.write(" \\\\\n")
            tex_f.write("\\hline\n")

    tex_f.write("\\hline\n")
    tex_f.write(r'''\end{tabular}
    \end{center}
    \end{table}
    \end{landscape}
    ''')

    aligners_stats_plots(aligners, in_dir)
    add_picture(tex_f, "plots/aln_aligned_reads.png", 0.5)
    add_picture(tex_f, "plots/aln_snp_indel_rate.png", 0.5)

def vc_single_plot(fig, plot, ylabel, title, var_caller_names, ind, aln_names, rects, i, saveto):
    plot.set_ylabel(ylabel)
    plot.set_title(title)
    plot.set_xticks(ind + width * len(var_caller_names) / 2)
    plot.set_xticklabels(aln_names)
    box = plot.get_position()
    plot.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    lgd = plot.legend(tuple([tmp[0] for tmp in rects[i::3]]), var_caller_names, loc='center left', bbox_to_anchor=(1, 0.5))
    fig.autofmt_xdate()
    fig.savefig(saveto, bbox_extra_artists=(lgd,), bbox_inches='tight')
 
def vc_single_plots(var_callers, in_dir, var_caller_names):
    snp_fig = plt.figure()
    snp_plot = snp_fig.add_subplot(111)
    indel_fig = plt.figure()
    indel_plot = indel_fig.add_subplot(111)
    time_fig = plt.figure()
    time_plot = time_fig.add_subplot(111)

    rects = []
    tmp_data = {'aln':[], 'snp':{}, 'indel':{}, 'time':{}}
    for vc in var_caller_names:
        tmp_data['snp'][vc] = [0] * len(var_callers.keys())
        tmp_data['indel'][vc] = [0] * len(var_callers.keys())
        tmp_data['time'][vc] = [0] * len(var_callers.keys())
    i = 0
    ind = np.arange(0, len(var_callers.keys()) * width * len(var_caller_names) - 0.0001, len(var_caller_names) * width)
    print len(var_caller_names)
    print len(var_callers.keys())
    print ind

    for aligner in sorted(var_callers.keys()):
        if (not len(var_callers[aligner])):
            continue
        tmp_data['aln'].append(aligner)
        for var_caller in sorted(var_callers[aligner], key = lambda vc: vc.name):
            if (not var_caller.single_stats.is_empty()):
                tmp_data['snp'][var_caller.name][i] = int(var_caller.single_stats.SN[0]["number of SNPs"])
                tmp_data['indel'][var_caller.name][i] = int(var_caller.single_stats.SN[0]["number of indels"])
                tmp_data['time'][var_caller.name][i] = var_caller.time.seconds
        i += 1

    i = 0
    for var_caller in var_caller_names:
        print var_caller
        print tmp_data['snp'][var_caller]
        rects.append(snp_plot.bar(ind + i * width, tmp_data['snp'][var_caller], width - 0.02, color = colors[i]))
        rects.append(indel_plot.bar(ind + i * width, tmp_data['indel'][var_caller], width - 0.02, color = colors[i]))
        rects.append(time_plot.bar(ind + i * width, tmp_data['time'][var_caller], width - 0.02, color = colors[i]))
        i += 1

    vc_single_plot(snp_fig, snp_plot, 'Number of SNP', 'Number of SNP', var_caller_names, ind, tmp_data['aln'], rects, 0, in_dir + '/test_stats/plots/vc_aln_snp.png')
    vc_single_plot(indel_fig, indel_plot, 'Number of indels', 'Number of indels', var_caller_names, ind, tmp_data['aln'], rects, 1, in_dir + '/test_stats/plots/vc_aln_indel.png')
    vc_single_plot(time_fig, time_plot, 'Time, s', 'Variant calling time', var_caller_names, ind, tmp_data['aln'], rects, 2, in_dir + '/test_stats/plots/vc_aln_time.png')

def vc_single_stats(var_callers, tex_f, in_dir, var_caller_names):

    print var_callers.keys()
    for aligner in sorted(var_callers.keys()):
        if (not len(var_callers[aligner])):
            continue
        tex_f.write("\\subsection{" + aligner + r'''}
    \begin{table}[H]
    \begin{center}
    \begin{tabular}{|p{1.5cm}|''' 
            + (("p{1.5cm}") * 6) + "|p{1.5cm}|}")

        tex_f.write(r'''
    \hline\hline
    Name & SNP & MNP & Indel & Other & Multiallelic sites & Ts/tv & Time \\ [0.5ex]
    \hline\hline
    ''')

        for var_caller in sorted(var_callers[aligner], key = lambda vc: vc.name):
            if (not var_caller.single_stats.is_empty()):
                tex_f.write(var_caller.name + " & " 
                    + ((var_caller.single_stats.SN[0]["number of SNPs"] + " & "
                    + var_caller.single_stats.SN[0]["number of MNPs"] + " & "
                    + var_caller.single_stats.SN[0]["number of indels"] + " & "
                    + var_caller.single_stats.SN[0]["number of others"] + " & "
                    + var_caller.single_stats.SN[0]["number of multiallelic sites"] + " & " 
                    + var_caller.single_stats.SN[0]["ts/tv"]) 
                    if (len(var_caller.single_stats.SN) != 0) 
                    else "0 & 0 & 0 & 0 & 0 & 0")
                    + " & "
                    + str(var_caller.time)  + "\\\\\n"
                    + "\\hline\n")

        tex_f.write("\\hline\n")
        tex_f.write(r'''
    \end{tabular}
    \end{center}
    \end{table}
    ''')
    vc_single_plots(var_callers, in_dir, var_caller_names)
    add_picture(tex_f, "plots/vc_aln_snp.png", 0.5)
    add_picture(tex_f, "plots/vc_aln_indel.png", 0.5)
    add_picture(tex_f, "plots/vc_aln_time.png", 0.5)

def vc_shared_plot(fig, plot, ylabel, title, var_caller_names, ind, aln_names, rects, i, saveto):
    plot.set_ylabel(ylabel)
    plot.set_title(title)
    plot.set_xticks(ind + (width*len(var_caller_names) / 2))
    plot.set_xticklabels(aln_names)
    box = plot.get_position()
    plot.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    lgd = plot.legend(tuple([tmp[0] for tmp in rects[i + 1::4]] + [rects[i][0]]), var_caller_names + ["FP"], loc='center left', bbox_to_anchor=(1, 0.5))

    fig.autofmt_xdate()
    fig.savefig(saveto, bbox_extra_artists=(lgd,), bbox_inches='tight')

def vc_shared_plots(var_callers, in_dir, var_caller_names):
    snp_fig = plt.figure()
    snp_plot = snp_fig.add_subplot(111)
    indel_fig = plt.figure()
    indel_plot = indel_fig.add_subplot(111)
    rects = []
    tmp_data = {'aln':[], 'snp':{}, 'indel':{}}

    for vc in var_caller_names:
        tmp_data['snp'][vc] = [[0] * len(var_callers.keys()), [0] * len(var_callers.keys())]
        tmp_data['indel'][vc] = [[0] * len(var_callers.keys()), [0] * len(var_callers.keys())]
    i = 0
    ind = np.arange(0, len(var_callers.keys()) * width * len(var_caller_names) - 0.0001, len(var_caller_names) * width)
    print ind
    print len(var_callers.keys()) * width * len(var_caller_names)
    print len(var_caller_names)

    for aligner in sorted(var_callers.keys()):
        if (not len(var_callers[aligner])):
            continue
        tmp_data['aln'].append(aligner)

        for var_caller in sorted(var_callers[aligner], key = lambda vc: vc.name):
            if (var_caller.shared_stats.is_empty() 
                    or len(var_caller.shared_stats.SN) < 3):
                continue

            for j in [2, 0, 1]:
                if (var_caller.shared_stats.ID[j] == "TP"):
                    tmp_data['snp'][var_caller.name][0][i] = int(var_caller.shared_stats.SN[j]["number of SNPs"])
                    tmp_data['indel'][var_caller.name][0][i] = int(var_caller.shared_stats.SN[j]["number of indels"])
                elif (var_caller.shared_stats.ID[j] == "FP"):
                    tmp_data['snp'][var_caller.name][1][i] = int(var_caller.shared_stats.SN[j]["number of SNPs"])
                    tmp_data['indel'][var_caller.name][1][i] = int(var_caller.shared_stats.SN[j]["number of indels"])
        i += 1

    i = 0
    for var_caller in var_caller_names:
        rects.append(snp_plot.bar(ind + (i * width), 
            [tmp_data['snp'][var_caller][0][j] + tmp_data['snp'][var_caller][1][j] for j in xrange(len(tmp_data['aln']))], width - 0.02, color = "#E0E0E0"))
        rects.append(snp_plot.bar(ind + (i * width), tmp_data['snp'][var_caller][0], width - 0.02, color = colors[i]))
        rects.append(indel_plot.bar(ind + (i * width), 
            [tmp_data['indel'][var_caller][0][j] + tmp_data['indel'][var_caller][1][j] for j in xrange(len(tmp_data['aln']))], width - 0.02, color = "#E0E0E0"))
        rects.append(indel_plot.bar(ind + (i * width), tmp_data['indel'][var_caller][0], width - 0.02, color = colors[i]))
        i += 1

    vc_shared_plot(snp_fig, snp_plot, 'Number of SNP', 'Number of SNP', var_caller_names, ind, tmp_data['aln'], rects, 0, in_dir + "/test_stats/plots/vc_snp.png")
    vc_shared_plot(indel_fig, indel_plot, 'Number of indels', 'Number of indels', var_caller_names, ind, tmp_data['aln'], rects, 2, in_dir + "/test_stats/plots/vc_indel.png")

def vc_shared_stats(var_callers, tex_f, in_dir, var_caller_names, golden_vcf):
    for aligner in sorted(var_callers.keys()):
        if (not len(var_callers[aligner])):
            continue
        tex_f.write("\\subsection{" + aligner + r'''}
\begin{table}[H]
\begin{center}
\begin{tabular}{|p{1.2cm}|p{1.2cm}|''' 
            + (("p{2.2cm}") * 5) + "p{0.8cm}" + "|}")
        tex_f.write(r'''
\hline\hline
Name & Category & SNP & MNP & Indel & Other& Multiallelic sites & Ts/tv \\ [0.5ex]
\hline\hline
''')

        for var_caller in sorted(var_callers[aligner], key = lambda vc: vc.name):
            if (var_caller.shared_stats.is_empty() 
                    or len(var_caller.shared_stats.SN) < 3):
                print "Shared stats is empty"
                continue

            tex_f.write("\\multirow{" + str(len(var_caller.shared_stats.SN)) + "}{*}{" + var_caller.name + "}")
            for j in [2, 0, 1]:
                tex_f.write(" & " + var_caller.shared_stats.ID[j] + " & "
                    + var_caller.shared_stats.SN[j]["number of SNPs"] 
                    + ("" if golden_vcf.single_stats.SN[0]["number of SNPs"] == "0" 
                        else ("(" + str(round(float(var_caller.shared_stats.SN[j]["number of SNPs"])
                        / float(golden_vcf.single_stats.SN[0]["number of SNPs"]) * 100, 2)) 
                    + "\%)")) + " & "
                    + var_caller.shared_stats.SN[j]["number of MNPs"] 
                    + ("" if golden_vcf.single_stats.SN[0]["number of MNPs"] == "0" 
                        else ("(" + str(round(float(var_caller.shared_stats.SN[j]["number of MNPs"])
                        / float(golden_vcf.single_stats.SN[0]["number of MNPs"]) * 100, 2)) 
                    + "\%)")) + " & "
                    + var_caller.shared_stats.SN[j]["number of indels"] 
                    + ("" if golden_vcf.single_stats.SN[0]["number of indels"] == "0"
                        else ("(" + str(round(float(var_caller.shared_stats.SN[j]["number of indels"])
                        / float(golden_vcf.single_stats.SN[0]["number of indels"]) * 100, 2))
                    + "\%)")) + " & "
                    + var_caller.shared_stats.SN[j]["number of others"] 
                    + ("" if golden_vcf.single_stats.SN[0]["number of others"] == "0"
                        else ("(" + str(round(float(var_caller.shared_stats.SN[j]["number of others"])
                        / float(golden_vcf.single_stats.SN[0]["number of others"]) * 100, 2))
                    + "\%)")) + " & "
                    + var_caller.shared_stats.SN[j]["number of multiallelic sites"] 
                    + ("" if golden_vcf.single_stats.SN[0]["number of multiallelic sites"] == "0"
                        else ("(" + str(round(float(var_caller.shared_stats.SN[j]["number of multiallelic sites"])
                        / float(golden_vcf.single_stats.SN[0]["number of multiallelic sites"]) * 100, 2))
                    + "\%)")) + " & "
                    + var_caller.shared_stats.SN[j]["ts/tv"] + "\\\\\n")

            tex_f.write("\\hline\n")

        tex_f.write("\\hline\n")
        tex_f.write(r'''
\end{tabular}
\end{center}
\end{table}
''')
    vc_shared_plots(var_callers, in_dir, var_caller_names)
    add_picture(tex_f, "plots/vc_snp.png", 0.7)
    add_picture(tex_f, "plots/vc_indel.png", 0.7)
