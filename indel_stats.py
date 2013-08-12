#!/usr/bin/python

from argparse import ArgumentParser
from datetime import datetime
import os
import re
import sys
import shutil
from lib import stats_writer
from operator import itemgetter, attrgetter
from Bio import SeqIO

class Position:
    
    def __init__(self, pos, chr, unified_ref, unified_alt, diff):
        self.pos = pos
        self.chr = chr
        self.variants = []
        self.unified_ref = unified_ref
        self.unified_alt = unified_alt
        self.diff = diff

    def add_variant(self, variant):
        self.variants.append(variant)

def is_indel(fields):
    if len(fields[3]) != 1 or len(fields[4]) != 1 or fields[3] == '.' or fields[4] == '.':
        return True
    else:
        return False

def find_borders(seq, rep, begin):
    rep_len = len(rep)
    r = begin - 1
    l = begin
    for i in xrange(begin - rep_len, 0, -rep_len):
        if seq[i:i + rep_len] != rep:
            break
        r = i - 1
    for i in xrange(begin, len(seq), rep_len):
        if seq[i:i + rep_len] != rep:
            break
        l = i + rep_len
#    print "Full variant for " + rep + " is " + str(seq[r:l])
    return (r, l)

def find_repeat(rep):
    for i in xrange(1, len(rep)):
        is_rep = True
        for j in xrange(i, len(rep), i):
            if rep[j:j + i] != rep[0:i]:
                is_rep = False
                break
        if is_rep:
#            print "Old diff is " + rep + ", new - " + rep[0:i]
            return rep[0:i]
    return rep

def unify(pos, ref, alt, fn, ref_seq):
    ref_len = len(ref)
    alt_len = len(alt)
    if ref[0] != alt[0]:
        print "Malformed variant in " + fn + " - " + str(pos) + ": " + ref + " " + alt
        return (None, None, None, None)
    diff = ""
    diff_begin = min([ref_len, alt_len])
    if ref.startswith(alt):
        diff = ref[alt_len:]
    elif alt.startswith(ref):
        diff = alt[ref_len:]
    else:
        diff_begin = next(i for i in xrange(ref_len) if ref[i] != alt[i])
        if ref_len > alt_len:
            diff = ref[diff_begin:diff_begin - alt_len] 
            if ref[diff_begin - alt_len:] != alt[diff_begin:]:
                print "Malformed variant in " + fn + " - " + str(pos) + ": " + ref + " " + alt + ", diff - " + diff
                return (None, None, None, None)
        else:
            diff = alt[diff_begin:diff_begin - ref_len] 
            if ref[diff_begin:] != alt[diff_begin - ref_len:]:
                print "Malformed variant in " + fn + " - " + str(pos) + ": " + ref + " " + alt + ", diff - " + diff
                return (None, None, None, None)
#    print "Difference for " + fn + " - " + str(pos) + ": " + ref + " " + alt + " - " + diff + ", begin " + str(diff_begin)
    short_diff = find_repeat(diff)
    (r, l) = find_borders(str(ref_seq.seq), short_diff, diff_begin + pos - 1)
    unified_ref = ref_seq.seq[r:l]
    if ref_len > alt_len:
        unified_alt = unified_ref[0:len(unified_ref) - len(diff)]
    else:
        unified_alt = unified_ref + diff
    pos = r + 1
#    print "Variant is: ref - " + unified_ref + ", alt - " + unified_alt + ", pos - " + str(pos)
    return (pos, unified_ref, unified_alt, ("+" + diff if alt_len > ref_len else "-" + diff))

def print_results(files, positions, out_f):
    out_f.write("Position;unified ref;unified alt;difference;")    
    for f in files:
        out_f.write(f + ";;;")
    out_f.write("\n")
    out_f.write(";;;;")
    for f in files:
        out_f.write("orig pos;ref;alt;")
    out_f.write("\n")

    for pos in positions:
#        print(str(pos.pos) + ";" + pos.unified_ref + ";" + pos.unified_alt + ";")
        out_f.write(str(str(pos.pos) + ";" + pos.unified_ref + ";" + pos.unified_alt + ";" + pos.diff + ";"))
        for f in files:
            var = next(([p, ref, alt] for (p, ref, alt, fn) in pos.variants if fn == f), ["", "", ""])
            out_f.write(str(var[0]) + ";" + str(var[1]) + ";" + str(var[2]) + ";")
        out_f.write("\n")

def convert_indels(ref, alt):
    if alt[0] == '-':
        return ref + alt[1:], ref
    if alt[0] == '+':
        return ref, ref + alt[1:]
    return ref, alt

def merge(pos, next_pos, shift):
    overlap_ref = pos.unified_ref[shift:]
    overlap_alt = pos.unified_alt[shift:]
    if not overlap_ref or not overlap_alt:
#        print "Failed to merge because of no overlap"
        return (None, None)
    if pos.unified_ref == next_pos.unified_ref:
#        print "Failed to merge because reference sequences are same"
        return (None, None)
    if not (overlap_ref.startswith(next_pos.unified_ref) or next_pos.unified_ref.startswith(overlap_ref)):
#        print "Failed to merge references " + pos.unified_ref + " and " + next_pos.unified_ref
        return (None, None)
    if not (overlap_alt.startswith(next_pos.unified_alt) or next_pos.unified_alt.startswith(overlap_alt)):
#        print "Failed to merge alternatives " + pos.unified_alt + " and " + next_pos.unified_alt
        return (None, None)
    merged_ref = pos.unified_ref[:shift] + (next_pos.unified_ref if len(next_pos.unified_ref) > len(overlap_ref) else overlap_ref)
    merged_alt = pos.unified_alt[:shift] + (next_pos.unified_alt if len(next_pos.unified_alt) > len(overlap_alt) else overlap_alt)
#    print "Old refs are " + pos.unified_ref + " and " + next_pos.unified_ref + ", new - " + merged_ref
#    print "Old alts are " + pos.unified_alt + " and " + next_pos.unified_alt + ", new - " + merged_alt
    print "Position " + str(pos.pos) + " and " + str(next_pos.pos) + " merged successfully"
    return (merged_ref, merged_alt)

if __name__ == "__main__":
    dirname = os.path.dirname(sys.argv[0])
    parser = ArgumentParser(description='Test mappers and snp callers')
    parser.add_argument("-l", action="store", dest="input", help="list of files to process", required=True)
    parser.add_argument("-r", action="store", dest="ref", help="reference file", required=True)

    args = parser.parse_args()

    if not os.path.exists(args.input):
        sys.stderr.write("Directory does not exists")
        sys.exit(-1)

    positions = []
    files = []

    fin = open(args.input)

    ref_dict = SeqIO.index(args.ref, "fasta")
    for f in fin:
        files.append(f.strip())
        if not os.path.exists(files[-1]):
            sys.stderr.write("File " + files[-1] + " not exist\n")
            sys.exit(-1)
        f_curr = open(files[-1])
        for l in f_curr:
            if l[0] == "#":
                continue
            fields = l.strip().split('\t')
            if "snver" in files[-1]: # dirty hack for snver
                fields[1] = int(fields[1]) - 1
            if is_indel(fields):
                ref, alt = convert_indels(fields[3], fields[4])
                pos, unified_ref, unified_alt, diff = unify(int(fields[1]), ref, alt, files[-1], ref_dict[fields[0]])
                if not pos:
                    continue
                position = next((i for i in positions if i.pos == pos), None)
                if not position:
                    positions.append(Position(pos, fields[0], unified_ref, unified_alt, diff))
                    position = positions[-1]
                position.add_variant((int(fields[1]), ref, alt, files[-1]))

    positions = sorted(positions, key=attrgetter('chr', 'pos'))

    to_del = []
    for i, position in enumerate(positions):
        if i in to_del:
            continue
        for j, next_pos in enumerate(positions[i + 1:]):
            if position.pos + len(position.unified_ref) < next_pos.pos:
                break
#            print "Trying to merge " + str(position.pos) + " - " + position.unified_ref + " and " + str(next_pos.pos) + " - " + next_pos.unified_ref
            (merged_ref, merged_alt) = merge(position, next_pos, next_pos.pos - position.pos)
            if not merged_ref:
                break
            position.unified_ref = merged_ref
            position.unified_alt = merged_alt
            if len(position.unified_ref) > len(position.unified_alt):
                position.diff = '-' + position.unified_ref[len(position.unified_alt):]
            else:
                position.diff = '+' + position.unified_alt[len(position.unified_ref):]
            position.variants += next_pos.variants
            to_del.append(j + i + 1)
    for i in to_del[::-1]:
        del positions[i]

    out_f = open("out.csv", "w")

    print_results(files, positions, out_f)

    shared = [ [0] * len(files) for i in xrange(len(files))]
    complete = [0] * len(files)

    for position in positions:
        curr = set([files.index(fn) for (pos, ref, alt, fn) in position.variants])
#        print "Position is " + str(position.pos)
        for i in curr:
            complete[i] += 1
        if len(curr) == 1:
            continue
        for i in list(curr)[:-1]:
            for j in list(curr)[i + 1:]:
#                print str(i) + " " + str(j)
#                print "Share " + files[i] + " and " + files[j]
                shared[i][j] += 1
                shared[j][i] += 1

    out_f.close()
    out_f = open("out_stats.csv", "w")

    out_f.write(";")
    for f in files:
        out_f.write(f + ';')
    out_f.write("\n")

    out_f.write("# of variants;")
    for i in xrange(len(files)):
        out_f.write(str(complete[i]) + ';')

    out_f.write('\n\nShared variants;')
    for f in files:
        out_f.write(f + ';')
    out_f.write("\n")
    for i in xrange(len(files)):
        out_f.write(files[i] + ';')
        for j in xrange(len(files)):
            out_f.write(str(shared[i][j]) + ';')
        out_f.write("\n")
    
