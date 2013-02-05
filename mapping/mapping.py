from argparse import ArgumentParser
import os
import re

def get_mappers_list(technology, is_paired):
    files = os.listdir(".")
    prog_pattern = re.compile("[^_]+\.sh")
    if (technology == "454"):
        if (is_paired):
            prog_pattern = re.compile(".+_454_paired\.sh")
        else:
            prog_pattern = re.compile(".+_454\.sh")
    if (technology == "Solid"):
        if (is_paired):
            prog_pattern = re.compile(".+_solid_paired\.sh")
        else:
            prog_pattern = re.compile(".+_solid\.sh")
    mappers = []
    for f in files:
        if (prog_pattern.match(f)):
            mappers.append(f)
    return mappers

parser = ArgumentParser(description='Test mappers')
parser.add_argument("-i", action="store", dest="reads")
parser.add_argument("-1", action="store", dest="reads1")
parser.add_argument("-2", action="store", dest="reads2")
parser.add_argument("-r", "--ref", action="store", dest="ref", required=True)
parser.add_argument("-t", "--threads", action="store", dest="threads", default=4)
parser.add_argument("-p", "--progs", action="store", dest="prog_list", default="*")
parser.add_argument("--hash", action="store", dest="hash_sz", default=10)
parser.add_argument("-o", "--out", action="store", dest="out", required=True)
parser.add_argument("-a", action="store", dest="additional")
parser.add_argument("--tech", action="store", dest="technology", default="Illumina", choices=['Illumina', '454', 'Solid'])
parser.add_argument("--test", action="store_true", dest="test")

args = parser.parse_args()
print args

base_name_pattern = re.compile("(\w+)(\.\w+)$")
ref_name = base_name_pattern.search(args.ref).group(1)

reads1 = (args.reads) ? args.reads : args.reads1
reads1_base = base_name_pattern.search(reads1).group(1)
if (not reads1):
    sys.stderr.write("No reads file is specified")
    sys.exit(-1)

is_paired = False
reads2 = ""
reads2_base = ""
if (args.reads2):
    reads2 = args.reads2
    reads2_base = base_name_pattern.search(reads2).group(1)
    is_paired = True

if (args.out and not os.path.exists(args.out)):
    os.makedirs(args.out)

progs = []

if (args.prog_list == "*"):
    progs = get_mappers_list(args.technology, is_paired)
else:
    progs = args.prog_list.split(" ")
    prog_pattern = re.compile(".*\.sh")
    for i in range(progs):
        if (not prog_pattern.match(progs[i])):
            progs[i] += ".sh"

for prog in progs:
    dir_path = args.out + "/" + prog
    if (not os.path.exists(dir_path)):
        os.makedirs(dir_path)

    if (is_paired):
        os.system("./" + prog + " " + args.ref + " " + ref_base + " " 
                + reads1 + " " + reads1_base + " "
                + reads2 + " " + reads2_base + " "
                + dir_path + " " + args.threads + " "
                + args.hash_sz + " " 
                + ((args.additional) ? args.additional : "") 
                + " | tee -a " + dir_path + "/mapping.log")
    else:
        os.system("./" + prog + " " + args.ref + " " + ref_base + " " 
                + reads1 + " " + reads1_base + " "
                + dir_path + " " + args.threads + " "
                + args.hash_sz + " " 
                + ((args.additional) ? args.additional : "") 
                + " | tee -a " + dir_path + "/mapping.log")

    if 
