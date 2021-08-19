import argparse
import tempfile
from os.path import join
import os
import asterid as ad
from random import sample
import treeswift as tsf
from copy import copy
import nj
ASTRALMPATH = "/home/baqiaol2/scratch/atsume/methods/FASTRAL/FASTRAL/ASTRAL-modified/astral.5.7.3.modified.jar"

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str,
                        help="Input tree list file", required=True)
parser.add_argument("-j", "--constraint", type=str,
                        help="Input constraint tree", required=True)
parser.add_argument("-o", "--output", type=str, required=True)
args = parser.parse_args()
with open(args.input) as files: genes = files.readlines()
with open(args.constraint) as j: constraint = tsf.read_tree_newick(j.read())

k = len(genes)
samplelen = []
samplesize = []
samplelen.append(1)
samplesize.append(k)
samplelen.append(10)
samplesize.append(round(k * 0.5))
samplelen.append(20)
samplesize.append(round(k * 0.25))
samplelen.append(20)
samplesize.append(round(k * 0.1))
trees = []
for l, s in zip(samplelen, samplesize):
    for i in range(l):
        tees = sample(genes, s)
        ts = ad.get_ts(tees)
        D = ad.mk_distance_matrix(ts, tees)
        trees.append(nj.treeresolve_lua(copy(constraint), ts, D).newick())
        # trees.append(ad.fastme_balme(ts, D, 1, 1))

ctreepath = args.output + ".j"
with open(ctreepath, "w+") as fh:
    for t in trees:
        fh.write(t)
        fh.write("\n")
import os
os.system(
    f"java -jar {ASTRALMPATH} -i {args.input} -s {ctreepath} -o {args.output}")

# print(len(trees))
# def consensus(trees, minfreq=0.5):
#     import dendropy
#     res = dendropy.TreeList()
#     for treenewick in trees:
#         res.read(data=treenewick, schema="newick", rooting='force-unrooted')
#     con = res.consensus(min_freq = minfreq)
#     con.is_rooted = False
#     return con.as_string(schema="newick")
# C = consensus(trees)
# with open(args.output, "w+") as fh:
#     fh.write(C[5:])