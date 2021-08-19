import argparse
import tempfile
from os.path import join
import os
import asterid as ad
from random import sample
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str,
                        help="Input tree list file", required=True)
parser.add_argument("-o", "--output", type=str, required=True)
args = parser.parse_args()
with open(args.input) as files: genes = files.readlines()

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
        trees.append(ad.fastme_balme(ts, D, 1, 1))
print(len(trees))
def consensus(trees, minfreq=0.5):
    import dendropy
    res = dendropy.TreeList()
    for treenewick in trees:
        res.read(data=treenewick, schema="newick", rooting='force-unrooted')
    con = res.consensus(min_freq = minfreq)
    con.is_rooted = False
    return con.as_string(schema="newick")
C = consensus(trees)
with open(args.output, "w+") as fh:
    fh.write(C[5:])