# Imports

import sys, os
import gseapy
import csv
import pandas as pd
from scipy.stats import fisher_exact
import collections
import numpy as np
import uuid

# Variables

PATH  = os.path.dirname(os.path.realpath(__file__))

min_size = 10

# Make directories

filename, geneset = sys.argv[1].split("---")

OUTPATH = PATH + "/sysmalvac_results_wgcna/" + filename + "/" + geneset + "/"

RNKPATH = PATH + "/data/wgcna_cor/genes_" + filename + ".rnk"
GENESETPATH = PATH + "/data/genesets/" + geneset + ".gmt"


# Read data

# Genes

rnk = pd.read_csv(RNKPATH, delimiter = "\t", header = None, names = ["gene", "fc"]).sort_values("fc", ascending = False)

# Read genesets

gs = collections.defaultdict(set)

with open(GENESETPATH, "r") as f:
    
    for r in csv.reader(f, delimiter = "\t"):
        gs[r[0]].update(r[2:])

# Universes

array_universe = set()
with open(PATH + "/data/probeset2genesymbol.tsv", "r") as f:
    for r in csv.reader(f, delimiter = "\t"):
        array_universe.update([r[1]])

samples_universe = set(rnk["gene"])

geneset_universe = set()
for k,v in gs.items():
    geneset_universe.update(v)

# GSEAPY

SCRATCH = "/aloy/scratch/mduran/sysmalvac/"

def gsea(rnk, full_universe):

    if full_universe:
        sur = "fulluniv"
        universe = array_universe.union(geneset_universe)
    else:
        sur = "reduniv"
        universe = samples_universe.intersection(geneset_universe)

    outpath = OUTPATH + "/gsea_%s/" % sur
    if not os.path.exists(outpath): os.mkdir(outpath)

    tmpfile = SCRATCH + "/" + str(uuid.uuid4())

    with open(tmpfile, "w") as f:

        genes = np.array(rnk["gene"])
        fcs   = np.array(rnk["fc"])

        for g_, fc_ in zip(genes, fcs):
            if g_ not in universe: continue
            f.write("%s\t%.3f\n" % (g_, fc_))

    res = gseapy.prerank(rnk = tmpfile, gene_sets = GENESETPATH, processes = 1, permutation_num = 10000, min_size = min_size, outdir = outpath, format = "pdf")

    with open(outpath + "/output.tsv", "w") as f:
        f.write("geneset\tsize\tmatched_size\tdirection\tes\tnes\tpval\tfdr\tleading_edge\n")
        for k,v in res.results.items():
            if v['nes'] < 0:
                direction = "neg"
                continue # I do not save negative enrichments (it is none-sense, as all correlations are positive)
            else:
                direction = "pos"
            f.write("%s\t%d\t%d\t%s\t%.3f\t%.3f\t%.3E\t%.3E\t%s\n" % (k, v["geneset_size"], v["matched_size"], direction, v["es"], v["nes"], v["pval"], v["fdr"], v["ledge_genes"]))


gsea(rnk, True)
gsea(rnk, False)
