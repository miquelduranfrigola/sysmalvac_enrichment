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

alpha = 0.05
top   = 100

min_size = 10

# Make directories

filename, geneset = sys.argv[1].split("---")

OUTPATH = PATH + "/sysmalvac_results/" + filename + "/" + geneset + "/"

RNKPATH = PATH + "/data/diff_expr/genes_" + filename + ".rnk"
GENESETPATH = PATH + "/data/genesets/" + geneset + ".gmt"


# Read data

# Genes

rnk = pd.read_csv(RNKPATH, delimiter = "\t", header = None, names = ["gene", "fc", "pval", "adj_pval"]).sort_values("fc", ascending = False)

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

# Fisher test

def fisher_test(rnk, full_universe = True, significance = True):

    if significance:
        sur1 = "sign05"
    else:
        sur1 = "top100"

    if full_universe:
        sur2 = "fulluniv"
        universe = array_universe.union(geneset_universe)
    else:
        sur2 = "reduniv"
        universe = samples_universe.intersection(geneset_universe)

    outpath = OUTPATH + "/fisher_%s_%s/" % (sur1, sur2)
    if not os.path.exists(outpath): os.mkdir(outpath)

    with open(outpath + "/output.tsv", "w") as f:

        f.write("geneset\tsize\tdirection\tA\tB\tC\tD\todds\tpval\tgenes\n")

        # Ensure they are in the Universe

        genes = []
        fcs   = []
        pvals = []
        for g_, fc_, p_ in zip(np.array(rnk["gene"]), np.array(rnk["fc"]), np.array(rnk["pval"])):
            if g_ not in universe: continue
            genes += [g_]
            fcs   += [fc_]
            pvals += [p_]
        genes = np.array(genes)
        fcs   = np.array(fcs)
        pvals = np.array(pvals)

        if significance:
            pos   = set(genes[np.logical_and(fcs > 0, pvals < alpha)])
            neg   = set(genes[np.logical_and(fcs < 0, pvals < alpha)])
        else:
            pos   = set(genes[:top])
            neg   = set(genes[-top:])


        def fisher(gs, hits, direction):

            for k,v in gs.items():

                v = v.intersection(universe)

                if len(v) < min_size: continue

                common = v.intersection(hits)
                A = len(common)
                B = len(hits) - A
                C = len(v) - A
                D = len(universe) - (A + B + C)

                try:
                    odds, pval = fisher_exact([[A, B], [C, D]], alternative = "greater")
                except:
                    odds, pval = -666, -666

                f.write("%s\t%d\t%s\t%d\t%d\t%d\t%d\t%.3f\t%.3E\t%s\n" % (k, len(v), direction, A, B, C, D, odds, pval, ";".join(common)))

        fisher(gs, pos, "pos")
        fisher(gs, neg, "neg")

    fd = pd.read_csv(outpath + "/output.tsv", delimiter = "\t")

    fd = fd.sort_values("odds", ascending = False)

    fd.to_csv(outpath + "/output.tsv", sep = "\t", index = False)

fisher_test(rnk, True, True)
fisher_test(rnk, True, False)
fisher_test(rnk, False, True)
fisher_test(rnk, False, False)

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
            else:
                direction = "pos"
            f.write("%s\t%d\t%d\t%s\t%.3f\t%.3f\t%.3E\t%.3E\t%s\n" % (k, v["geneset_size"], v["matched_size"], direction, v["es"], v["nes"], v["pval"], v["fdr"], v["ledge_genes"]))


gsea(rnk, True)
gsea(rnk, False)