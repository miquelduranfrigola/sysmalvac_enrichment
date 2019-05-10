# See libraries at http://amp.pharm.mssm.edu/Enrichr/#stats

# Imports

import sys, os
import gseapy
import csv
import pandas as pd
from scipy.stats import fisher_exact
import collections
import numpy as np
import uuid
import shutil

# Variables

PATH  = os.path.dirname(os.path.realpath(__file__))

alpha = 0.05
top   = 100

# Make directories

filename = sys.argv[1]

OUTPATH = PATH + "/sysmalvac_results/" + filename + "/enrichr/"
if os.path.exists(OUTPATH):
    shutil.rmtree(OUTPATH)
os.mkdir(OUTPATH)

RNKPATH = PATH + "/data/diff_expr/genes_" + filename + ".rnk"

# Gene sets

gene_sets = [
"KEGG_2016",
"BioCarta_2016",
"WikiPathways_2016",
"Panther_2016"
]


# Genes

rnk = pd.read_csv(RNKPATH, delimiter = "\t", header = None, names = ["gene", "fc", "pval", "adj_pval"]).sort_values("fc", ascending = False)

def enrichr_test(direction = "pos", significance = True):

    genes = np.array(rnk["gene"])
    pvals = np.array(rnk["pval"])
    fcs   = np.array(rnk["fc"])

    if direction == "pos":
        if significance:
            hits = set(genes[np.logical_and(fcs > 0, pvals < alpha)])
        else:
            hits = set(genes[:top])
    else:
        if significance:
            hits = set(genes[np.logical_and(fcs < 0, pvals < alpha)])
        else:
            hits = set(genes[-top:])

    if significance:
        sur = "sign05"
    else:
        sur = "top100"

    outpath = OUTPATH + "/enrichr_%s_%s/" % (direction, sur)
    if not os.path.exists(outpath): os.mkdir(outpath)

    enr = gseapy.enrichr(gene_list=list(hits),
                         description="%s_%s_%s" % (filename, direction, sur),
                         gene_sets=gene_sets,
                         outdir=outpath,
                         cutoff=0.5 # only used for plotting.
                         )

enrichr_test("pos", True)
enrichr_test("neg", True)
enrichr_test("pos", False)
enrichr_test("neg", False)