# Imports

import sys, os
import math
import argparse
import uuid
import subprocess

# Variables

SCRATCH = "/aloy/scratch/mduran/sysmalvac/"

PATH = os.path.dirname(os.path.realpath(__file__))

CMD = "/usr/local/bin/singularity exec %s/singularity/sysmalvac.simg python" % PATH

SCRIPT = PATH + "/enrichment.py"

granularity = 1

# Create folder structure

filenames = [
    "CSP_I7_C1",
    "CSP_RTS_CONTROL",
    "iRBC_I7_C1",
    "iRBC_RTS_CONTROL",
    "malaria_CSP_CONTROL",
    "malaria_CSP_RTS",
    "malaria_iRBC_CONTROL",
    "malaria_iRBC_RTS",
    "protected_CSP_C1",
    "protected_CSP_I7",
    "protected_iRBC_C1",
    "protected_iRBC_I7"
]

genesets = [
    "btm",
    "modules",
    "biocarta",
    "kegg",
    "reactome"    
]

for filename in filenames:
    fld1 = PATH + "/sysmalvac_results/" + filename
    if not os.path.exists(fld1): os.mkdir(fld1)
    for geneset in genesets:
        fld2 = fld1 + "/" + geneset
        if not os.path.exists(fld2): os.mkdir(fld2)

# Prepare input files

os.chdir(SCRATCH)

TMPFILE = SCRATCH + "/" + str(uuid.uuid4())

with open(TMPFILE, "w") as f:
    S = 0
    for x in filenames:
        for y in genesets:
            f.write("%s---%s\n" % (x,y))
            S += 1
t = math.ceil(float(S)/granularity)

subprocess.Popen("./setupArrayJob.py -x -N malenrich -l %s -t %d %s %s" % (TMPFILE, t, CMD, SCRIPT), shell = True).wait()
s = subprocess.check_output("qsub -q all.q,fast.q job-malenrich.sh", shell = True)
job_id = int(s.split("job-array ")[1].split(".")[0])
