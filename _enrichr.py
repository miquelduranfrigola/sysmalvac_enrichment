import sys, os
import math
import argparse
import uuid
import subprocess

# Variables

SCRATCH = "/aloy/scratch/mduran/sysmalvac/"

PATH = os.path.dirname(os.path.realpath(__file__))

CMD = "/usr/local/bin/singularity exec %s/singularity/sysmalvac.simg python" % PATH

SCRIPT = PATH + "/enrichr.py"

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


# Prepare input files

os.chdir(SCRATCH)

TMPFILE = SCRATCH + "/" + str(uuid.uuid4())

with open(TMPFILE, "w") as f:
    S = 0
    for x in filenames:
        f.write("%s\n" % (x))
        S += 1
t = math.ceil(float(S)/granularity)

subprocess.Popen("./setupArrayJob.py -x -N malenrichr -l %s -t %d %s %s" % (TMPFILE, t, CMD, SCRIPT), shell = True).wait()
s = subprocess.check_output("qsub -q all.q,fast.q job-malenrichr.sh", shell = True)
job_id = int(s.split("job-array ")[1].split(".")[0])