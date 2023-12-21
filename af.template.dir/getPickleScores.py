import sys
import pandas as pd
import glob
import os
import os.path

outputDir = sys.argv[1]
outputFile = os.path.join(outputDir, "scores.csv")

pickleFiles = glob.glob(os.path.join(outputDir, "result_*.pkl"))
print ("\n".join(pickleFiles))

#fields = ("ptm", "iptm")
allDats = [pd.read_pickle(pf) for pf in pickleFiles]
ptms = [dat["ptm"] for dat in allDats]
iptms = [dat["iptm"] for dat in allDats]


with open(outputFile, "w") as fp:
    # header
    fp.write (f"model,ptm,iptm\n")
    # data 
    for pf,ptm,iptm in zip([os.path.basename(pf) for pf in pickleFiles], ptms, iptms):
        fp.write(f"{pf},{ptm},{iptm}\n")
        


