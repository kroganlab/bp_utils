import sys
import glob
import os.path
import json
import pandas as pd

# sys.argv should either be a single directory or
# if first arg is --files then a space separated list of pkl files

outputDir = sys.argv[1]
if outputDir == "--files":
  pickleFiles = sys.argv[2:]
else:
  pickleFiles = glob.glob(os.path.join(outputDir, "result_*.pkl"))


print ("\n".join(pickleFiles))

for picklePath in pickleFiles:
  print (picklePath)
  outputFile = os.path.splitext(picklePath)[0] + ".PAE.json"
  dat = pd.read_pickle(picklePath)
  with open(outputFile, "w") as fp:
    json.dump(dat["predicted_aligned_error"].tolist(), fp)


        


