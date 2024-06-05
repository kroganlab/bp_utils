import sys
import pandas as pd
import glob
import os
import os.path
import json


outputDir = sys.argv[1]
outputFile = os.path.join(outputDir, "scores.csv")

pickleFiles = glob.glob(os.path.join(outputDir, "result_*.pkl"))
print ("\n".join(pickleFiles))

#fields = ("ptm", "iptm")
allDats = [pd.read_pickle(pf) for pf in pickleFiles]
ptms = [dat["ptm"] for dat in allDats]
iptms = [dat["iptm"] for dat in allDats]


with open(outputFile, "a") as fp:
    # header
    fp.write (f"model,ptm,iptm\n")
    # data 
    for pf,ptm,iptm in zip([os.path.basename(pf) for pf in pickleFiles], ptms, iptms):
        fp.write(f"{pf},{ptm},{iptm}\n")
        


# Get other items from pkl files:
# available are:
# ['distogram', 'experimentally_resolved', 'masked_msa', 'num_recycles', 'predicted_aligned_error', 'predicted_lddt', 'structure_module', 'plddt', 'aligned_confidence_probs', 'max_predicted_aligned_error', 'ptm', 'iptm', 'ranking_confidence']
featuresOI = ('predicted_aligned_error',)
for pickleFile, dat in zip(pickleFiles,allDats):
  for featureName in featuresOI:
    try:
      singleDat = dat[featureName]
    except KeyError:
      print(f"Feature {featureName} not found in {pickleFile}")
      next
    outputFile = os.path.splitext(pickleFile)[0] + f".{featureName}.json"
    print (f"Writing {featureName} to {outputFile}")
    with open(outputFile, "w") as fp:
      json.dump(singleDat.tolist(), fp)


