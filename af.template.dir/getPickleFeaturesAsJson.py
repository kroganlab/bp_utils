import sys
import glob
import os.path
import json
import pandas as pd


# Opens files like features.pkl and outputs a subset of feature objects as json

# sys.argv should either be a single directory or
# if first arg is --files then a space separated list of features.pkl files

outputDir = sys.argv[1]
if outputDir == "--files":
  pickleFiles = sys.argv[2:]
else:
  pickleFiles = glob.glob(os.path.join(outputDir, "features.pkl"))


print ("\n".join(pickleFiles))

# ['aatype', 'residue_index', 'seq_length', 'msa', 'num_alignments', 'template_aatype', 'template_all_atom_mask', 'template_all_atom_positions', 'asym_id', 'sym_id', 'entity_id', 'deletion_matrix', 'deletion_mean', 'all_atom_mask', 'all_atom_positions', 'assembly_num_chains', 'entity_mask', 'num_templates', 'cluster_bias_mask', 'bert_mask', 'seq_mask', 'msa_mask']
featuresOI = ("aatype", "residue_index", "msa", "template_aatype", "entity_id")

for picklePath in pickleFiles:
  print (picklePath)
  features = pd.read_pickle(picklePath)
  for featureName in featuresOI:
    try:
      dat = features[featureName]
    except KeyError:
      print(f"Feature {featureName} not found in {picklePath}")
      next
     
    outputFile = os.path.splitext(picklePath)[0] + f".{featureName}.json"
    with open(outputFile, "w") as fp:
      json.dump(dat.tolist(), fp)
        
  


        


