import sys
import glob
import os.path
import json
import pandas as pd


# Opens files like features.pkl and writes out the msa as single letter amino acid codes

# sys.argv should either be a single directory or
# if first arg is --files then a space separated list of features.pkl files

outputDir = sys.argv[1]
if outputDir == "--files":
  pickleFiles = sys.argv[2:]
else:
  pickleFiles = glob.glob(os.path.join(outputDir, "features.pkl"))


if len(pickleFiles) > 0:
  print ("\n".join(pickleFiles))
else:
  raise ValueError("no pickle files found in input")

restypes = [
  'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P',
  'S', 'T', 'W', 'Y', 'V',
  'X', '-']

featureName = "msa"

def convertRow(row):
  return ("".join([restypes[i] for i in row]))


for picklePath in pickleFiles:
  print (picklePath)
  features = pd.read_pickle(picklePath)
  try:
    dat = features[featureName]
    outputFile = os.path.splitext(picklePath)[0] + f".{featureName}.aln"
    print (outputFile)
    with open(outputFile, "w") as fp:
      for row in dat:
        fp.write (convertRow(row))
        fp.write ("\n")
  except KeyError:
    print(f"Feature {featureName} not found in {picklePath}")

     
        
  


        


