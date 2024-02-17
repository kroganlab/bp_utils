import os
import sys
import re
from Bio import AlignIO


## Copied from alphafold code
# Sequences coming from UniProtKB database come in the
# `db|UniqueIdentifier|EntryName` format, e.g. `tr|A0A146SKV9|A0A146SKV9_FUNHE`
# or `sp|P0C2L1|A3X1_LOXLA` (for TREMBL/Swiss-Prot respectively).
_UNIPROT_PATTERN = re.compile(
    r"""
    ^
    # UniProtKB/TrEMBL or UniProtKB/Swiss-Prot
    (?:tr|sp)
    \|
    # A primary accession number of the UniProtKB entry.
    (?P<AccessionIdentifier>[A-Za-z0-9]{6,10})
    # Occasionally there is a _0 or _1 isoform suffix, which we ignore.
    (?:_\d)?
    \|
    # TREMBL repeats the accession ID here. Swiss-Prot has a mnemonic
    # protein ID code.
    (?:[A-Za-z0-9]+)
    _
    # A mnemonic species identification code.
    (?P<SpeciesIdentifier>([A-Za-z0-9]){1,5})
    # Small BFD uses a final value after an underscore, which we ignore.
    (?:_\d+)?
    #$
    """,
    re.VERBOSE)



def RemoveSpeciesFromUniprotMSA():
  a = AlignIO.read("uniprot_hits.sto", "stockholm")
  ids = [rec.id for rec in a]
  species = []
  for i,ID in enumerate (ids):
    m = re.search(_UNIPROT_PATTERN, ID)
    if m:
      species.append (m.group('SpeciesIdentifier'))
      newSpecies = ("ZZ%s" % m.group('SpeciesIdentifier'))[0:5] # HUMAN becomes ZZHUM
      newID = ID.replace(m.group('SpeciesIdentifier'), newSpecies)
      a[i].id = newID
    else:
      species.append(None)

  # confirm here that no ZZ... species exist
  # https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/speclist.txt
  newSpecies = [ ("ZZ%s" % s)[0:5] for s in species]

  # make sure I havne't just created a name that matches an existing species (in current MSA at least).
  assert len(set(species).intersection(set(newSpecies))) == 0

  with open("uniprot_hits_notSpecies.sto", "w") as fp:
    #fp.write (a.format("stockholm"))
    fp.write (format(a, "stockholm"))



if __name__ == "__main__":
  try:
    os.chdir(sys.argv[1])
    RemoveSpeciesFromUniprotMSA()
  except IndexError as e:
    print (e)
    print ("Usage: RemoveSpeciesToBlockMSAPairing.py output/DYRK1A__DCAF7/msas/A/")
    print ("   only argument is path that contains uniprot_hits.sto")
    print ("   recommended to only do one sequence per dimer")
    
