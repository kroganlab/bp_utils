#!/usr/bin/python3
#
#$ -S /usr/bin/python3
#$ -q gpu.q
#$ -N alphafold 
#$ -cwd
###$ -l h_rt=24:00:00
#$ -l h_rt=48:00:00
#$ -l mem_free=60G
#$ -l scratch=50G
#$ -l compute_cap=80,gpu_mem=40G
#
# Compute cap for A100 GPU is 8.0 (40 or 80 GB), for A40 GPU is 8.6 (48 GB).
#
# Adapted from alphafold/docker/run_alphafold.py script.
# Original version runs AlphaFold using a docker image.
# This adapted version uses a singularity image with defaults
# set for the UCSF Wynton cluster.
#

"""Singularity launch script for Alphafold."""

def parse_args():
  import argparse

  parser = argparse.ArgumentParser(description='Run AlphaFold structure prediction using singularity image.')
  
  # BP added arguments to manage multiple jobs running in a task array
  parser.add_argument(
    '--job_id', type = int, required = True,
    help = 'job id to look up sequence names in AlphaFoldJobList')

  parser.add_argument(
    '--master_fasta', default = "masterFasta.fasta",
    help = "fasta file with sequence info relevant to names in AlphaFoldJobList")
    
  parser.add_argument(
    '--jobTable', default= "AlphaFoldJobList.csv",
    help = "path to csv formatte file with columns ID, seq1.name, seq2.name")

  parser.add_argument(
    '--alignmentRepo', default  =  "/wynton/group/gladstone/users/bpolacco/AF_msas",
    help = "path to directory containing previously generated MSAs")
    
  parser.add_argument(
    '--run_alpha_fold', default = True, type = str_to_bool,
    help = "Helpful to disable and get only the setup, score extraction, msa archive steps")
    
  parser.add_argument(
    '--setup_job', type = str_to_bool, default = True,
    help = "Disable to skip the setup steps (msa copying, pair fasta writing)")
    
  parser.add_argument(
    '--prevent_alphafold_output', type = str_to_bool, default = False,
    help = "attempts to lock files in output directory to kill the alphafold job early (good for MSA generation on non-GPU nodes)")
    
  parser.add_argument(
    '--postprocess_job', type = str_to_bool, default = True,
    help = "Disable to skip the job postprocessing (score extraction from pickle files and copying MSA)")

  parser.add_argument(
          '--check_if_completed', type = str_to_bool, default = True,
          help = "Check if the run is completed already. If so, don't rerun.")

  # /END BP added arguments

  parser.add_argument(
    '--fasta_paths', required=False, default = "notUsed",
    help='BP: NOT USED in this version. '
    'Paths to FASTA files, each containing a prediction '
    'target that will be folded one after another. If a FASTA file contains '
    'multiple sequences, then it will be folded as a multimer. Paths should be '
    'separated by commas. All FASTA paths must have a unique basename as the '
    'basename is used to name the output directories for each prediction.')

  parser.add_argument(
    '--use_gpu', type=str_to_bool, default=True,
    help='Enable NVIDIA runtime to run with GPUs.')

  import os
  parser.add_argument(
    '--gpu_devices', default=os.environ.get('SGE_GPU', '0'),
    help='Comma separated list GPU identifiers to set environment variable CUDA_VISIBLE_DEVICES.')

  parser.add_argument(
    '--run_relax', type=str_to_bool, default=True,
    help='Whether to do OpenMM energy minimization of each predicted structure.')

  parser.add_argument(
    '--use_gpu_relax', type=str_to_bool, default=True,
    help='Whether to do OpenMM energy minimization using GPU.')

  parser.add_argument(
    '--output_dir', default='output',
    help='Path to a directory that will store the results.')

  parser.add_argument(
    '--data_dir', default='/wynton/home/krogan/bpolacco/gladstoneHome/af.data/databases/af_C14_v230',
    help='Path to directory with supporting data: AlphaFold parameters and genetic '
    'and template databases. Set to the target of download_all_databases.sh.')

  parser.add_argument(
    '--mount_data_dir', default='/wynton/group/databases',
    help='Path to directory where databases reside. On UCSF Wynton '
    'some of the databases are symbolic links to various locations in this directory '
    'and singularity needs to mount this directory to see them.')

  parser.add_argument(
    '--singularity_image_path', default='/wynton/home/ferrin/goddard/alphafold_singularity/alphafold231.sif',
    help='Path to the AlphaFold singularity image.')

  parser.add_argument(
    '--max_template_date', default='2100-01-01',
    help='Maximum template release date to consider (ISO-8601 format: YYYY-MM-DD). '
    'Important if folding historical test sets.')

  parser.add_argument(
    '--db_preset', default='full_dbs', choices=['full_dbs', 'reduced_dbs'],
    help='Choose preset MSA database configuration - smaller genetic database '
    'config (reduced_dbs) or full genetic database config (full_dbs)')

  parser.add_argument(
    '--model_preset', default='monomer_ptm',
    choices=['monomer', 'monomer_casp14', 'monomer_ptm', 'multimer'],
    help='Choose preset model configuration - the monomer model, the monomer model '
    'with extra ensembling, monomer model with pTM head, or multimer model')

  parser.add_argument(
      '--num_multimer_predictions_per_model', default=1,
      help='How many predictions (each with a different random seed) will be '
      'generated per model. E.g. if this is 2 and there are 5 '
      'models then there will be 10 predictions per input. '
      'Note: this FLAG only applies if model_preset=multimer')

  parser.add_argument(
    '--benchmark', default=False,
    help='Run multiple JAX model evaluations to obtain a timing that excludes the '
    'compilation time, which should be more indicative of the time required '
    'for inferencing many proteins.')

  parser.add_argument(
    '--use_precomputed_msas', default=True,
    help='Whether to read MSAs that have been written to disk. WARNING: This will '
    'not check if the sequence, database or configuration have changed.')

  args = parser.parse_args()
  return args

def str_to_bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        import argparse
        raise argparse.ArgumentTypeError('Boolean value expected.')

def main():

  args = parse_args()
  
  ############################
  fastaPath = BP_setupJob(args)
  args.fasta_paths = fastaPath
  ###########################
  

  # You can individually override the following paths if you have placed the
  # data in locations other than the parser.data_dir.

  # Path to the Uniref90 database for use by JackHMMER.
  import os.path
  uniref90_database_path = os.path.join(
      args.data_dir, 'uniref90', 'uniref90.fasta')

  # Path to the Uniprot database for use by JackHMMER.
  uniprot_database_path = os.path.join(
      args.data_dir, 'uniprot', 'uniprot.fasta')

  # Path to the MGnify database for use by JackHMMER.
  mgnify_database_path = os.path.join(
      args.data_dir, 'mgnify', 'mgy_clusters_2022_05.fa')

  # Path to the BFD database for use by HHblits.
  bfd_database_path = os.path.join(
      args.data_dir, 'bfd',
      'bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt')

  # Path to the Small BFD database for use by JackHMMER.
  small_bfd_database_path = os.path.join(
      args.data_dir, 'small_bfd', 'bfd-first_non_consensus_sequences.fasta')

  # Path to the UniRef30 database for use by HHblits.
  uniref30_database_path = os.path.join(
      args.data_dir, 'uniref30', 'UniRef30_2023_02')

  # Path to the PDB70 database for use by HHsearch.
  pdb70_database_path = os.path.join(args.data_dir, 'pdb70', 'pdb70')

  # Path to the PDB seqres database for use by hmmsearch.
  pdb_seqres_database_path = os.path.join(
      args.data_dir, 'pdb_seqres', 'pdb_seqres.txt')

  # Path to a directory with template mmCIF structures, each named <pdb_id>.cif.
  template_mmcif_dir = os.path.join(args.data_dir, 'pdb_mmcif', 'mmcif_files')

  # Path to a file mapping obsolete PDB IDs to their replacements.
  obsolete_pdbs_path = os.path.join(args.data_dir, 'pdb_mmcif', 'obsolete.dat')

  mounts = []
  command_args = []

  # FASTA paths
  command_args.append(f'--fasta_paths={args.fasta_paths}')

  database_paths = [
      ('uniref90_database_path', uniref90_database_path),
      ('mgnify_database_path', mgnify_database_path),
      ('data_dir', args.data_dir),
      ('template_mmcif_dir', template_mmcif_dir),
      ('obsolete_pdbs_path', obsolete_pdbs_path),
  ]

  if args.model_preset == 'multimer':
    database_paths.append(('uniprot_database_path', uniprot_database_path))
    database_paths.append(('pdb_seqres_database_path',
                           pdb_seqres_database_path))
  else:
    database_paths.append(('pdb70_database_path', pdb70_database_path))

  if args.db_preset == 'reduced_dbs':
    database_paths.append(('small_bfd_database_path', small_bfd_database_path))
  else:
    database_paths.append(('uniref30_database_path', uniref30_database_path))
    database_paths.append(('bfd_database_path', bfd_database_path))

  for name, path in database_paths:
    if path:
      command_args.append(f'--{name}={path}')

  command_args.extend([
      f'--output_dir={args.output_dir}',
      f'--max_template_date={args.max_template_date}',
      f'--db_preset={args.db_preset}',
      f'--model_preset={args.model_preset}',
      f'--num_multimer_predictions_per_model={args.num_multimer_predictions_per_model}',
      f'--run_relax={args.run_relax}',
      f'--use_gpu_relax={args.use_gpu_relax}',
      f'--benchmark={args.benchmark}',
      f'--use_precomputed_msas={args.use_precomputed_msas}',
      '--logtostderr',
  ])

  env_vars = {
          'CUDA_VISIBLE_DEVICES': args.gpu_devices,
          'NVIDIA_VISIBLE_DEVICES': args.gpu_devices,
          # The following flags allow us to make predictions on proteins that
          # would typically be too long to fit into GPU memory.
          'TF_FORCE_UNIFIED_MEMORY': '1',
          'XLA_PYTHON_CLIENT_MEM_FRACTION': '4.0',
  }
  env_vals = ','.join('%s=%s' % (key,value) for key,value in env_vars.items())

  # AlphaFold uses Python tempfile which uses TMPDIR env variable
  # which is /scratch/job-id-string on wynton.  Otherwise Python will use /tmp
  # which is only 4-8 GB on wynton and will cause write errors on large sequences.
  import os
  tempdir = os.environ.get('TMPDIR', '/scratch')

  singularity_args = ['singularity',
          'run',
          '--nv',  # Use Nvidia container library to use CUDA
          '-B "%s"' % args.mount_data_dir,    # Mount AlphaFold databases
          '-B "/wynton/group/gladstone/users/bpolacco/af.data/databases"',
          '-B "%s"' % os.getcwd(),	# Mount current directory for sequence
          '-B "%s"' % tempdir,		# Mount scratch directory
          '--env %s' % env_vals, 
          args.singularity_image_path
        ] + command_args
  cmd = ' '.join(singularity_args)
  
  print("\n# AF run command:")
  print (cmd)
  print()

  from subprocess import run
  import sys
  if args.run_alpha_fold and BP_notYetCompleted(args):
    run('module load cuda/11.0 ; %s' % cmd,
        stdout = sys.stdout, stderr = sys.stderr,
        shell = True,  # module command is a csh alias on Wynton
        executable = '/bin/csh',
        check = True)
  else:
    if not args.run_alpha_fold:
      print ("Not running alphafold as requested by run_alpha_fold=False")
    else:
      print ("Alphafold appears to be completed already, not running")
  
  ########### clean up tasks ###########
  if args.prevent_alphafold_output:
    alphaFoldLockFiles( os.path.split(args.fasta_paths)[0], lock = False)

  if args.postprocess_job:
    # scores
    BP_processScores(args)  
    # archive the MSAS
    BP_archiveMSAs(args)
  ######################################

######################### BP Inserted functions ###############
import sys
import os.path
import collections
import hashlib
import shutil

# helper functions
def nameFromHeader (header):
  pipeParts = header[1:].split("|")
  if(pipeParts[0] in ("sp", "tr", "up")):
    return(pipeParts[1])
  else:
    return(pipeParts[0])

def read_fasta(path):
  allSeqs = collections.defaultdict(str)
  namesInOrder = []
  with open(path) as fp:
    currentName = ""
    for line in fp.readlines():
      line = line.strip()
      if line[0] == ">":
        currentName = nameFromHeader(line)
        if len(namesInOrder) == 0 or namesInOrder[-1] != currentName:
          namesInOrder.append (currentName)
      else:
        assert currentName != "", "Unexpected format, empty sequence name or sequence data before first >"
        allSeqs[currentName] = allSeqs[currentName] + line
  return allSeqs, namesInOrder



def getMSADirectoryForSequence(sequence, alignmentRepo):
  seqHash = hashlib.md5(sequence.upper().encode('utf-8')).hexdigest()
  subDir = seqHash[0:2]
  dir = os.path.join(alignmentRepo, subDir, seqHash)

  # if directory exists, raise an error if seqeunces mismatch
  try:
    fastaPaths = []
    fastaPaths = [p for p in os.listdir(dir) if p.split(".")[-1] == "fasta"]
  except:
    pass
  if len(fastaPaths) > 0:
    print (fastaPaths)
    for p in fastaPaths:
        for name,seq in read_fasta(os.path.join(dir,p))[0].items():
            if seq != sequence:
                raise Exception("Sequence Mismatch Error", f"{dir}:{sequence}")
    
  return dir
  

def alphaFoldRunOutputDirectory(seq1, seq2, outputPath):
  return os.path.join(outputPath, f"{seq1}__{seq2}")


def alphaFoldLockFiles (outDir, lock = True):
  for fn in ("relaxed_model_1_multimer_v3_pred_0.pdb", "result_model_1_multimer_v3_pred_0.pkl", "unrelaxed_model_1_multimer_v3_pred_0.pdb"):
    fullName = os.path.join(outDir,fn)
    if not (os.path.isfile(fullName)):
      with open(fullName, "a") as fp:
        pass
    os.chmod(fullName, 0o444)
  

def BP_setupJob(args):
  jobID = args.job_id
  with open(args.jobTable) as fp:
    for line in fp.readlines():
      #print(line.strip())
      j,seq1,seq2 = [word.strip() for word in line.strip().split(",")]
      if(int(j) == jobID):
        break
  assert int(j) == jobID, f"job_id {jobID} not found in AlphaFoldJobList at {args.jobTable}"
  
  outDir = alphaFoldRunOutputDirectory(seq1, seq2, args.output_dir)
  fastaPath = os.path.join(outDir, f"{os.path.split(outDir)[-1]}.fasta")
  
  if (args.setup_job):  
    seqs, inOrder = read_fasta(args.master_fasta)
  
    # setup output directory
    try:
      os.makedirs(os.path.join(outDir, "msas")) # need one subdirectory
    except FileExistsError:
        pass
    print (outDir)

    # create fasta for input
    print (fastaPath)
    with open(fastaPath, "w") as fp:
      fp.write(f">{seq1}\n{seqs[seq1]}\n>{seq2}\n{seqs[seq2]}\n")


    # find pre-built MSAs in alignmentRepo, and copy to output dir
    # we assume first sequence will be A and second will be B
    for chainID,seqName in dict(A = seq1, B = seq2).items():
      msaDir = getMSADirectoryForSequence(seqs[seqName], args.alignmentRepo)
      localDir = os.path.join(outDir, "msas", chainID)

      if os.path.isdir(msaDir) and not os.path.isdir(localDir):
        print (f"copying from {msaDir} for {chainID}:{seqName} to {localDir}")
        shutil.copytree(msaDir, localDir)
      else:
          print (f"Not copying MSA to local {chainID}:{seqName} at {msaDir} because {'no archive MSA' if not os.path.isdir(msaDir) else 'local msa exists'}")
  
  if args.prevent_alphafold_output:
    alphaFoldLockFiles(outDir, lock = True)
        
  return fastaPath
def BP_notYetCompleted(args):
  if (args.check_if_completed == False):
    return True
  outDir = os.path.split(args.fasta_paths)[0] # I store the fasta in the otuput directory:  ./output/A123_B456/A123_B456.fasta
  if os.path.isfile(os.path.join(outDir, "result_model_5_multimer_v3_pred_0.pkl")):
    return False 
  return True


def BP_archiveMSAs(args):
  outDir = os.path.split(args.fasta_paths)[0] # I store the fasta in the otuput directory:  ./output/A123_B456/A123_B456.fasta
  # see if there are succesful msas to copy over...
  # we use the presence of file features.pkl as evidence that msas completed.
  # this is probably over-stringent in that it only appears after all MSAs are completed.
  if not os.path.isfile(os.path.join(outDir, "features.pkl")):
    print ("MSAs appear not to have finished, not copying to msa repository")
    print (os.path.join(outDir, "features.pkl"))
  else:
    dimerSeqs,namesInOrder = read_fasta(args.fasta_paths)
    print (namesInOrder)
    chainIDs = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    for i,name in enumerate(namesInOrder):
      chainID = chainIDs[i]
      msaRepoDir = getMSADirectoryForSequence(dimerSeqs[name], args.alignmentRepo)
      msaRunDir = os.path.join(outDir, "msas", chainID)
      if (os.path.isdir(msaRepoDir)):
          repoFiles = os.listdir(msaRepoDir)
      else:
          repoFiles = []
      if (os.path.isdir(msaRunDir)):
          runFiles = os.listdir(msaRunDir)
      else:
          runFiles = []
      # if any files are in run and not in repo, copy whole run directory over
      if len(set(runFiles) - set(repoFiles)) > 0:
        print (f"Copying {msaRunDir} to {msaRepoDir}; {chainID}:{name}")
        shutil.copytree(msaRunDir, msaRepoDir)
        # also write a fasta with current seq as a record...
        with open(os.path.join(msaRepoDir, f"{name}.fasta"), "w") as fp:
          fp.write(f">{name}\n{dimerSeqs[name]}\n")
      else:
        print (f"NOT Copying {msaRunDir} to {msaRepoDir}; {chainID}:{name}")

def BP_processScores(args):
  outputDir = os.path.split(args.fasta_paths)[0] # I store the fasta in the otuput directory:  ./output/A123_B456/A123_B456.fasta

  singularity_args = ['singularity',
          'exec',
          '--nv',  # Use Nvidia container library to use CUDA
          #'-B "%s"' % args.mount_data_dir,    # Mount AlphaFold databases
          '-B "%s"' % os.getcwd(),	# Mount current directory for sequence
          #'-B "%s"' % tempdir,		# Mount scratch directory
          #'--env %s' % env_vals, 
          args.singularity_image_path,
          "python3",
          "getPickleScores.py",
          outputDir
        ]
  cmd = ' '.join(singularity_args)
  print ("\nScore Summarization Command")
  print (cmd)
  print ()

  from subprocess import run
  import sys
  run(#'module load cuda/11.0 ; %s' % 
      cmd,
      stdout = sys.stdout, stderr = sys.stderr,
      shell = True,  # module command is a csh alias on Wynton
      executable = '/bin/csh',
      check = True)




###############################################################



if __name__ == '__main__':
  main()
