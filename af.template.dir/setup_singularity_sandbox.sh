#!/bin/bash


# create a sandbox image that we can modify and run from
singularity build --sandbox ../alphafoldSandbox /wynton/home/ferrin/goddard/alphafold_singularity/alphafold231.sif


# copy the modified model config files
cp ./alt_model_configs/config* ../alphafoldSandbox/app/alphafold/alphafold/model/

# copy the modified startup script
# this startup script looks for a --startModel= parameter
cp ./alt_model_configs/run_alphafold.sh   ../alphafoldSandbox/app/run_alphafold.sh 

