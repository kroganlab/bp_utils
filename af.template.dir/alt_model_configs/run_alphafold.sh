#!/bin/bash
ldconfig


# find, handle, and remove the startModel argument
declare -a ARGS
for var in "$@"; do
    # look for model argument
    if [[ "$var" = '--startModel'* ]]; then
	startModel=${var:13:1}
	ln -sf /app/alphafold/alphafold/model/config$startModel.py /app/alphafold/alphafold/model/config.py
        ls -l  /app/alphafold/alphafold/model/config*
	continue
    fi
    ARGS[${#ARGS[@]}]="$var"
done


echo
echo "${ARGS[@]}"
echo

python /app/alphafold/run_alphafold.py "${ARGS[@]}" 
