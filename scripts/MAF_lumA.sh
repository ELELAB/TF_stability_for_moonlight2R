#!/bin/bash

# Activate Python Virtual Environment
#. /usr/local/envs/py37/bin/activate
module load python/3.10
# Define the command you want to run
python3 tf_maf.py ../raw_data/trrust_rawdata.human.tsv ../raw_data/mutations.csv


