#!/bin/bash -l

##$ -t 1-2900

##$ -l h_rt=60:00:00

##$ -pe omp 28

#$ -m ea

#$ -N sim_256

#$ -j y

#$ -o sim_256.qlog

module load miniconda
conda activate brienv
cd /projectnb/crc-nak/brpp/Aussel_2023/model_files
git checkout same-object
#python FEF_and_LIP_argument.py $SGE_TASK_ID $1 /projectnb/crc-nak/brpp/Aussel_2023/simulation_results
python FEF_and_LIP_no_cache.py 256 45 /projectnb/crc-nak/brpp/Aussel_2023/simulation_results