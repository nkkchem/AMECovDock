#!/bin/bash
#SBATCH -A mtcovid
#SBATCH -N 1
#SBATCH -J Refinement
#SBATCH -t 48:00:00
#SBATCH -n 8

module load python/anaconda3.2019.3
source /share/apps/python/anaconda3.2019.3/etc/profile.d/conda.sh
conda activate my-rdkit-env
module load gcc/6.1.0
module load openbabel/2.4.1

# change autodock to have the correct files
./AMECovDock_ADFR.sh
