#!/bin/bash
### PBDB triassic ichthyosaur occurrences
### 10 replicates array job
### name of job
#SBATCH --job-name=Triassic_ichthyosaurs
#SBATCH --output=Triassic_ichthyosaurs_out
### time to stop job
#SBATCH --time=0-24:00:00
### number of nodes/cpus *for each job*
### in this case each run of PyRate uses 1 CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=4gb
#SBATCH --cpus-per-task=1
#SBATCH --partition=innovation
### index range for subjobs
### this is how many runs of PyRate we want
#SBATCH --array=1-10
### add Python 3.8 module
module add lang/python/anaconda/3.8.3-2020-math
python3 /user/home/glbcm/PyRate/PyRate.py -n 100000000 -s 1000 -j $SLURM_ARRAY_TASK_ID -edgeShift 260 200 -out triassic_ichthyosaurs /user/home/glbcm/ichthy_pyrate/triassic_species_PyRate.py -mG -A 2
