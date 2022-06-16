#!/bin/bash

### ichthyosaur occurrences
### 10 replicates array job
### name of job
#SBATCH --job-name=TPP-subages
#SBATCH --output=TPP-subages-out
### time to stop job
#SBATCH --time=0-48:00:00
### number of nodes/cpus *for each job*
### in this case each run of PyRate uses 1 CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=4gb
#SBATCH --cpus-per-task=1
#SBATCH --partition=Innovation

### index range for subjobs
### this is how many runs of PyRate we want
#SBATCH --array=1-24

### add Python 3.8 module
module add lang/python/anaconda/3.8.3-2020-math

### export useful things
export RUNDIR="/user/home/glbcm/ichthy_pyrate"
export PYRATESCRIPT="/user/home/glbcm/PyRate/PyRate.py"
export PYRATEINPUT="/user/home/glbcm/ichthy_pyrate/ichthyosaurs_PyRate.py"

cd $RUNDIR

### run main PyRate analysis
python3 $PYRATESCRIPT $PYRATEINPUT -n 100000000 -s 5000 -j $SLURM_ARRAY_TASK_ID -out TPP-subages -A 4 -qShift ./substages -mG -pP 1.5 0

### check the likelihood
python3 $PYRATESCRIPT -mProb './pyrate_mcmc_logs/ichthyosaur_'$SLURM_ARRAY_TASK_ID'TPP-subages_G_mcmc.log' -b 200 > $SLURM_ARRAY_TASK_ID'TPP-subages-prob'
