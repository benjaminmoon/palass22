#!/bin/bash
### combine TPP subages logs
#SBATCH --job-name=combine-TPP-subages
#SBATCH --output=combine-TPP-subages-out
#SBATCH --time=0-02:00:00
# 1 node, 1 cpu, 8gb mem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8gb
#SBATCH --cpus-per-task=1
#SBATCH --partition=Innovation

# load python and R modules
module add lang/python/anaconda/3.8.3-2020-math
module add lang/r/3.6.1

export RAPPLICATION="/sw/lang/R-3.6.1-intel/bin/Rscript"
export ESSSCRIPT="/user/home/glbcm/ichthy_pyrate/bpbc_array_ESS.R"
export RUNDIR="/user/home/glbcm/ichthy_pyrate"
export PYRATESCRIPT="/user/home/glbcm/PyRate/PyRate.py"

cd $RUNDIR

# combine PyRate logs
# python3 $PYRATESCRIPT -combLogRJ ./pyrate_mcmc_logs/ -b 100 -tag TPP-subages

# create plot
python3 $PYRATESCRIPT -plotRJ ./pyrate_mcmc_logs/ -grid_plot 0.1 -tag combined_10TPP-subages

# calculate ESS and plot traces
Rscript $ESSSCRIPT
