#!/bin/bash 

#SBATCH --job-name=leggcon0


#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=24
#SBATCH --mem=2gb
#SBATCH --cpus-per-task=1
#SBATCH --time=13-24:00:00
 
cd /user/work/bv20692/bayesian100522/   

module load apps/mrbayes/3.3.7a 

mpiexec -n 16 mb Nature2013_con0.nex 
