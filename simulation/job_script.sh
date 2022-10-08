#!/usr/bin/bash

#SBATCH -p lowpriority
#SBATCH -D /users/ostmeyer/volatile2/random-field-ising/simulation
#SBATCH -J NAME
#SBATCH --output=/users/ostmeyer/volatile2/random-field-ising/simulation/NAME.out --error=/users/ostmeyer/volatile2/random-field-ising/simulation/NAME.err
#SBATCH -t 1-0
#SBATCH -n 1 -c 40 --exclusive # The number of cores we need...

module purge
module load apps/R/3.6.3

export OMP_NUM_THREADS=40
export OMP_PROC_BIND=close
export OMP_PLACES=cores

echo "print SLURM_JOB_NODELIST = $SLURM_JOB_NODELIST"
pwd

time Rscript $1 "new"

exit 0
