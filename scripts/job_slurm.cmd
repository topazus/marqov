#!/bin/bash

# --- MARQOV job script for SLURM ---

# set job name
#SBATCH -J MARQOV

# set queue
#SBATCH -p standard

# number of nodes
#SBATCH -N 8

# number of tasks
#SBATCH -n 8

# cpus per task
#SBATCH -c 32

# tasks per node
# set to one since MARQOV is a hybrid code
#SBATCH --ntasks-per-node=1

# set max wallclock time
#SBATCH -t 2-00:00:00

# logfiles
#SBATCH -o ./output.%a.out
#SBATCH -e ./output.%a.err

# run
srun ./build/src/main
