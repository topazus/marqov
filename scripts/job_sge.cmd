#!/bin/bash
#$ -S /bin/bash

# --- MARQOV SGE job script for the ITPA cluster ---
# move to build directory and use qsub to submit

# start job in current working directory
#$ -cwd  	

# Pretend to require so much memory that the SGE is
# forced to spread the MPI process across nodes
#$ -l h_vmem=48G	# reduce for mpi16 queue

# Choose queue and number of MPI processes
#$-pe mpi32* 8

# Job name
#$-N MARQOV

#$ -e error.txt
#$ -o output.txt

echo "Starting MARQOV on $NSLOTS nodes!"
mpirun ./build/src/main
