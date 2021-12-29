#!/bin/bash
#$ -S /bin/bash

# --- MARQOV job script for the Sun Grid Engine (SGE) ---

# start job in current working directory
#$ -cwd  	

# Pretend to require so much memory that the SGE is
# forced to spread the MPI process across nodes
#$ -l h_vmem=48G	# adjust this value depending on your hardware

# Choose queue/parallel environment(PE) and number of MPI processes
#$-pe mpi32* 8

# Job name
#$-N MARQOV

#$ -e error.txt
#$ -o output.txt

echo "Starting MARQOV on $NSLOTS nodes!"
mpirun ./build/src/main
