#!/bin/bash

# Name job 'poisson'
#PBS -N poisson

# Allocate two nodes with 12 processors from the default resources
#PBS -lnodes=3:ppn=12:default

# Expect to run up to 1 hour
#PBS -lwalltime=01:00:00

# Memory per process
#PBS -lpmem=2000MB

# Run on the freecycle account
#PBS -A freecycle

# Run in the optimist queue by default
#PBS -q optimist

# Join stdout and stderr output to one file
#PBS -j oe

# Change directory to dir with the job script
cd ${PBS_O_WORKDIR}

# Load needed modules
module load intelcomp
module load openmpi/1.4.3-intel

# Set thread affinity
KMP_AFFINITY="granularity=fine,compact"

# Run with 8 MPI processes, each with 3 threads
OMP_NUM_THREADS=2 mpirun -npernode 6 poisson 16384 10
