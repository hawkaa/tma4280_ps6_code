#!/bin/bash

# Name job 'poisson'
#PBS -N poisson_speedup_96

# Allocate two nodes with 12 processors from the default resources
#PBS -lnodes=8:ppn=12:default

# Expect to run up to 1 hour
#PBS -lwalltime=00:30:00

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

OMP_NUM_THREADS=1 mpirun -npernode 12 poisson_speedup 10 3 14
