#!/bin/sh
#PBS -P ACD112153
#PBS -N out
#PBS -q ctest
#PBS -l select=4:ncpus=32:mpiprocs=32
#PBS -l place=scatter
#PBS -l walltime=00:01:00
#PBS -j n
module purge
module load mpi/openmpi-3.0.0/gcc485
cd $PBS_O_WORKDIR
echo $PBS_O_WORKDIR
date
mpicxx 9-10.cpp -o 9-10 -lm
mpirun ./9-10
