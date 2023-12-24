#!/bin/sh
#PBS -P ACD112153
#PBS -N out
#PBS -q ctest
#PBS -l select=1:ncpus=8:mpiprocs=2
#PBS -l place=scatter
#PBS -l walltime=00:01:00
#PBS -j n
module purge
module load mpi/openmpi-3.0.0/gcc485
cd $PBS_O_WORKDIR
echo $PBS_O_WORKDIR
date
mpirun ./6-8
