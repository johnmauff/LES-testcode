#!/bin/bash -l
### Job Name
#PBS -N test_fftw
#PBS -A UNDM0006
#PBS -l select=1:ncpus=128:mpiprocs=128
#PBS -l walltime=01:00:00
#PBS -q main@gusched01
#PBS -j oe

module purge
module load oneapi/2023.0.0
module load cray-mpich/8.1.21
module load fftw

which mpiexec
mpiexec -n 128 ./testx > test_fftw.test.out
