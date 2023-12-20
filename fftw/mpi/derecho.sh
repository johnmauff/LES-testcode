#!/bin/bash -l
### Job Name
#PBS -N test_fftw
#PBS -A NTDD0004
#PBS -l select=1:ncpus=8:mpiprocs=8
#PBS -l walltime=01:00
#PBS -q main
#PBS -j oe

module load cray-mpich/8.1.25

mpirun -n 8  ./testx > test_fftw.test.out
