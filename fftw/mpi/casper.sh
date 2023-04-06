#!/bin/bash -l
### Job Name
#PBS -N test_fftw
#PBS -A NTDD0004
#PBS -l select=2:ncpus=128:mpiprocs=128
#PBS -l walltime=01:00
#PBS -q casper
#PBS -j oe

module load cray-mpich/8.1.21

mpirun -n 256  ./testx > test_fftw.test.out
