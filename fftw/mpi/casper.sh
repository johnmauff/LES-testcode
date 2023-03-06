#!/bin/bash -l
### Job Name
#PBS -N test_fftw
#PBS -A NTDD0004
#PBS -l select=2:ncpus=32:mpiprocs=32
#PBS -l walltime=01:00
#PBS -q casper
#PBS -j oe

mpirun -n 64 ./a.out > test_fftw.test.out
