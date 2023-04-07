#!/bin/bash -l
### Job Name
#PBS -N test_fftw
#PBS -A UNDM0006
#PBS -l select=4:ncpus=128:mpiprocs=128
#PBS -l walltime=01:00:00
#PBS -q main@gusched01
#PBS -j oe

module --force purge
module load ncarenv/23.03
module load intel/2023.0.0
module load ncarcompilers/0.8.0
module load cray-mpich/8.1.25
module load fftw/3.3.10

mpiexec -n 512 ./testx >& test_fftw.test.log
