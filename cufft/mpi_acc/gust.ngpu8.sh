#!/bin/bash
#PBS -N cm1
#PBS -l select=2:ncpus=4:mpiprocs=4:ngpus=4:mem=100gb
#PBS -l walltime=00:10:00
#PBS -j oe
#PBS -q main@gusched01
#PBS -A UNDM0006

#module purge
#module load nvhpc cuda cray-mpich

export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED=1
export PALS_NRANKS=8
export PALS_PPN=4
export PALS_DEPTH=32
export PALS_CPU_BIND=depth

mpiexec -n 8 -ppn 4 get_local_rank ./testx >& 8a100.test_cufft.32.A.log
# nsys profile -o cufft_test --trace openacc,cuda,mpi mpiexec -n 8 ./a.out


