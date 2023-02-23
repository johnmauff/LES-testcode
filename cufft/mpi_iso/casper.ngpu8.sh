#!/bin/bash -l
#PBS -N test_cufft
#PBS -l select=2:ncpus=4:mpiprocs=4:ngpus=4
#PBS -l gpu_type=v100
#PBS -l walltime=02:00:00
#PBS -j oe
#PBS -A NTDD0004
#PBS -q casper

module purge ; module load ncarenv/1.3 nvhpc/22.2 openmpi/4.1.1 netcdf/4.8.1 cuda/11.6

export LD_LIBRARY_PATH=${NCAR_ROOT_CUDA}/lib64:${LD_LIBRARY_PATH}
export UCX_MEMTYPE_CACHE=n
export UCX_TLS=rc,sm,cuda_copy,cuda_ipc
export OMPI_MCA_pml=ucx
export OMPI_MCA_btl=self,vader,tcp,smcuda
export UCX_RNDV_SCHEME=put_zcopy
export UCX_RNDV_THRESH=2


mpirun -n 8 ./a.out >& test_cufft.out

