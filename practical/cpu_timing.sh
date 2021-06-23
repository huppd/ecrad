#!/bin/bash -l
#SBATCH --job-name="ecrad_gpu"
#SBATCH --account="s83"
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=debug
#SBATCH --constraint=gpu
#SBATCH --hint=nomultithread
#SBATCH --exclusive 

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export CRAY_CUDA_MPS=1

srun ../bin/ecrad config_cpu.nam era5slice_2500.nc tt.nc
# srun ../bin/ecrad config_cpu.nam era5slice_6000.nc tt.nc
# srun ../bin/ecrad config_cpu.nam era5slice_10000.nc tt.nc
