#!/bin/bash
#SBATCH --job-name=SlantwiseTest1
#SBATCH --output=slurm.out
#SBATCH --error=slurm.err
#SBATCH --time=0:59:59
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --exclusive


module unload intel intelmpi
module unload netcdf
module load netcdf/4.4.1+intel-16.0
module load netcdf-fortran/4.4.4+intel-16.0
module load intelmpi/5.1+intel-16.0

mpirun -n 1 ./mitgcmuv





