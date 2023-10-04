#!/bin/bash
#SBATCH -J Name_of_job
#SBATCH -o job.o%j
#SBATCH -t 24:00:00
#SBATCH -N 10 --tasks-per-node=42 
#SBATCH	--mail-user=email@edu.com
#SBATCH --mail-type=all
#SBATCH -A Account
#SBATCH --mem=90GB

module add FFTW HDF5

mpirun -np 420 ./boutnp

    
