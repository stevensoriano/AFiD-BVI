#!/bin/bash
#SBATCH -J ForceVsTime
#SBATCH -o python_job.o%j
#SBATCH -t 12:00:00 
#SBATCH --mail-user=sasoriano2@uh.edu
#SBATCH --mail-type=all
#SBATCH -A ostilla
#SBATCH --mem=100GB

module load python/3.8
source /home/sasoria2/miniconda3/bin/activate

python ForceOverTime.py
