#!/bin/bash
#SBATCH --job-name=           # Job name
#SBATCH --mail-user=dcox67@gatech.edu # E-mail address for notifications
#SBATCH --mail-type=END,FAIL         # Mail preferences
#SBATCH --nodes=1
#SBATCH --tasks=24
#SBATCH --mem-per-cpu=1gb
#SBATCH --time=00:30:00
#SBATCH --output=output_log/testing.out

## MODIFY FOLLOWING SECTION ##

# define sbatch nodes, ntasks-per-node, mem-per-cpu, time, output
module load anaconda3
conda activate openmc-env
pip install -e .
export OMP_NUM_THREADS=$(nproc)
echo $OMP_NUM_THREADS
python arc-standard.py Li4SiO4 0 beryllium 0 lead 0 2 1 #need to chagne for the new input parameters
date

