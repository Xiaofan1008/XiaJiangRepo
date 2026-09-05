#!/bin/bash
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=sabrina.meikle@monash.edu
#SBATCH --job-name=Loop_sigchn
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=30000
module load matlab
matlab -nodisplay -nojvm -nosplash < loopsigchn.m