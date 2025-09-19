#!/bin/bash
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sabrina.meikle@monash.edu
#SBATCH --job-name=DeleteNonOriginalDatafiles
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8000
module load matlab
matlab -nodisplay -nojvm -nosplash < Delete_non_originalFilesMASSIVE.m
