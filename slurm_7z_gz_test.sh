#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 10
#SBATCH --mem=20G
#SBATCH -J "7zip verify"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de

module load p7zip

7z t -mmt20 $1
