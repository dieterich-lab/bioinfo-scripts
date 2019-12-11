#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=40G
#SBATCH -J "BAM -> BW"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de

samtools index ${1}
bamCoverage -b ${1} -o ${1}.bw

#bamTobw.sh -b ${1}
