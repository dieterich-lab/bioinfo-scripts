#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @Project: University Hospital Heidelberg, Section of Bioinformatics and Systems Cardiology
# @Last modified by:   tjakobi
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=64G
#SBATCH -J "macs2"
#SBATCH -c 2
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de

macs2 $@
