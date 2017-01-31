#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Date:   Friday, May 6, 2016 4:10 PM
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @Project: University Hospital Heidelberg, Section of Bioinformatics and Systems Cardiology
# @Last modified by:   tjakobi
# @Last modified time: Friday, May 6, 2016 4:22 PM
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --mem=100G
#SBATCH -J "Stringtie"

module load stringtie

# check if we have 3 arguments
if [ ! $# == 3 ]; then
  echo "Usage: $0 [BAM file] [GTF file] [target dir e.g. /tmp/]"
  exit
fi

# $1 -> BAM file
# $2 -> GTF file
# $3 -> output directory

# create the target directory
mkdir $3/ballgown/

stringtie $1 -v -p 20 -G $2 -e -B -o $3/ballgown/ballgown.gtf
