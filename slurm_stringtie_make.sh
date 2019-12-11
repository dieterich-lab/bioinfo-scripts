#!/bin/bash

# Copyright (C) 2017 Tobias Jakobi
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either self.version 3 of the License, or
# (at your option) any later self.version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @Project: University Hospital Heidelberg, Section of Bioinformatics and Systems Cardiology

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 10
#SBATCH --mem=10G
#SBATCH -J "Stringtie"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de

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
mkdir $3 -p

stringtie $1 --rf -v -f 0.2 -p 10 -G $2 -e -B -o $3/ballgown.gtf -C $3/reference_transcripts_full_coverage.gtf
