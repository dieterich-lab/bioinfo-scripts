#!/bin/bash

# Copyright (C) 2020 Tobias Jakobi
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
#SBATCH -c 40
#SBATCH --mem=50G
#SBATCH -J "build training set" 
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de

# check if we have 8 arguments
if [ ! $# == 8 ]; then
  echo "Usage: $0 [Target folder] [Read FASTA] [READ BAM] [# unmod reads] [# mod reads] [mod read list] [unmod read list] [FAST5 dir]"
  exit
fi

target_dir=${1}__${4}_unmod_${5}_mod
fasta=$2
bam=$3
num_unmod=$4
num_mod=$5
mod_reads=$6
unmod_reads=$7
fast5=${8}

module load seqtk

# create the target directory
mkdir $target_dir -pv

cd $target_dir

ln -s $fasta .
ln -s $bam .

shuf -n $num_mod $mod_reads > ${num_mod}_random_curlcake_reads.list
shuf -n $num_unmod $unmod_reads > ${num_unmod}_random_reads.list

seqtk subseq $fasta  ${num_mod}_random_curlcake_reads.list >  ${num_mod}_random_curlcake_reads.fasta
seqtk subseq $fasta  ${num_unmod}_random_reads.list >  ${num_unmod}_random_reads.fasta

sed 's/A/Y/g' ${num_mod}_random_curlcake_reads.fasta > ${num_mod}_random_curlcake_reads_Y.fasta

cat ${num_mod}_random_curlcake_reads.list ${num_unmod}_random_reads.list > merged_reads.list
cat ${num_mod}_random_curlcake_reads_Y.fasta ${num_unmod}_random_reads.fasta > merged_reads.fasta

fast5_subset -i $fast5 -l merged_reads.list -t 40 -r -s training_set_fast5

