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
#SBATCH -c 10
#SBATCH --mem=400G
#SBATCH -J "Nanocompore"
#SBATCH -p general
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de

module load nanocompore

# check if we have 4 arguments
if [ ! $# == 4 ]; then
  echo "Usage: $0 [Guppy folder] [Transcript file fasta] [Basecalled FASTQ file] [target dir e.g. /tmp/]"
  exit
fi

# create the target directory
 mkdir $4 -p

# Steps from https://nanocompore.rna.rocks/data_preparation/

# align to reference
minimap2 -ax map-ont -L ${2} ${3} | samtools view -bh -F 2324 -q 10 | samtools sort -O bam > $4/aligned.bam

# build index
samtools index $4/aligned.bam

# index first with nanopolish index
nanopolish index -s ${1}/sequencing_summary.txt -d ${1}/workspace ${3}

# realign raw signal to the expected reference sequence
nanopolish eventalign --threads 10 --reads ${3} --bam $4/aligned.bam --genome ${2} --samples --print-read-names --scale-events --samples | gzip -c --best > ${4}/eventalign_reads.tsv.gz

# data has to be collapsed per kmer and indexed by NanopolishComp Eventalign_collapse
zcat ${4}/eventalign_reads.tsv.gz | NanopolishComp Eventalign_collapse --threads 10 -o ${4}/
