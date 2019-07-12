#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @Project: University Hospital Heidelberg, Section of Bioinformatics and Systems Cardiology
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 40 
#SBATCH --mem=200G
#SBATCH -J "guppy basecaller"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de
#SBATCH -p gpu

module load guppy

# check if we have 2 arguments
if [ ! $# == 4 ]; then
  echo "Usage: $0 [Input folder (recursive)] [Output folder] [Kit] [Flow cell type]"
  exit
fi

guppy_basecaller	--qscore_filtering \
			--verbose_logs \
			--compress_fastq \
			--fast5_out \
			-r \
			-i ${1} \
			-s ${2} \
			--flowcell ${4} \
			--kit ${3} \
			--num_callers 4 \
			--gpu_runners_per_device 3 \
			--cpu_threads_per_caller 1 \
			--device "cuda:all"
