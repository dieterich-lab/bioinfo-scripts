#!/bin/bash
#SBATCH --job-name=rpl
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de
#SBATCH --mail-type=FAIL,END
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G



parallel -j1 --xapply rpl {1} {2} ${1} :::: gene_ids.txt :::: mappings.txt 
