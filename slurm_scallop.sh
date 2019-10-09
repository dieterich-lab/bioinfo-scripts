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
#SBATCH -c 40
#SBATCH --mem=500G
#SBATCH -J "scallop"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de
#SBATCH -p long

module load scallop
module load gffread
module load rnaseqtools

# check if we have 6 arguments
if [ ! $# == 4 ]; then
  echo "Usage: $0 [Input BAM] [Output GTF file] [Reference genome GTF] [Reference genome FASTA]"
  exit
fi

scallop -i ${1} -o ${2} --min_transcript_coverage 24 --verbose 1

# Step 3: Use gffcompare to evaluate the assembled transcripts using a reference annotation:
gffcompare -o ${2}.all -r ${3} ${2}

# Step 4: Union the assembled transcripts with the reference transcriptome.
# Specifically, First, use our tool gtfcuff to fetch the transcripts that are only in scallop.gtf:
gtfcuff puniq ${2}.all.${2}.tmap ${2} ${3} ${2}.unique.gtf

# The uniquely expressed transcripts (i.e., those are in scallop.gtf but not in reference.gtf) will be written to unique.gtf.
# Second, extract the cDNA sequences of the transcripts in unique.gtf from a reference genome using tool gffread:
gffread ${2}.unique.gtf -g ${4} -w ${2}.unique.fa

# where genome is the reference genome, for example ensembl reference genome.
# The cDNA sequences of the uniquely assembled transcripts (i.e., those in unique.gtf) will be written to unique.fa. 
# Finally, merge unique.fa and the reference transcriptome to obtained the extended transcriptome:
# cat ${2}.unique.fa ${4} > ${2}.union.fa

# Step 5: Run Salmon to quantify with respect to the above extended transcriptome. First, create Salmon index:
#salmon index -t ${2}.union.fa -i salmon.index -p 40

# quantify
#salmon quant -i salmon.index -1 ${5} -2 ${6} -p 4
