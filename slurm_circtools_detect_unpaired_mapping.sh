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
#SBATCH --mem=200G
#SBATCH -J "circtools alignment"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de

module load star
module load samtools

# check if we have 6 arguments
if [ ! $# == 4 ]; then
  echo "Usage: $0 [STAR index] [Read 1 file] [target dir e.g. /tmp/] [GTF file]"
  exit
fi

# $1 -> Genome index
# $2 -> Read 1
# $3 -> Target directory

# remove the file extension and potential "R1" markings
# (works for double extension, e.g. .fastq.gz)
target=`expr ${2//} : '\(.*\)\..*\.'`

# create the target directory, STAR will not do that for us
mkdir -pv $3/$target

STAR --runThreadN 20 \
       --genomeDir $1 \
       --readFilesIn $2 \
       --readFilesCommand zcat \
       --outFileNamePrefix $3/$target/ \
       --outSJfilterOverhangMin 5 5 5 5 \
       --alignSJoverhangMin 5 \
       --alignSJDBoverhangMin 5 \
       --outFilterMultimapNmax 20 \
       --outFilterScoreMin 1 \
       --outFilterMatchNmin 1 \
       --outFilterMismatchNmax 2 \
       --chimSegmentMin 5 \
       --chimScoreMin 5 \
       --chimScoreSeparation 10 \
       --chimJunctionOverhangMin 5 \
       --chimOutType Junctions SeparateSAMold\
       --sjdbGTFfile $4 

cd $3/$target

awk 'BEGIN {OFS="\t"} {split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}' Aligned.out.sam > Aligned.noS.sam

awk 'BEGIN {OFS="\t"} {split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}' Chimeric.out.sam > Chimeric.noS.sam

grep "^@" Aligned.out.sam > header.txt

#rm -f Aligned.out.sam
#rm -f Chimeric.out.sam

#rm -f -r _STARgenome
#rm -f -r _STARpass1

samtools view -bS Aligned.noS.sam | samtools sort -@ 10 -m 2G -T tempo -o Aligned.noS.bam /dev/stdin
samtools reheader header.txt Aligned.noS.bam > Aligned.noS.tmp
mv Aligned.noS.tmp Aligned.noS.bam
samtools index Aligned.noS.bam

samtools view -bS Chimeric.noS.sam | samtools sort -@ 10 -m 2G -T tempo -o Chimeric.noS.bam /dev/stdin
samtools reheader header.txt Chimeric.noS.bam > Chimeric.noS.tmp
mv Chimeric.noS.tmp Chimeric.noS.bam
samtools index Chimeric.noS.bam

#rm -f Aligned.noS.sam
#rm -f Chimeric.noS.sam


