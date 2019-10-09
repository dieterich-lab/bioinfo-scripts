#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @Project: University Hospital Heidelberg, Section of Bioinformatics and Systems Cardiology
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 40
#SBATCH --mem=250G
#SBATCH -J "circtools alignment"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de

#module load star

# check if we have 5 arguments
if [ ! $# == 5 ]; then
  echo "Usage: $0 [STAR index] [Read 1 file] [target dir e.g. /tmp/] [Read 1 marker, e.g. R1] [GTF file]"
  exit
fi

# $1 -> Genome index
# $2 -> Read 1
# $3 -> Read 2
# $4 -> Target directory

# remove the file extension and potential "R1" markings
# (works for double extension, e.g. .fastq.gz)
target=$4

# create the target directory, STAR will not do that for us
mkdir -pv $3/$target

# create random string
TMP_RND=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 10 | head -n 1`

OLD_PATH=`pwd`

# main mapping part

STAR	--runThreadN 40\
	--genomeDir $1\
	--genomeLoad NoSharedMemory\
	--outTmpDir /scratch/global_tmp/${TMP_RND}_${target}/\
	--readFilesIn $2 \
	--readFilesCommand zcat\
	--outFileNamePrefix $3/$target/\
	--outSAMattributes NH   HI   AS   nM   NM   MD   jM   jI   XS\
	--outSJfilterOverhangMin 15   15   15   15\
	--outFilterMultimapNmax 20\
	--chimMultimapNmax 20\
	--outFilterScoreMin 1\
	--outFilterMatchNminOverLread 0.7\
	--outFilterMismatchNmax 999\
	--outFilterMismatchNoverLmax 0.05\
	--alignIntronMin 20\
	--alignIntronMax 1000000\
	--alignMatesGapMax 1000000\
	--alignSJoverhangMin 15\
	--alignSJDBoverhangMin 10\
	--alignSoftClipAtReferenceEnds No\
	--chimSegmentMin 15\
	--chimScoreMin 15\
	--chimScoreSeparation 10\
	--chimJunctionOverhangMin 15\
	--sjdbGTFfile $5\
	--quantMode GeneCounts\
	--twopassMode Basic\
	--chimOutType Junctions

cd $3/$target

awk 'BEGIN {OFS="\t"} {split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}' Aligned.out.sam > Aligned.noS.sam

grep "^@" Aligned.out.sam > header.txt

rm -f Aligned.out.sam

rm -f -r _STARgenome
rm -f -r _STARpass1

samtools view -bS Aligned.noS.sam | samtools sort -@ 10 -m 2G -T tempo -o Aligned.noS.bam /dev/stdin
samtools reheader header.txt Aligned.noS.bam > Aligned.noS.tmp
mv Aligned.noS.tmp Aligned.noS.bam
samtools index Aligned.noS.bam

rm -f Aligned.noS.sam

cd $OLD_PATH

# remove tmp dirs

rm /scratch/global_tmp/${TMP_RND}_${target}/ -rf
