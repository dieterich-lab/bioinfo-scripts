#!/bin/bash
#SBATCH --job-name=blastn
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de
#SBATCH --mail-type=FAIL,END
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=50G
 
blastn -query ${1} -db /scratch/tjakobi/ncbi_star/mouse_blast/Mus_musculus.GRCm38.cdna.all.fa -outfmt 6 -out ${1}.csv -dust no -task blastn -num_threads 10  -max_target_seqs 1
