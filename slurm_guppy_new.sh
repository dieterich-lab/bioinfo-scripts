#!/bin/bash

# Settings to be modified:

# gres setup:
# - gpu:tesla is currently fixed as we only have GPUs as extract resource
# - Moreover we only have cards of the tesla type
# - The number after the clon specifies the number of GPUs required,
# e.g. something between 1 and 4

#SBATCH --gres=gpu:tesla:1

# Modifiy other SLURM variables as needed

#SBATCH --partition=gpu
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=40G

# Template for SLURM GPU handling scripts
# From https://techfak.net/gpu-cluster

# please do not modify the following parts

################################################################################################
echo "==== Start of GPU information ===="

CUDA_DEVICE=$(echo "$CUDA_VISIBLE_DEVICES," | cut -d',' -f $((SLURM_LOCALID + 1)) );
T_REGEX='^[0-9]$';
if ! [[ "$CUDA_DEVICE" =~ $T_REGEX ]]; then
        echo "error no reserved gpu provided"
        exit 1;
fi

# Print debug information

echo -e "SLURM job:\t$SLURM_JOBID"
echo -e "SLURM process:\t$SLURM_PROCID"
echo -e "SLURM GPU ID:\t$SLURM_LOCALID"
echo -e "CUDA_DEVICE ID:\t$CUDA_DEVICE"
echo -e "CUDA_VISIBLE_DEVICES:\t$CUDA_VISIBLE_DEVICES"
echo "Device list:"
echo "$(nvidia-smi --query-gpu=name,gpu_uuid --format=csv -i $CUDA_VISIBLE_DEVICES | tail -n +2)"
echo "==== End of GPU information ===="
echo ""
#################################################################################################

# actual script calls start here
# $CUDA_VISIBLE_DEVICES contains the GPU IDs assigned by SLURM separated with comma
# e.g. "0,1,2"

module load guppy

# check if we have 2 arguments
if [ ! $# == 5 ]; then
  echo "Usage: $0 [Input folder (recursive)] [Output folder] [Kit] [Flow cell type] [model]"
  exit
fi

guppy_basecaller	--qscore_filtering \
                        --trim_strategy 'none'\
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
			--device auto
