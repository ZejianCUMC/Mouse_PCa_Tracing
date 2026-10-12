#!/bin/bash
#SBATCH --job-name cellbender_MJZ021
#SBATCH --output=/users/zw2994/OE/cellbender_MJZ021.out
#SBATCH --error=/users/zw2994/OE/cellbender_MJZ021.err
#SBATCH --partition=gpu
#SBATCH --ntasks=1
#SBATCH --gpus-per-task=1
#SBATCH --mem-per-gpu=48G
#SBATCH --time=0-12:10:00

source /groups/ms3625_gp/zw2994/mamba/etc/profile.d/conda.sh
export PATH=/groups/ms3625_gp/zw2994/mamba/envs/cellbender/bin:$PATH
module load cuda

cd /groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-counts/MJZ021/outs
cellbender remove-background \
      --cuda \
      --input raw_feature_bc_matrix.h5 \
      --output cellbender_output_file.h5
