#!/bin/bash
#SBATCH --job-name scvi_npp53
#SBATCH --output=/users/zw2994/OE/scvi_npp53_%j.out
#SBATCH --error=/users/zw2994/OE/scvi_npp53_%j.err
#SBATCH --partition=gpu
#SBATCH --ntasks=1
#SBATCH --gpus-per-task=1
#SBATCH --mem-per-gpu=64G
#SBATCH --time=0-12:00:00
source /groups/ms3625_gp/zw2994/mamba/etc/profile.d/conda.sh
export PATH=/groups/ms3625_gp/zw2994/mamba/envs/scvi/bin:$PATH
module load cuda
python /groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/scripts/scvi_npp53.py
