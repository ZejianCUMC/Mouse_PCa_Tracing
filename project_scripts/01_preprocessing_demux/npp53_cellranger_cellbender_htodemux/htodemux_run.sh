#!/bin/bash
#SBATCH --job-name htodemux
#SBATCH --output=/users/zw2994/OE/htodemux_%x_%j.out
#SBATCH --error=/users/zw2994/OE/htodemux_%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=0-01:30:00
source /groups/ms3625_gp/zw2994/mamba/etc/profile.d/conda.sh
export PATH=/groups/ms3625_gp/zw2994/mamba/envs/ST/bin:$PATH
Rscript /groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/scripts/htodemux_$1.R
