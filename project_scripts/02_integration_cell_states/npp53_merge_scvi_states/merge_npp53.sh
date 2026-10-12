#!/bin/bash
#SBATCH --job-name merge_npp53
#SBATCH --output=/users/zw2994/OE/merge_npp53_%j.out
#SBATCH --error=/users/zw2994/OE/merge_npp53_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=0-08:00:00
source /groups/ms3625_gp/zw2994/mamba/etc/profile.d/conda.sh
export PATH=/groups/ms3625_gp/zw2994/mamba/envs/ST/bin:$PATH
ulimit -s unlimited
Rscript /groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/scripts/merge_npp53.R
