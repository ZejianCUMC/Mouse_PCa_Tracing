#!/bin/bash
#SBATCH --job-name DF
#SBATCH --output=/users/zw2994/OE/DF_%x_%j.out
#SBATCH --error=/users/zw2994/OE/DF_%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=0-06:00:00
source /groups/ms3625_gp/zw2994/mamba/etc/profile.d/conda.sh
export PATH=/groups/ms3625_gp/zw2994/mamba/envs/ST/bin:$PATH
ulimit -s unlimited
Rscript --max-ppsize=500000 /groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/scripts/doubletfinder_per_sample.R $1
