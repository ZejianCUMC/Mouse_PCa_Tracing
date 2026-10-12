#!/bin/bash
#SBATCH --job-name annotate_npp53
#SBATCH --output=/users/zw2994/OE/annotate_npp53_%j.out
#SBATCH --error=/users/zw2994/OE/annotate_npp53_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=0-03:00:00
source /groups/ms3625_gp/zw2994/mamba/etc/profile.d/conda.sh
export PATH=/groups/ms3625_gp/zw2994/mamba/envs/ST/bin:$PATH
ulimit -s unlimited
Rscript /groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/scripts/annotate_npp53.R
