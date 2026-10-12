#!/bin/bash
#SBATCH --job-name plots_npp53
#SBATCH --output=/users/zw2994/OE/plots_npp53_%j.out
#SBATCH --error=/users/zw2994/OE/plots_npp53_%j.err
#SBATCH --nodes=1 --cpus-per-task=4 --mem=64G --time=0-01:30:00
source /groups/ms3625_gp/zw2994/mamba/etc/profile.d/conda.sh
export PATH=/groups/ms3625_gp/zw2994/mamba/envs/ST/bin:$PATH
Rscript /groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/scripts/plots_npp53.R
