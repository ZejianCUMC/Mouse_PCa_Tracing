#!/bin/bash
#SBATCH --job-name lock_scvi
#SBATCH --output=/users/zw2994/OE/lock_scvi_%j.out
#SBATCH --error=/users/zw2994/OE/lock_scvi_%j.err
#SBATCH --nodes=1 --cpus-per-task=4 --mem=48G --time=0-00:20:00
source /groups/ms3625_gp/zw2994/mamba/etc/profile.d/conda.sh
/groups/ms3625_gp/zw2994/mamba/envs/ST/bin/Rscript /groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/scripts/set_scvi_primary.R
