#!/bin/bash
#SBATCH --job-name relabel_scvi
#SBATCH --output=/users/zw2994/OE/relabel_scvi_%j.out
#SBATCH --error=/users/zw2994/OE/relabel_scvi_%j.err
#SBATCH --nodes=1 --cpus-per-task=6 --mem=80G --time=0-01:00:00
source /groups/ms3625_gp/zw2994/mamba/etc/profile.d/conda.sh
export PATH=/groups/ms3625_gp/zw2994/mamba/envs/ST/bin:$PATH
ulimit -s unlimited
Rscript /groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/scripts/relabel_scvi_native.R
