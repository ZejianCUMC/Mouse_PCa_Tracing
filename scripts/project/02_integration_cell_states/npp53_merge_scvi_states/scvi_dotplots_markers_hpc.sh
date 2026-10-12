#!/bin/bash
#SBATCH --job-name scvi_de
#SBATCH --output=/users/zw2994/OE/scvi_de_%j.out
#SBATCH --error=/users/zw2994/OE/scvi_de_%j.err
#SBATCH --nodes=1 --cpus-per-task=16 --mem=128G --time=0-01:30:00
source /groups/ms3625_gp/zw2994/mamba/etc/profile.d/conda.sh
export PATH=/groups/ms3625_gp/zw2994/mamba/envs/ST/bin:$PATH
mkdir -p /groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/npp53_plots
ulimit -s unlimited
Rscript /groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/scripts/scvi_dotplots_markers_hpc.R
