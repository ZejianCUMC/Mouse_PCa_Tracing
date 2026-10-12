#!/bin/bash
#SBATCH --job-name scvi_dense
#SBATCH --output=/users/zw2994/OE/scvi_dense_%j.out
#SBATCH --error=/users/zw2994/OE/scvi_dense_%j.err
#SBATCH --nodes=1 --cpus-per-task=6 --mem=64G --time=0-00:30:00
source /groups/ms3625_gp/zw2994/mamba/etc/profile.d/conda.sh
/groups/ms3625_gp/zw2994/mamba/envs/ST/bin/Rscript /groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/scripts/scvi_dense_umap.R
