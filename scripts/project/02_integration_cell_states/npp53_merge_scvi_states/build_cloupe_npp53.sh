#!/bin/bash
#SBATCH --job-name cloupe_npp53
#SBATCH --output=/users/zw2994/OE/cloupe_npp53_%j.out
#SBATCH --error=/users/zw2994/OE/cloupe_npp53_%j.err
#SBATCH --nodes=1 --cpus-per-task=6 --mem=96G --time=0-01:30:00
source /groups/ms3625_gp/zw2994/mamba/etc/profile.d/conda.sh
ulimit -s unlimited
SDIR=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/scripts
/groups/ms3625_gp/zw2994/mamba/envs/scvi/bin/python $SDIR/export_scvi_umap.py
/groups/ms3625_gp/zw2994/mamba/envs/ST/bin/Rscript   $SDIR/build_cloupe_npp53.R
