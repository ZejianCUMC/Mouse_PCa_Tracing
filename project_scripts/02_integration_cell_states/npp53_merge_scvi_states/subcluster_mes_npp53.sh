#!/bin/bash
#SBATCH --job-name subcluster_mes
#SBATCH --output=/users/zw2994/OE/subcluster_mes_%j.out
#SBATCH --error=/users/zw2994/OE/subcluster_mes_%j.err
#SBATCH --nodes=1 --cpus-per-task=6 --mem=80G --time=0-01:00:00
source /groups/ms3625_gp/zw2994/mamba/etc/profile.d/conda.sh
SDIR=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/scripts
/groups/ms3625_gp/zw2994/mamba/envs/scvi/bin/python $SDIR/export_scvi_latent.py
/groups/ms3625_gp/zw2994/mamba/envs/ST/bin/Rscript   $SDIR/subcluster_mes_npp53.R
