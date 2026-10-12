#!/bin/bash
#SBATCH --job-name jac_04_08
#SBATCH --output=/users/zw2994/OE/jac_04_08_%j.out
#SBATCH --error=/users/zw2994/OE/jac_04_08_%j.err
#SBATCH --nodes=1 --cpus-per-task=16 --mem=128G --time=0-01:00:00
source /groups/ms3625_gp/zw2994/mamba/etc/profile.d/conda.sh
export PATH=/groups/ms3625_gp/zw2994/mamba/envs/ST/bin:$PATH
ulimit -s unlimited
Rscript /groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/scripts/percluster_jaccard_res04_08.R
