#!/bin/bash
#SBATCH --job-name rebuild_npp53
#SBATCH --output=/users/zw2994/OE/rebuild_npp53_%j.out
#SBATCH --error=/users/zw2994/OE/rebuild_npp53_%j.err
#SBATCH --nodes=1 --cpus-per-task=8 --mem=96G --time=0-02:00:00
source /groups/ms3625_gp/zw2994/mamba/etc/profile.d/conda.sh
export PATH=/groups/ms3625_gp/zw2994/mamba/envs/ST/bin:$PATH
ulimit -s unlimited
Rscript /groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/scripts/rebuild_npp53.R
