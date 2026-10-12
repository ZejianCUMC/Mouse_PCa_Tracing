#!/bin/bash
#SBATCH --job-name rerender_sub
#SBATCH --output=/users/zw2994/OE/rerender_sub_%j.out
#SBATCH --error=/users/zw2994/OE/rerender_sub_%j.err
#SBATCH --nodes=1 --cpus-per-task=4 --mem=48G --time=0-00:20:00
source /groups/ms3625_gp/zw2994/mamba/etc/profile.d/conda.sh
/groups/ms3625_gp/zw2994/mamba/envs/ST/bin/Rscript /groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/scripts/rerender_substate_umap.R
