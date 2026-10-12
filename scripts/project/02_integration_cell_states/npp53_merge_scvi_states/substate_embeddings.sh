#!/bin/bash
#SBATCH --job-name substate_emb
#SBATCH --output=/users/zw2994/OE/substate_emb_%j.out
#SBATCH --error=/users/zw2994/OE/substate_emb_%j.err
#SBATCH --nodes=1 --cpus-per-task=6 --mem=64G --time=0-00:40:00
source /groups/ms3625_gp/zw2994/mamba/etc/profile.d/conda.sh
/groups/ms3625_gp/zw2994/mamba/envs/ST/bin/Rscript /groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/scripts/substate_embeddings.R
