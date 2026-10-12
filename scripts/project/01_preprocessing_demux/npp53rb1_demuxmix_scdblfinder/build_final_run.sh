#!/bin/bash
#SBATCH --job-name=build_final_singlets
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=0-01:00:00
#SBATCH --output=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53Rb1_Tracing/logs/%x_%j.out
#SBATCH --error=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53Rb1_Tracing/logs/%x_%j.err
source /groups/ms3625_gp/zw2994/mamba/etc/profile.d/conda.sh
export PATH=/groups/ms3625_gp/zw2994/mamba/envs/ST/bin:$PATH
for lib in MJZ025 MJZ026 MJZ027 MJZ028; do python3 /groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53Rb1_Tracing/scripts/build_final_mtx.py $lib; done
Rscript /groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53Rb1_Tracing/scripts/build_final_seurat.R
