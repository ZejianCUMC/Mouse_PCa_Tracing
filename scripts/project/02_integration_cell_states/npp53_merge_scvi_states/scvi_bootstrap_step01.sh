#!/bin/bash
#SBATCH --job-name boot_step01
#SBATCH --output=/users/zw2994/OE/boot_step01_%j.out
#SBATCH --error=/users/zw2994/OE/boot_step01_%j.err
#SBATCH --nodes=1 --cpus-per-task=16 --mem=128G --time=0-02:30:00
source /groups/ms3625_gp/zw2994/mamba/etc/profile.d/conda.sh
export PATH=/groups/ms3625_gp/zw2994/mamba/envs/ST/bin:$PATH
ulimit -s unlimited
Rscript /groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/scripts/scvi_bootstrap_step01.R
