#!/bin/bash
#SBATCH --job-name=build_scdbl_env
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53Rb1_Tracing/logs/%x_%j.out
#SBATCH --error=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53Rb1_Tracing/logs/%x_%j.err
source /groups/ms3625_gp/zw2994/mamba/etc/profile.d/conda.sh
mamba create -y -p /groups/ms3625_gp/zw2994/mamba/envs/scdbl -c conda-forge -c bioconda \
    bioconductor-scdblfinder bioconductor-dropletutils r-matrix
/groups/ms3625_gp/zw2994/mamba/envs/scdbl/bin/Rscript -e 'suppressPackageStartupMessages({library(scDblFinder);library(DropletUtils)}); cat("scDblFinder", as.character(packageVersion("scDblFinder")), "R", as.character(getRversion()), "OK\n")'
