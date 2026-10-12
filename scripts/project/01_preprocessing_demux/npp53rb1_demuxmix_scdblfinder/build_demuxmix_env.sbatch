#!/bin/bash
#SBATCH --job-name=build_demuxmix_env
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53Rb1_Tracing/logs/%x_%j.out
#SBATCH --error=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53Rb1_Tracing/logs/%x_%j.err
source /groups/ms3625_gp/zw2994/mamba/etc/profile.d/conda.sh
# let bioconda resolve a consistent R (do NOT pin r-base -> avoids the pinned-old-R breakage)
mamba create -y -p /groups/ms3625_gp/zw2994/mamba/envs/demuxmix -c conda-forge -c bioconda \
    bioconductor-demuxmix r-matrix r-ggplot2 r-gridextra
echo "=== verify ==="
/groups/ms3625_gp/zw2994/mamba/envs/demuxmix/bin/Rscript -e 'suppressPackageStartupMessages(library(demuxmix)); cat("demuxmix", as.character(packageVersion("demuxmix")), "R", as.character(getRversion()), "OK\n")'
