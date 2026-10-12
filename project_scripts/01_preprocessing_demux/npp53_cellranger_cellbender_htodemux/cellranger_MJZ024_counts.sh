#!/bin/bash
#SBATCH --mem=96G
#SBATCH --time=55:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --job-name=MJZ024_align
#SBATCH --error=/users/zw2994/OE/cellranger_MJZ024.err
#SBATCH --output=/users/zw2994/OE/cellranger_MJZ024.out

cd /groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-counts

/groups/ms3625_gp/zw2994/tools/cellranger-10.0.0/cellranger count --id=MJZ024 \
                   --libraries=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/cellranger_MJZ19_22_prep/MJZ024_library.csv \
                   --transcriptome=/groups/ms3625_gp/zw2994/0-database/GRCm39_tracing \
                   --feature-ref=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/cellranger_MJZ19_22_prep/MJZ024_feature_ref.csv \
                   --localcores=12 \
                   --localmem=94 \
                   --create-bam true
