#!/bin/bash
#SBATCH --job-name annot_detailed
#SBATCH --output=/users/zw2994/OE/annot_detailed_%j.out
#SBATCH --error=/users/zw2994/OE/annot_detailed_%j.err
#SBATCH --nodes=1 --cpus-per-task=4 --mem=64G --time=0-00:30:00
source /groups/ms3625_gp/zw2994/mamba/etc/profile.d/conda.sh
/groups/ms3625_gp/zw2994/mamba/envs/ST/bin/Rscript /groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/scripts/annotate_detailed_npp53.R
