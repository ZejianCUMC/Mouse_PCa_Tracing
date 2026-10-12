#!/bin/bash
#SBATCH --job-name phaseA_edited
#SBATCH --output=/users/zw2994/OE/phaseA_edited_%j.out
#SBATCH --error=/users/zw2994/OE/phaseA_edited_%j.err
#SBATCH --nodes=1 --cpus-per-task=8 --mem=96G --time=0-02:00:00
source /groups/ms3625_gp/zw2994/mamba/etc/profile.d/conda.sh
export PATH=/groups/ms3625_gp/zw2994/mamba/envs/ST/bin:$PATH
ulimit -s unlimited
Rscript /groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/scripts/phaseA_edited_substates.R
