#!/bin/bash
# =============================================================================
# create_envs.sh — build the conda/mamba environments for this tutorial on
# your own Linux machine (or your own HPC space).
#
#   bash envs/create_envs.sh            # all four envs
#   bash envs/create_envs.sh scvi       # just one: cellbender | seurat | scvi | cassiopeia
#
# Each envs/tracing_<name>.json is a normal conda environment spec written as
# JSON (JSON is valid YAML). mamba wants a .yml extension, so we copy it first.
# Versions match the lab's c2b2 envs used to validate the scripts:
#   tracing_cellbender  = HPC env cellbender     (CellBender 0.3.0, python 3.7)   step 02
#   tracing_seurat      = HPC env ST             (Seurat 5.1, DoubletFinder)      step 03
#   tracing_scvi        = HPC env scvi           (scanpy 1.11, scvi-tools 1.2)    step 04
#   tracing_cassiopeia  = HPC env cassiopeia_env (Cassiopeia 2.1 @ 48ea342)       steps 05, 06
# Cell Ranger (step 01) is not a conda package — see the note at the end.
# Tips: needs mamba (Miniforge). GPU is optional everywhere; CellBender on CPU
# works but is slow, so run step 02 on a GPU node for real data.
# =============================================================================
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
MAMBA=$(command -v mamba || command -v conda)
TMP=$(mktemp -d)
WANT=${1:-"cellbender seurat scvi cassiopeia"}

for n in ${WANT}; do
    spec=${HERE}/tracing_${n}.json
    [ -f "${spec}" ] || { echo "no spec ${spec}"; exit 1; }
    cp "${spec}" "${TMP}/tracing_${n}.yml"
    echo ">>> creating tracing_${n}"
    ${MAMBA} env create -y -f "${TMP}/tracing_${n}.yml"

    if [ "${n}" = "seurat" ]; then
        # DoubletFinder is GitHub-only; pin the commit used in the lab env (v2.0.4)
        ${MAMBA} run -n tracing_seurat Rscript -e \
          'remotes::install_github("chris-mcginnis-ucsf/DoubletFinder", ref = "03e9f37f891ef76a23cc55ea69f940c536ae8f9f", upgrade = "never")'
    fi
done
rm -rf "${TMP}"

cat <<'NOTE'

Done. Check each env, e.g.:
  mamba run -n tracing_scvi       python -c "import scanpy, scvi; print(scanpy.__version__, scvi.__version__)"
  mamba run -n tracing_cassiopeia python -c "import cassiopeia; print(cassiopeia.__version__)"
  mamba run -n tracing_seurat     Rscript -e "library(Seurat); library(DoubletFinder); library(DropletUtils)"
  mamba run -n tracing_cellbender cellbender --version

Cell Ranger 10.0.0 (step 01): download from 10x Genomics (free, needs a licence click):
  https://www.10xgenomics.com/support/software/cell-ranger/downloads
The GRCm39_tracing reference (mouse genome + GFP/mCherry/Cas9/Cre transgenes) and
refOct1.fa (recorder reference) are lab resources on c2b2 — ask Zejian for a copy;
they are deliberately not in this public repository.

To run the tutorial with your own envs instead of the lab ones, set in config/config.sh:
  ENVS=<your conda prefix>/envs   and rename the env folders used by the sbatch files
  (cellbender, ST, scvi, cassiopeia_env) — or edit the PATH line in each sbatch.
NOTE
