#!/bin/bash
# opt5 meta-analysis rounds (Zejian 2026-10-09: "submit these tests; if still no OS-related module, keep optimizing and re-submit").
# Usage: bash run_meta_round.sh <arm>   (the arm folder must already hold fsu_lists.json from 22_meta_lists.py)
# 04 Hotspot inputs per cell set in the lists -> upload -> c2b2 Hotspot array (05h_hotspot_array.sbatch: 5 CPU / 4 GB / 30 min per
# task; HPC jobs < 10 h are pre-approved) -> fetch -> 05c profile modules, 06 candidates, km_scan (SU2C + WCDT merged / JCO),
# 08 select, 09 random-set null, 14 proliferation correlation, 16 Zejian's rule (+ Castration / Cas+ENZA), 17a/17b figures.
set -euo pipefail
cd "$(dirname "$0")"; ARM=$1
PY=/home/jrchristin/zejian/miniforge3/envs/moscot/bin/python; RS=/home/jrchristin/miniconda3/envs/rMAST/bin/Rscript
H=/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp_tumor_merged/opt5_corrected_20261009
export PYTHONHASHSEED=42 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NCORES=12
CS=$($PY -c "import json; print(' '.join(sorted({d['cells'] for d in json.load(open('$ARM/fsu_lists.json')).values()})))")
for cs in $CS; do $PY 04_hotspot_inputs.py $ARM $cs; done
KS=${KS:-30}; TOPS=${TOPS:-500}
$PY -c "import json; open('hotspot_tasks_$ARM.txt','w').write(''.join(f'$ARM {l} {k} {t}\n' for l in json.load(open('$ARM/fsu_lists.json')) for k in '$KS'.split() for t in '$TOPS'.split()))"
N=$(wc -l < hotspot_tasks_$ARM.txt); echo "$(date +%T) $ARM: $N Hotspot tasks; cell sets: $CS"
ssh -o BatchMode=yes c2b2 "mkdir -p $H/$ARM/hotspot_inputs $H/$ARM/hotspot_out $H/logs"
scp -q 05_hotspot_run.py 05h_hotspot_array.sbatch hotspot_tasks_$ARM.txt c2b2:$H/; scp -q $ARM/fsu_lists.json c2b2:$H/$ARM/   # scripts re-uploaded every round
for cs in $CS; do scp -q $ARM/hotspot_inputs/${cs}_* c2b2:$H/$ARM/hotspot_inputs/; done
JOB=$(ssh -o BatchMode=yes c2b2 "cd $H && sbatch --parsable --export=ALL,TASKS=hotspot_tasks_$ARM.txt --array=1-$N%80 05h_hotspot_array.sbatch")
echo "$(date +%T) submitted c2b2 job $JOB"
until ssh -o BatchMode=yes c2b2 "sacct -j $JOB -X -n -o State --parsable2 | grep -qE 'PENDING|RUNNING'" 2>/dev/null; [ $? -eq 1 ]; do sleep 30; done
echo "$(date +%T) job $JOB: $(ssh -o BatchMode=yes c2b2 "sacct -j $JOB -X -n -o State --parsable2 | sort | uniq -c | tr '\n' ' '")"
mkdir -p $ARM/hotspot_out; scp -q -r "c2b2:$H/$ARM/hotspot_out/*" $ARM/hotspot_out/; echo "fetched $(ls $ARM/hotspot_out/*_k30.json | wc -l) Hotspot outputs"
$PY 05c_profile_modules.py $ARM; $PY 06_candidates.py $ARM
$RS km_scan.R $ARM/candidates_human.json $ARM/km_scan.csv 1 sel; $PY 08_select.py $ARM; $RS 09_null.R $ARM; $RS 14_prolif_corr.R $ARM
$PY 16_user_rule.py $ARM
if [ "$($PY -c "import pandas as pd; print(len(pd.read_csv('$ARM/user_rule_passers.csv')))" 2>/dev/null || echo 0)" != "0" ]; then $PY 17a_mouse_scores.py $ARM && $RS 17b_plot_user_rule.R $ARM; fi
echo "$(date +%T) ROUND $ARM DONE"
