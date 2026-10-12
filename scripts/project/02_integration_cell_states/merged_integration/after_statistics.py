"""Queue descriptive evidence once statistics succeed; never select states."""
import json
import subprocess
import sys
import time
from pathlib import Path

project = Path(__file__).resolve().parents[1]
out = project / 'analysis/merged_20260923_seed42'
states = out / 'states_seed42'
scripts = project / 'scripts'
python = '/home/jrchristin/miniconda3/envs/ST/bin/python'
status = out / 'DEPENDENT_EVIDENCE_STATUS.json'
assert not status.exists(), 'Refuse duplicate evidence worker'


def record(stage, **kw):
    value = {'stage': stage, 'time_unix': time.time(), **kw}
    status.write_text(json.dumps(value, indent=2))
    print(value, flush=True)


def run(name, command):
    subprocess.run([python, str(scripts / 'run_local.py'), '--log',
                    str(out / 'logs' / f'{name}.log'), '--hours', '.5',
                    *command], check=True)


record('WAITING_FOR_ALL_FIVE_STATISTICAL_TESTS')
deadline = time.monotonic() + 16 * 3600
while True:
    if time.monotonic() > deadline:
        record('STOPPED_WAIT_DEADLINE')
        sys.exit(1)
    try:
        upstream = json.loads((out / 'DEPENDENT_STATS_STATUS.json').read_text())
    except json.JSONDecodeError:
        # The upstream status writer may be between truncate and write.
        # Re-observing its status is not retrying a scientific job.
        time.sleep(1)
        continue
    stage = upstream['stage']
    if stage == 'SCSHC_TOP5_COMPLETE_REQUIRES_BIOLOGICAL_REVIEW':
        break
    if stage.startswith('STOPPED') or stage == 'NEEDS_RESOLUTION_RANGE_REVIEW':
        record('STOPPED_FOR_UPSTREAM_REVIEW', upstream=upstream)
        sys.exit(0)
    time.sleep(30)

try:
    record('AUDITING_CANDIDATE_PARTITIONS')
    run('review_candidates', [python, str(scripts / 'review_candidates.py'), str(states)])
    import pandas as pd
    rank = pd.read_csv(states / 'scSHC_candidates.csv')
    for row in rank.itertuples():
        tag = f'r{row.resolution:.1f}'
        record('PROFILING_FULL_GENE_EVIDENCE', candidate=tag, AMI_rank=int(row.rank))
        run(f'profile_scshc_{tag}', [python, str(scripts / 'profile_states.py'),
            str(out), str(states / f'scSHC_{tag}_labels_PRIVATE.csv'),
            'scSHC_cluster', f'scshc_{tag}'])
    record('EVIDENCE_READY_REQUIRES_BIOLOGICAL_REVIEW')
except Exception as error:
    record('STOPPED_EVIDENCE_FAILURE', error=str(error))
    raise
