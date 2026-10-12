"""Honor seed42-only instruction; preserve, but exclude, superseded computations."""
from pathlib import Path
import json,shutil
import numpy as np,pandas as pd
p=Path(__file__).resolve().parents[1];d=p/'analysis/merged_20260923_seed42/restored_old_c18_seed42/recluster_r1_seed42/approved_merges/reviewed10/pathway_statistics'
r=json.loads((d/'STARTED.json').read_text());meta=pd.read_csv(d/'metadata_PRIVATE.csv',index_col=0)
x=np.load(d/'aucell_seed42_PRIVATE.npy',mmap_mode='r');s=pd.read_csv(d/'statistics_seed42.csv');paired=pd.read_csv(d/'sample_paired_scores_seed42.csv')
assert x.shape==(55133,1619) and np.isfinite(x).all() and x.min()>=0 and x.max()<=1.000001
assert len(s)==16190 and set(s.seed)=={42} and set(paired.seed)=={42}
assert s.groupby('cluster').signature.nunique().eq(1619).all()
assert json.loads((d/'R_AUCELL_VERIFIED.json').read_text())['status']=='PASS'
unused=d/'unused_before_seed42_instruction';unused.mkdir(exist_ok=True);moved=[]
for seed in [19,888]:
 for stem in [f'aucell_seed{seed}_PRIVATE.npy',f'statistics_seed{seed}.csv',f'sample_paired_scores_seed{seed}.csv',f'sample_balanced_auc_seed{seed}.csv']:
  if (d/stem).exists():
   assert not (unused/stem).exists();shutil.move(d/stem,unused/stem);moved.append(stem)
r.update(status='PASS',seeds=[42],full_cohort_rescored=True,seed42_only_user_instruction=True,extra_seed_run_interrupted=True,unused_outputs_preserved=moved,multiple_testing='BH across all 10 states x all 1619 pathways for seed42, including untestable p=1',selection_before_results='seed42 q<0.05, positive median delta>=0.005 and >=70% positive samples; per-collection cap and Jaccard redundancy filter')
(d/'SCORING_COMPLETE.json').write_text(json.dumps(r,indent=2));print('PASS: seed42-only results validated; extra-seed outputs preserved outside active outputs')
