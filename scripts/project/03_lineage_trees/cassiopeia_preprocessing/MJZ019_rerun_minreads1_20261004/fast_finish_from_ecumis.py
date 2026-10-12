#!/usr/bin/env python
"""[MJZ019 rerun, fast branch 2026-10-04] Steps 9-10 of run_cassiopeia_pp_minreads1.py, started from the main run's
error_correct_umis checkpoint, with the two validated changes of run_cassiopeia_pp_faithful_wrapper.py
(byte-identical to an uninterrupted run on MJZ027, job 713908):
 (1) checkpoint read with keep_default_na=False (empty allele/r1-r3 stay '' as in memory);
 (2) error_correct_intbc replaced by ec_intbc_faithful (original ngs_tools N-aware distance, original order, parallel).
All filter/lineage parameters identical to run_cassiopeia_pp_minreads1.py (min_reads_per_umi=1).
Usage: fast_finish_from_ecumis.py LIB ECUMIS_CSV OUTDIR NTHREADS"""
import os, sys, time, functools
import pandas as pd
import cassiopeia as cas
from cassiopeia.preprocess import utilities as _u
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ec_intbc_faithful import error_correct_intbc_faithful
LIB, ECUMIS, OUTDIR, NTHREADS = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
_u.error_correct_intbc = functools.partial(error_correct_intbc_faithful, n_jobs=NTHREADS)
os.makedirs(OUTDIR, exist_ok=True)
def ck(n): return os.path.join(OUTDIR, f"{LIB}_{n}.csv")
def log(m): print(f"[{time.strftime('%F %T')}] {LIB}: {m}", flush=True)
assert not os.path.exists(ck("allele_table")), "refusing to overwrite"
cas.pp.setup(OUTDIR, verbose=True)
log(f"load {ECUMIS} (keep_default_na=False)"); umi = pd.read_csv(ECUMIS, keep_default_na=False); log(f"loaded {umi.shape}")
umi = cas.pp.filter_molecule_table(umi, output_directory=OUTDIR,
          min_umi_per_cell=10, min_avg_reads_per_umi=2.0, min_reads_per_umi=1,
          intbc_prop_thresh=0.5, intbc_umi_thresh=10, intbc_dist_thresh=1,
          doublet_threshold=0.35, allow_allele_conflicts=False, plot=False)
umi.to_csv(ck("filtered_molecule_table"), index=False); log(f"filter_molecule_table done  {umi.shape}")
allele = cas.pp.call_lineage_groups(umi, output_directory=OUTDIR,
             min_umi_per_cell=10, min_avg_reads_per_umi=2.0, min_cluster_prop=0.005,
             min_intbc_thresh=0.05, inter_doublet_threshold=0.35, kinship_thresh=0.25, plot=False)
allele.to_csv(ck("allele_table"), index=False)
log(f"DONE  allele_table {allele.shape}  cells={allele['cellBC'].nunique()}  intBCs={allele['intBC'].nunique()}")
