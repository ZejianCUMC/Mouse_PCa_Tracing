# =============================================================================
# JZ201-205 trees from the MJZ019 re-run allele table (min_reads_per_umi=1)
# =============================================================================
# Description : Per-sample tree inputs in the v8 per_sample format (same valid tumour cells as v8; error intBCs removed
#               with the <=2-edit / <=30%-cells rule), then Greedy trees with the exact v8 call
#               (convert_alleletable_to_character_matrix(allele_rep_thresh=0.99) on Ionocyte/Cycling-free rows,
#               VanillaGreedySolver). Original v8 trees untouched.
# Input       : cassiopeia_MJZ019_rerun_minreads1_20261004/fast_faithful/MJZ019/MJZ019_allele_table.csv (md5 = main job)
# Output      : inputs/<S>_allele_table.csv, greedy/<S>_tree_greedy.nwk, INPUT_SUMMARY.csv
# Conda env   : /home/jrchristin/miniconda3/envs/Cassiopeia_v2/bin/python (local Cassiopeia 2.0; v8 used HPC 2.1)
# Validated   : 2026-10-05
# =============================================================================
import sys, os, time
from pathlib import Path
import numpy as np, pandas as pd
T = Path("/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing")
N = T / "1-mouse_based/Npp53_tracing"; OUT = N / "lineage_tracing/MJZ019_rerun_trees_20261005"
DEFAULTS = dict(rerun=T / "cassiopeia_MJZ019_rerun_minreads1_20261004/fast_faithful/MJZ019/MJZ019_allele_table.csv",
                obs=N / "data/npp53_merge/obs_metadata.csv", samples=["JZ201", "JZ202", "JZ203", "JZ204", "JZ205"],
                exclude={"Ionocyte", "Cycling"}, max_lev=2, max_ratio=0.30)
def lev(a, b):
    d = list(range(len(b) + 1))
    for i, x in enumerate(a, 1):
        prev, d[0] = d[0], i
        for j, y in enumerate(b, 1): prev, d[j] = d[j], min(d[j] + 1, d[j - 1] + 1, prev + (x != y))
    return d[-1]
def errors(a, p):
    det = a.groupby("intBC").cellBC.nunique().sort_values(ascending=False); out = {}
    for g, n in det.items():
        c = [h for h, m in det.items() if m > n and n <= p["max_ratio"] * m and lev(g, h) <= p["max_lev"]]
        if c: out[g] = c[0]
    return out
def main(stage):
    p = DEFAULTS
    if stage == "inputs":
        (OUT / "inputs").mkdir(parents=True, exist_ok=True)
        obs = pd.read_csv(p["obs"], index_col=0); obs = obs[obs.run == "MJZ019"]
        obs["bc"] = obs.index.str.replace(r"^MJZ019_", "", regex=True).str.replace(r"-\d+$", "", regex=True)
        a = pd.read_csv(p["rerun"], dtype=str, keep_default_na=False); rows = []
        for s in p["samples"]:
            o = obs[obs.sample_id == s]; m = dict(zip(o.bc, o.index))
            x = a[a.cellBC.isin(set(o.bc))].copy()
            x["cell_full"] = x.cellBC.map(m); x["sample_id"] = s; x["run"] = "MJZ019"
            x["treatment"] = x.cell_full.map(o.treatment); x["putative_state"] = x.cell_full.map(o.putative_state)
            e = errors(x[~x.putative_state.isin(p["exclude"])], p)
            y = x[~x.intBC.isin(set(e))]
            f = OUT / "inputs" / f"{s}_allele_table.csv"; assert not f.exists(); y.to_csv(f, index=False)
            rows.append(dict(sample=s, rna_cells=len(o), cells=y.cellBC.nunique(), rows=len(y), intBCs=y.intBC.nunique(), error_intBCs_removed=len(e)))
            print(rows[-1], flush=True)
        pd.DataFrame(rows).to_csv(OUT / "INPUT_SUMMARY.csv", index=False)
    else:
        import cassiopeia as cas
        (OUT / "greedy").mkdir(exist_ok=True)
        for s in p["samples"]:
            f = OUT / "greedy" / f"{s}_tree_greedy.nwk"; assert not f.exists()
            at = pd.read_csv(OUT / "inputs" / f"{s}_allele_table.csv")          # default read, as v8
            at = at[~at.putative_state.isin(p["exclude"])]
            t0 = time.time(); cm, priors, _ = cas.pp.convert_alleletable_to_character_matrix(at, allele_rep_thresh=0.99)
            tree = cas.data.CassiopeiaTree(character_matrix=cm, priors=priors)
            cas.solver.VanillaGreedySolver().solve(tree, logfile="/dev/null")
            f.write_text(tree.get_newick(record_branch_lengths=True)); print(s, cm.shape, f"{time.time()-t0:.0f}s", flush=True)
if __name__ == "__main__": main(sys.argv[1])
