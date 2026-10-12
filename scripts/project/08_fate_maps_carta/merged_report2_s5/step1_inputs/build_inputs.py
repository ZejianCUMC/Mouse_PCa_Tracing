#!/usr/bin/env python
# =============================================================================
# Step 1: CARTA inputs for 47 tumors x {c10 included, c10 excluded}
# =============================================================================
# Description : Mirrors Npp53 carta_fine12/build_fine12_labels.py + prune_trees_fine12.py
#               (same >=20 cells per state and >=3 states per tumor rules, ete3 pruning with
#               branch lengths kept), but on the merged 10-state labels and the current trees
#               used by report 2 Appendix D. Cycling cells = c10_Mes-Proliferation (ZJ 2026-10-12);
#               the "exc" version drops every c10 cell before the per-state filter.
# Input       : reviewed10/fsu_m3_merged_20261008/cell_inputs_PRIVATE.csv (cell, sample_id, new_state, bc16),
#               .../metaNorm_20261011_all_cycout_bcr/tree_manifest.csv (current Greedy tree per tumor)
# Output      : step1_inputs/{inc,exc}/{T}_labels.tsv + {T}_tree_greedy.nwk, states_{inc,exc}.txt,
#               samples.txt, carta_tasks.tsv, state_crosswalk.csv, build_summary.tsv, md5sums.txt
#               One CARTA states file (all 10 states) serves both versions, as the original carried absent
#               states; where excluding c10 leaves the inputs byte-identical (c10 < 20 cells in that tumor),
#               the "exc" CARTA result is the "inc" result and no second run is queued.
# Conda env   : /home/jrchristin/miniconda3/envs/Cassiopeia/bin/python (ete3 3.1.3, pandas 2.2.2)
# Validated   : 2026-10-11
# Parameters  : See DEFAULTS dict below
# =============================================================================

# ── IMPORTS ──────────────────────────────────────────────────────────────────
import hashlib, json
from pathlib import Path
import pandas as pd
from ete3 import Tree

# ── DEFAULTS (all tunable parameters) ────────────────────────────────────────
T1 = Path("/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based")
RV10 = T1 / "Npp_tumor_merged/analysis/merged_20260923_seed42/restored_old_c18_seed42/recluster_r1_seed42/approved_merges/reviewed10"
DEFAULTS = {
    "cells": RV10 / "fsu_m3_merged_20261008/cell_inputs_PRIVATE.csv",          # labels + bc16 per cell (55,133 cells)
    "manifest": RV10 / "fsu_m3_merged_20261008/opt5_corrected_20261009/metaNorm_20261011_all_cycout_bcr/tree_manifest.csv",  # 47 current trees
    "cycling_state": "c10_Mes-Proliferation",   # ZJ 2026-10-12: all c10 = cycling cells, no re-identification
    "min_state_cells": 20,                      # per-tumor state filter, same as the original (>=20)
    "min_states": 3,                            # tumor kept only if >=3 states survive, same as the original
}
OUT = Path(__file__).resolve().parent

# ── HELPER FUNCTIONS ─────────────────────────────────────────────────────────
def carta_name(state):
    """CARTA-safe state name: parentheses dropped (c6_Partial-EMT(Lum) -> c6_Partial-EMT-Lum)."""
    return state.replace("(", "-").replace(")", "")

def bc16(name):
    return str(name).split("_")[-1].split("-")[0]

def md5(p):
    return hashlib.md5(Path(p).read_bytes()).hexdigest()

# ── MAIN PIPELINE FUNCTION ───────────────────────────────────────────────────
def run_pipeline(output_dir=OUT, **kwargs):
    P = {**DEFAULTS, **kwargs}
    cells = pd.read_csv(P["cells"], usecols=["cell", "model", "sample_id", "treatment", "new_state", "bc16", "LBI"])
    man = pd.read_csv(P["manifest"])
    states_all = sorted(cells.new_state.unique(), key=lambda s: int(s.split("_")[0][1:]))
    xw = pd.DataFrame({"new_state": states_all, "carta_state": [carta_name(s) for s in states_all]})
    xw.to_csv(output_dir / "state_crosswalk.csv", index=False)
    versions = {"inc": states_all, "exc": [s for s in states_all if s != P["cycling_state"]]}
    for v in versions:
        (output_dir / v).mkdir(exist_ok=True)
    (output_dir / "states.txt").write_text("\n".join(sorted(carta_name(s) for s in states_all)) + "\n")
    rows, samples = [], []
    for r in man.itertuples():
        t_id = r.tumour
        sub = cells[cells.sample_id == t_id]
        tree = Tree(str(T1 / r.tree), format=1)
        leaves = [l.name for l in tree.get_leaves()]
        lf = pd.DataFrame({"leaf": leaves, "bc16": [bc16(x) for x in leaves]})
        n_dup = int(lf.bc16.duplicated().sum())
        lf = lf.drop_duplicates("bc16")                       # same rule as 01_prepare.py
        m = lf.merge(sub[["bc16", "new_state", "model", "treatment", "LBI"]], on="bc16", how="inner")
        assert m.bc16.is_unique, t_id
        n_lbi = int(sub.LBI.notna().sum())                    # Appendix D tree-matched cells, must agree
        samples.append(t_id)
        for v, sts in versions.items():
            mv = m[m.new_state.isin(sts)]
            cnt = mv.new_state.value_counts()
            keep = sorted(cnt[cnt >= P["min_state_cells"]].index, key=states_all.index)
            mk = mv[mv.new_state.isin(keep)]
            ok = len(keep) >= P["min_states"]
            if ok:
                lab = pd.DataFrame({"cellBC": mk.leaf, "cell_state": mk.new_state.map(carta_name)})
                lab.to_csv(output_dir / v / f"{t_id}_labels.tsv", sep="\t", index=False)
                tp = tree.copy("newick")
                tp.prune(list(mk.leaf), preserve_branch_length=True)
                (output_dir / v / f"{t_id}_tree_greedy.nwk").write_text(tp.write(format=1))
                chk = {l.name for l in Tree(str(output_dir / v / f"{t_id}_tree_greedy.nwk"), format=1).get_leaves()}
                assert chk == set(mk.leaf), (t_id, v)
            rows.append(dict(tumor=t_id, model=m.model.iloc[0], treatment=m.treatment.iloc[0], version=v,
                             tree=r.tree, tree_leaves=len(leaves), dup_bc16_leaves=n_dup, matched=len(m),
                             appendixD_matched=n_lbi, c10_cells=int((m.new_state == P["cycling_state"]).sum()),
                             cells_kept=len(mk), n_states=len(keep), runnable=ok,
                             states=";".join(keep), counts=json.dumps({k: int(cnt[k]) for k in keep})))
    (output_dir / "samples.txt").write_text("\n".join(samples) + "\n")
    summ = pd.DataFrame(rows)
    tasks = []
    for t_id in samples:
        tasks.append(("inc", t_id))
        same = all(md5(output_dir / "inc" / f) == md5(output_dir / "exc" / f) for f in (f"{t_id}_labels.tsv", f"{t_id}_tree_greedy.nwk"))
        summ.loc[(summ.tumor == t_id) & (summ.version == "exc"), "exc_same_as_inc"] = same
        if not same: tasks.append(("exc", t_id))
    pd.DataFrame(tasks, columns=["version", "tumor"]).to_csv(output_dir / "carta_tasks.tsv", sep="\t", index=False)
    summ.to_csv(output_dir / "build_summary.tsv", sep="\t", index=False)
    files = sorted(p for p in output_dir.rglob("*") if p.is_file() and p.suffix in (".tsv", ".nwk", ".txt", ".csv") and p.name != "md5sums.txt")
    (output_dir / "md5sums.txt").write_text("".join(f"{md5(p)}  {p.relative_to(output_dir)}\n" for p in files))
    return {"output_path": str(output_dir), "params_used": {k: str(v) for k, v in P.items()},
            "summary": summ.groupby(["model", "version"]).agg(tumors=("tumor", "size"), runnable=("runnable", "sum"),
                                                              cells=("cells_kept", "sum"), c10=("c10_cells", "sum"))}

# ── CLI ENTRY POINT ──────────────────────────────────────────────────────────
if __name__ == "__main__":
    res = run_pipeline()
    s = pd.read_csv(OUT / "build_summary.tsv", sep="\t")
    print(res["summary"])
    print("matched == Appendix D tree-matched cells for every tumor:", bool((s.matched == s.appendixD_matched).all()))
    print("duplicate-barcode leaves:", int(s.drop_duplicates("tumor").dup_bc16_leaves.sum()))
    print(s.pivot(index="tumor", columns="version", values="n_states").to_string())
