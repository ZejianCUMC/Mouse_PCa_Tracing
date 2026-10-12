#!/usr/bin/env python
# =============================================================================
# Step 4: per-cell PhyloVelo pseudotime for the exact current tree of each tumor
# =============================================================================
# Description : Npp53: <tree folder>/phylovelo_out/{T}_phylovelo_pseudotime.csv (UMAP-based run, the same file
#               family build_combined_pseudotime.py used). Npp53Rb1: <viewer results>/{T}/PhyloVelo/pseudotime.csv
#               (written only when valid; UMAP-based). Folder per tumor mirrors the LBI source map in
#               reviewed10/survival_signature_proposal_20261008/check_fitness_mapping.py, i.e. the same tree.
#               Validity rule as the original: >= 20 distinct values after rounding to 4 decimals.
#               JZ270: the run with its one UMAP-isolated cell excluded (PhyloVelo_excl, approved 2026-09-25).
#               Fallback (as build_combined_pseudotime.py did): if the current-tree run is degenerate, use the older
#               cassiopeia_phylovelo_out/{T}_pseudotime_scaled.newick (root-to-tip distance) when it is valid.
# Input       : step3_c10_tree_position/labels_all/{T}.tsv (leaf, new_state), pseudotime CSVs above
# Output      : pseudotime_percell.csv (tumor, model, treatment, leaf, state, pt_raw), pseudotime_sources.tsv
# Conda env   : /home/jrchristin/miniconda3/envs/Cassiopeia/bin/python (pandas + ete3)
# Validated   : 2026-10-11
# =============================================================================
from pathlib import Path
import pandas as pd
from ete3 import Tree
T1 = Path("/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based")
N = T1 / "Npp53_tracing/lineage_tracing"; R = T1 / "Npp53Rb1_Tracing/lineage_tracing"
HERE = Path(__file__).resolve().parent
LAB = HERE.parent / "step3_c10_tree_position/labels_all"
def src(t, model):
    if model == "Npp53":
        if t == "JZ136": return N / "JZ136_intbcfix_20261004/phylovelo_out/JZ136_phylovelo_pseudotime.csv"
        if t == "JZ201": return N / "JZ201_excl_outlier_20261006/phylovelo_out/JZ201_phylovelo_pseudotime.csv"
        if t in ("JZ202", "JZ203", "JZ204", "JZ205"): return N / f"MJZ019_rerun_trees_20261005/phylovelo_out/{t}_phylovelo_pseudotime.csv"
        return N / f"cassiopeia_trees_v8_rebuilt/phylovelo_out/{t}_phylovelo_pseudotime.csv"
    tag = "JZ135_n3000" if t == "JZ135" else t
    if tag in ("JZ137", "JZ138", "JZ135_n3000"): return R / f"intbcfix_20261004_seed42/results/{tag}/PhyloVelo/pseudotime.csv"
    if tag in ("JZ256", "JZ259", "JZ261", "JZ262", "JZ263"): return R / f"tree_viewer_20260923_seed42/results/{tag}/PhyloVelo/pseudotime.csv"
    if tag in ("JZ267", "JZ269", "JZ274"): return R / f"state_min20_rebuild_20261006/{tag}/viewer/results/{tag}/PhyloVelo/pseudotime.csv"
    if tag == "JZ272": return R / "JZ272_state_min20_test_20261006/viewer/results/JZ272/PhyloVelo/pseudotime.csv"
    if tag == "JZ270": return R / "tree_viewer_mjz026_028_20260925_seed42/results/JZ270/PhyloVelo_excl/pseudotime.csv"
    return R / f"tree_viewer_mjz026_028_20260925_seed42/results/{tag}/PhyloVelo/pseudotime.csv"
def old_newick(t):
    f = N / f"cassiopeia_phylovelo_out/{t}_pseudotime_scaled.newick"
    if not f.exists(): return None, f
    tr = Tree(str(f), format=1)
    return pd.Series({l.name: tr.get_distance(l) for l in tr.get_leaves()}, name="pt_raw").to_frame(), f
bc16 = lambda x: str(x).split("_")[-1].split("-")[0]
tum = pd.read_csv(HERE.parent / "step3_c10_tree_position/tumors.tsv", sep="\t")
recs, rows = [], []
for r in tum.itertuples():
    f = src(r.tumor, r.model)
    lab = pd.read_csv(LAB / f"{r.tumor}.tsv", sep="\t"); lab["bc16"] = lab.leaf.map(bc16)
    row = dict(tumor=r.tumor, model=r.model, source=str(f.relative_to(T1)), exists=f.exists())
    pt = pd.read_csv(f, index_col=0).iloc[:, 0].rename("pt_raw").to_frame() if f.exists() else None
    if pt is not None and pt.pt_raw.round(4).nunique() < 20 and r.model == "Npp53":
        alt, fa = old_newick(r.tumor)
        if alt is not None and alt.pt_raw.round(4).nunique() >= 20:
            pt = alt; row.update(source=str(fa.relative_to(T1)), fallback=True)
    if pt is not None:
        pt["bc16"] = pt.index.map(bc16); pt = pt.drop_duplicates("bc16")
        m = lab.merge(pt, on="bc16", how="left")
        nun = int(m.pt_raw.round(4).nunique())
        row.update(n_tree_cells=len(lab), n_with_pt=int(m.pt_raw.notna().sum()), n_distinct_4dp=nun, valid=nun >= 20)
        m = m[m.pt_raw.notna()].assign(tumor=r.tumor, model=r.model, treatment=r.treatment)
        recs.append(m[["tumor", "model", "treatment", "leaf", "new_state", "pt_raw"]].rename(columns={"new_state": "state"}))
    rows.append(row)
src_tab = pd.DataFrame(rows); src_tab.to_csv(HERE / "pseudotime_sources.tsv", sep="\t", index=False)
pd.concat(recs).to_csv(HERE / "pseudotime_percell.csv", index=False)
print(src_tab.to_string(index=False))
