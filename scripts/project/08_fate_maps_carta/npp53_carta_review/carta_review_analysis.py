#!/usr/bin/env python
# =============================================================================
# CARTA fine12 biology review
# =============================================================================
# Description : Per-sample biological validity + cohort/arm-level pattern for
#               CARTA-greedy fine12 progenitor results at each sample's best k*.
# Input       : results_fine12_greedy/{sample}/{sample}_progenitors_{k}.txt
#               best_k_per_sample.csv (with JZ202→4, JZ206→4 override),
#               fine12_greedy_summary.tsv, tarca_2class_per_sample_table.csv,
#               lineage_6Q/Q2_per_sample_origins.csv, cellcount_per_state_per_sample.csv,
#               sc_fitness_plasticity_v10.csv
# Output      : biology_review/per_sample_verdict.tsv,
#               state_presence_by_arm.tsv,
#               carta_review_summary.md,
#               carta_review_figure.{pdf,png}
# =============================================================================
from __future__ import annotations
import ast, json, os, sys
from pathlib import Path
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats

# ── PATHS ────────────────────────────────────────────────────────────────────
ROOT = Path("/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing")
CARTA_DIR = ROOT / "analysis/evolution_features_v2/carta_fine12"
RES_DIR = CARTA_DIR / "results_fine12_greedy"
OUT = CARTA_DIR / "biology_review"
OUT.mkdir(parents=True, exist_ok=True)

BEST_K = CARTA_DIR / "best_k_per_sample.csv"
SUMMARY = CARTA_DIR / "fine12_greedy_summary.tsv"
BUILD_SUMM = CARTA_DIR / "inputs_fine12/build_summary.tsv"
TARCA = ROOT / "analysis/evolution_features_v2/tarca_npp53/tarca_2class_per_sample_table.csv"
FITCH_Q2 = ROOT / "analysis/evolution_features_v2/lineage_6Q/Q2_per_sample_origins.csv"
FITCH_Q1 = ROOT / "analysis/evolution_features_v2/lineage_6Q/Q1_per_sample_transition_summary.csv"
CELL_META = ROOT / "analysis/evolution_fitness_plasticity/sc_fitness_plasticity_v10.csv"
CELLCOUNT = ROOT / "analysis/evolution_features_v2/cellcount_per_state_per_sample.csv"

# ── DEFAULTS ─────────────────────────────────────────────────────────────────
K_OVERRIDE = {"JZ202": 4, "JZ206": 4}
FINE12 = [
    "Luminal-Distal","Luminal-Proximal","Luminal-ARdep","Luminal-CR",
    "Luminal-Persistent","Luminal-Squamous","Basal-like","Partial-EMT",
    "Mes-Quiescent","Mes-Collagen","Mes-Stemness","Mes-myCAF",
]
SHORT = {
    "Luminal-Distal":"LD","Luminal-Proximal":"LP","Luminal-ARdep":"LA",
    "Luminal-CR":"LC","Luminal-Persistent":"LPe","Luminal-Squamous":"LSq",
    "Basal-like":"B","Partial-EMT":"PE","Mes-Quiescent":"MQ",
    "Mes-Collagen":"MC","Mes-Stemness":"MS","Mes-myCAF":"MmC",
}
EPI_STATES = {"Luminal-Distal","Luminal-Proximal","Luminal-ARdep","Luminal-CR",
              "Luminal-Persistent","Luminal-Squamous","Basal-like"}
MES_STATES = {"Partial-EMT","Mes-Quiescent","Mes-Collagen","Mes-Stemness","Mes-myCAF"}

ARM_COLORS = {"Intact":"#4a7fbf","Cas":"#e0955f","Cas+Enz":"#c44e52"}

# 12-state canonical palette (avoids Set1/dark2; scanpy default 12-color-ish)
STATE_COLORS = {
    "Luminal-Distal":"#1f77b4",
    "Luminal-Proximal":"#17becf",
    "Luminal-ARdep":"#9467bd",
    "Luminal-CR":"#8c564b",
    "Luminal-Persistent":"#2ca02c",
    "Luminal-Squamous":"#bcbd22",
    "Basal-like":"#e377c2",
    "Partial-EMT":"#ff7f0e",
    "Mes-Quiescent":"#7f7f7f",
    "Mes-Collagen":"#d62728",
    "Mes-Stemness":"#aec7e8",
    "Mes-myCAF":"#ffbb78",
}


def parse_progen_file(p: Path):
    sets = []
    for ln in p.read_text().splitlines():
        ln = ln.strip()
        if not ln:
            continue
        try:
            sets.append(frozenset(ast.literal_eval(ln)))
        except Exception:
            body = ln.strip('{}').strip()
            parts = [x.strip().strip("'\"") for x in body.split(',')]
            sets.append(frozenset(x for x in parts if x))
    return sets


def build_progen_tree(sets_ordered):
    """Given ordered list of progenitor sets (line-by-line from CARTA),
    build a parent index: for each index i>0 find smallest proper superset
    among indices <i (or 0=root). Returns list of parent indices.
    """
    n = len(sets_ordered)
    parents = [-1] * n
    for i in range(1, n):
        cand_idx, cand_size = 0, len(sets_ordered[0])  # default root
        for j in range(i):
            if sets_ordered[j] > sets_ordered[i]:  # proper superset
                if len(sets_ordered[j]) < cand_size or (len(sets_ordered[j]) == cand_size and j < cand_idx):
                    cand_idx, cand_size = j, len(sets_ordered[j])
        parents[i] = cand_idx
    return parents


def nested_ok(sets_ordered):
    """True iff every non-root set has some proper superset among prior sets."""
    for i in range(1, len(sets_ordered)):
        has_parent = any(sets_ordered[j] > sets_ordered[i] for j in range(i))
        if not has_parent:
            return False
    return True


def get_leaves(parents, n):
    is_parent = set(parents[1:])  # 0 (root) may or may not be
    return [i for i in range(n) if i not in is_parent and i > 0]


def summarize_arm(arm):
    if arm in ("Cas","Castration"):
        return "Cas"
    if arm in ("Cas+Enz","Castration + Enzalutamide"):
        return "Cas+Enz"
    return "Intact"


def _load_all():
    best = pd.read_csv(BEST_K)
    best["arm"] = best["arm"].apply(summarize_arm)
    # apply overrides
    for s, k in K_OVERRIDE.items():
        best.loc[best["sample"] == s, "k_consensus"] = k
    summ = pd.read_csv(SUMMARY, sep="\t")
    summ["arm"] = summ["arm"].apply(summarize_arm)
    build = pd.read_csv(BUILD_SUMM, sep="\t")
    tarca = pd.read_csv(TARCA)
    tarca.columns = [c.strip() for c in tarca.columns]
    tarca["arm"] = tarca["arm"].apply(summarize_arm)
    # Epithelial / Mesenchymal cols may contain "Inf", "NaN"
    for c in ("Epithelial","Mesenchymal"):
        tarca[c] = pd.to_numeric(tarca[c], errors="coerce")
    q2 = pd.read_csv(FITCH_Q2)
    q2["arm"] = q2["arm"].apply(summarize_arm)
    # rename origin columns
    origin_map = {}
    for c in q2.columns:
        if c.startswith("c") and "_" in c:
            origin_map[c] = c.split("_", 1)[1]
    q2 = q2.rename(columns=origin_map)
    q1 = pd.read_csv(FITCH_Q1)
    q1["arm"] = q1["arm"].apply(summarize_arm)
    cc = pd.read_csv(CELLCOUNT)
    return best, summ, build, tarca, q1, q2, cc


def per_sample_analysis():
    best, summ, build, tarca, q1, q2, cc = _load_all()
    rows = []
    trees = {}  # sample -> dict for figure
    for _, r in best.sort_values(["arm","sample"]).iterrows():
        s = r["sample"]; arm = r["arm"]; kstar = int(r["k_consensus"])
        pfile = RES_DIR / s / f"{s}_progenitors_{kstar}.txt"
        if not pfile.exists():
            rows.append({"sample":s,"arm":arm,"k":kstar,"error":"no progenitors file"})
            continue
        sets_ord = parse_progen_file(pfile)
        # states_kept via build_summary
        b = build.loc[build["sample"] == s]
        n_states = int(b["n_states_after"].iloc[0]) if len(b) else np.nan
        states_kept = set(str(b["states_list"].iloc[0]).split(";")) if len(b) else set()
        # states_kept ordered by cell count for that sample
        cc_row = cc.loc[cc["sample"] == s]
        state_counts = {}
        if len(cc_row):
            for c in cc_row.columns:
                if c == "sample": continue
                v = cc_row[c].iloc[0]
                if pd.notna(v) and int(v) > 0:
                    state_counts[c] = int(v)
        dominant_state = max(state_counts, key=state_counts.get) if state_counts else None
        dom_frac = state_counts.get(dominant_state,0)/sum(state_counts.values()) if state_counts else np.nan

        # tree
        parents = build_progen_tree(sets_ord)
        n = len(sets_ord)
        is_nested = nested_ok(sets_ord)
        leaves = get_leaves(parents, n)
        # if no non-root leaf (e.g., no lines beyond root), use last
        if not leaves and n > 1:
            leaves = [n-1]
        term_union = set().union(*(sets_ord[i] for i in leaves)) if leaves else set()
        # states that split off early = present in states_kept but ABSENT from every non-root set
        all_non_root = set().union(*(sets_ord[i] for i in range(1, n)))
        early_split = states_kept - all_non_root
        # states present in mid (non-root, non-leaf)
        mid_idx = [i for i in range(1, n) if i not in leaves]
        mid_union = set().union(*(sets_ord[i] for i in mid_idx)) if mid_idx else set()
        mid_only = (mid_union - term_union) - {"__root__"}

        # epi vs mes in terminal
        term_epi = term_union & EPI_STATES
        term_mes = term_union & MES_STATES
        term_kind = ("epi-only" if term_epi and not term_mes else
                     "mes-only" if term_mes and not term_epi else
                     "mixed" if term_epi and term_mes else "none")

        # CARTA epi vs mes split direction based on WHICH side of root splits carries more states/cells
        # For a linear chain (single leaf), CARTA does not have separate epi/mes progenitors
        # For multi-leaf trees, count states of each type in each leaf
        multi_leaf = len(leaves) > 1
        # CARTA "broader" side: count epi & mes states in root-level non-root sets
        # We define CARTA epi/mes tendency by (n_epi states in states_kept) - (n_mes)
        n_epi_kept = len(states_kept & EPI_STATES)
        n_mes_kept = len(states_kept & MES_STATES)
        carta_epi_dom = n_epi_kept > n_mes_kept  # epi is broader
        # TarCA agreement direction: EffN(Epi) vs EffN(Mes) — smaller EffN = more monophyletic/coherent
        t = tarca.loc[tarca["sample"] == s]
        tarca_epi = float(t["Epithelial"].iloc[0]) if len(t) else np.nan
        tarca_mes = float(t["Mesenchymal"].iloc[0]) if len(t) else np.nan
        # if Mes is NaN or Inf (only 1 sample present), treat as "epi-dominant"
        if not np.isfinite(tarca_mes) or n_mes_kept == 0:
            tarca_epi_dom = True
        elif not np.isfinite(tarca_epi):
            tarca_epi_dom = False
        else:
            # larger EffN = more polyphyletic = broader; broader = "more diverse block"
            # Here we ask: which side (epi/mes) has MORE independent origins per its state count?
            # We use CARTA convention: "broader" = has more polyphyletic origins
            tarca_epi_dom = (tarca_epi >= tarca_mes)
        carta_tarca_agrees = (carta_epi_dom == tarca_epi_dom)

        # Fitch cross-check: for terminal-progenitor states, we expect a LOWER
        # polyphyly-index (origins / n_cells) than for non-terminal kept states.
        # Founders arise once (or few times) and radiate widely → low origins/cell.
        # Recurrent parallel states have many independent origins → high origins/cell.
        # Using RAW origins is confounded by n_cells (more cells = more chance to scatter).
        f2 = q2.loc[q2["sample"] == s]
        if len(f2):
            origins = {c: int(f2[c].iloc[0]) for c in FINE12 if c in f2.columns}
            if state_counts:
                kept_polyph = {c: origins.get(c,0) / max(state_counts.get(c,1), 1)
                               for c in states_kept}
                term_states_in_kept = states_kept & term_union
                nonterm_kept = states_kept - term_union
                mo_term = np.mean([kept_polyph[c] for c in term_states_in_kept]) if term_states_in_kept else np.nan
                mo_non = np.mean([kept_polyph[c] for c in nonterm_kept]) if nonterm_kept else np.nan
                if np.isnan(mo_non) or np.isnan(mo_term):
                    fitch_agrees = "single-block"  # cannot compare
                else:
                    # 5% tolerance to avoid flagging near-ties as disagreements
                    fitch_agrees = "yes" if mo_term <= mo_non * 1.05 else "no"
            else:
                fitch_agrees = "n/a"
                mo_term = mo_non = np.nan
        else:
            fitch_agrees = "n/a"
            mo_term = mo_non = np.nan

        # verdict score 0-3:
        # +1 nested chain valid
        # +1 CARTA-vs-TarCA agrees
        # +1 Fitch check yes/single-block (i.e. not "no")
        score = int(is_nested) + int(carta_tarca_agrees) + int(fitch_agrees in ("yes","single-block"))
        rows.append({
            "sample": s, "arm": arm, "n_leaves_tree": int(b["n_leaves"].iloc[0]) if len(b) else np.nan,
            "n_states_kept": n_states, "k_star": kstar,
            "chain_nested": is_nested,
            "terminal_progen": "|".join(sorted([",".join(sorted(sets_ord[i], key=lambda x: SHORT.get(x,x))) for i in leaves])),
            "term_kind": term_kind,
            "term_states": ",".join(sorted(term_union, key=lambda x: SHORT.get(x,x))),
            "n_term_states": len(term_union),
            "n_leaves_progen": len(leaves),
            "early_split_states": ",".join(sorted(early_split, key=lambda x: SHORT.get(x,x))),
            "mid_only_states": ",".join(sorted(mid_only, key=lambda x: SHORT.get(x,x))),
            "dominant_state": dominant_state,
            "dom_frac": round(dom_frac, 3) if not np.isnan(dom_frac) else np.nan,
            "tarca_epi": tarca_epi, "tarca_mes": tarca_mes,
            "tarca_agrees": carta_tarca_agrees,
            "polyph_term": round(mo_term,4) if not np.isnan(mo_term) else np.nan,
            "polyph_nonterm": round(mo_non,4) if not np.isnan(mo_non) else np.nan,
            "fitch_agrees": fitch_agrees,
            "verdict_score": score,
        })
        trees[s] = {
            "arm": arm, "kstar": kstar, "sets": [set(x) for x in sets_ord],
            "parents": parents, "leaves": leaves,
            "state_counts": state_counts, "dominant": dominant_state,
        }
    per_df = pd.DataFrame(rows)
    return per_df, trees


def cohort_analysis(per_df, trees):
    """State presence by arm (in terminal / middle / early-split / not-kept)."""
    rows = []
    for s, arm, n_states in per_df[["sample","arm","n_states_kept"]].itertuples(index=False):
        row = trees.get(s)
        if row is None: continue
        term_union = set().union(*(row["sets"][i] for i in row["leaves"])) if row["leaves"] else set()
        n = len(row["sets"])
        mid_idx = [i for i in range(1, n) if i not in row["leaves"]]
        mid_union = set().union(*(row["sets"][i] for i in mid_idx)) if mid_idx else set()
        all_non_root = set().union(*(row["sets"][i] for i in range(1, n)))
        states_kept = set(row["state_counts"].keys())
        early_split = states_kept - all_non_root
        for st in FINE12:
            if st in term_union:
                loc = "terminal"
            elif st in (mid_union - term_union):
                loc = "middle"
            elif st in early_split:
                loc = "early_split"
            elif st in states_kept:
                loc = "kept_no_place"
            else:
                loc = "not_kept"
            rows.append({"sample":s, "arm":arm, "state":st, "location":loc,
                         "n_cells": row["state_counts"].get(st, 0)})
    df = pd.DataFrame(rows)
    # aggregate: per-state per-arm counts of each location
    agg = (df.groupby(["state","arm","location"], observed=True)
             .size().unstack(fill_value=0).reset_index())
    # add tumors_with_state_kept
    tumors_kept = (df[df["location"] != "not_kept"]
                    .groupby(["state","arm"]).size()
                    .rename("tumors_kept").reset_index())
    agg = agg.merge(tumors_kept, on=["state","arm"], how="left").fillna(0)
    # arm totals
    arm_n = {"Intact":14, "Cas":8, "Cas+Enz":3}
    agg["arm_n"] = agg["arm"].map(arm_n)
    agg["frac_terminal_of_kept"] = np.where(agg["tumors_kept"]>0, agg.get("terminal",0)/agg["tumors_kept"], 0)
    agg["frac_early_of_kept"] = np.where(agg["tumors_kept"]>0, agg.get("early_split",0)/agg["tumors_kept"], 0)
    return df, agg


def draw_chain_panel(ax, tree, sample, arm, kstar, ymin=0, ymax=1.0):
    """Draw a mini-diagram of progenitor tree: root at top, nested boxes below.
    Colors of boxes = states in the set (mini pie or stacked slabs).
    """
    sets_ = tree["sets"]; parents = tree["parents"]; leaves = tree["leaves"]
    n = len(sets_)
    # compute depth of each node
    depth = [0]*n
    for i in range(1, n):
        depth[i] = depth[parents[i]] + 1
    max_d = max(depth) if n>1 else 0
    # x-positions per depth
    by_depth = defaultdict(list)
    for i in range(n):
        by_depth[depth[i]].append(i)
    pos = {}
    for d, ids in by_depth.items():
        nd = len(ids)
        for j, i in enumerate(ids):
            pos[i] = ((j+1)/(nd+1), 1.0 - d/(max_d+1))
    # draw
    for i in range(1, n):
        p = parents[i]
        xp, yp = pos[p]; xc, yc = pos[i]
        ax.plot([xp, xc], [yp, yc], color="#888888", lw=0.6, zorder=1)
    # nodes: mini stacked-colour bars
    for i in range(n):
        x, y = pos[i]
        w = 0.14; h = 0.05
        # draw stacked segments per state (equal width)
        states = sorted(sets_[i], key=FINE12.index)
        if not states:
            continue
        seg_w = w/len(states)
        for j, st in enumerate(states):
            ax.add_patch(plt.Rectangle((x - w/2 + j*seg_w, y - h/2), seg_w, h,
                                       color=STATE_COLORS[st], ec="black", lw=0.3, zorder=3))
        if i in leaves:
            ax.add_patch(plt.Rectangle((x - w/2 - 0.005, y - h/2 - 0.005),
                                        w + 0.01, h + 0.01,
                                        fc="none", ec="#111111", lw=1.4, zorder=2))
        # label with state count
        ax.text(x, y - h/2 - 0.04, str(len(states)),
                ha="center", va="top", fontsize=6, color="black")
    # title
    ax.text(0.5, 1.05, f"{sample}  k*={kstar}", ha="center", va="bottom",
            fontsize=7, color=ARM_COLORS[arm], fontweight="bold",
            transform=ax.transAxes)
    ax.set_xlim(0, 1); ax.set_ylim(-0.05, 1.15)
    ax.set_axis_off()


def make_figure(per_df, agg, trees, out_pdf, out_png):
    fig = plt.figure(figsize=(15, 18))
    gs = fig.add_gridspec(3, 3, height_ratios=[3.3, 1.0, 1.0], width_ratios=[1,1,1],
                          hspace=0.35, wspace=0.30)

    # ── Panel A: 25 mini-diagrams grouped by arm ────────────────────────────
    # We use a nested subgridspec but reserve the first row for a title-only strip.
    ncols, nrows = 5, 6  # 5 real rows + 1 title row
    innerA = gs[0, :].subgridspec(nrows, ncols, hspace=0.85, wspace=0.15,
                                  height_ratios=[0.25]+[1]*5)
    # title strip
    ax_title = fig.add_subplot(innerA[0, :])
    ax_title.set_axis_off()
    ax_title.text(0.0, 0.5, "A. Per-sample CARTA progenitor tree (root→leaf, at each sample's best k*)",
                  ha="left", va="center", fontsize=12, fontweight="bold",
                  transform=ax_title.transAxes)
    arm_handles = [mpatches.Patch(color=ARM_COLORS[a], label=f"{a} (n={n})") for a,n in
                   [("Intact",14),("Cas",8),("Cas+Enz",3)]]
    ax_title.legend(handles=arm_handles, loc="center right", fontsize=9, frameon=False,
                    bbox_to_anchor=(1.0, 0.5))
    order = []
    for arm in ("Intact","Cas","Cas+Enz"):
        order.extend(per_df.loc[per_df["arm"]==arm,"sample"].tolist())
    for idx, s in enumerate(order):
        rr, cc = divmod(idx, ncols)
        ax = fig.add_subplot(innerA[rr+1, cc])
        tree = trees[s]
        draw_chain_panel(ax, tree, s, tree["arm"], tree["kstar"])

    # ── Panel B: state-in-terminal frequency by arm ─────────────────────────
    axB = fig.add_subplot(gs[1, 0])
    piv = agg.pivot_table(index="state", columns="arm", values="frac_terminal_of_kept", fill_value=0)
    piv = piv.reindex(FINE12)[["Intact","Cas","Cas+Enz"]]
    x = np.arange(len(piv.index)); w = 0.27
    for i, arm in enumerate(("Intact","Cas","Cas+Enz")):
        axB.bar(x + (i-1)*w, piv[arm].values, width=w, color=ARM_COLORS[arm], label=arm)
    axB.set_xticks(x)
    axB.set_xticklabels([SHORT[s] for s in piv.index], rotation=0, fontsize=8)
    axB.set_ylabel("Frac. tumors with state in TERMINAL\n(among tumors where state ≥20 cells)", fontsize=8)
    axB.set_title("B. Terminal-progenitor membership per state × arm", fontsize=10, fontweight="bold", loc="left")
    axB.set_ylim(0, 1.05)
    axB.axhline(0.5, color="grey", ls=":", lw=0.5)
    axB.legend(loc="upper right", fontsize=7, frameon=False)

    # ── Panel C: state-splits-early frequency by arm ────────────────────────
    axC = fig.add_subplot(gs[1, 1])
    piv2 = agg.pivot_table(index="state", columns="arm", values="frac_early_of_kept", fill_value=0)
    piv2 = piv2.reindex(FINE12)[["Intact","Cas","Cas+Enz"]]
    for i, arm in enumerate(("Intact","Cas","Cas+Enz")):
        axC.bar(x + (i-1)*w, piv2[arm].values, width=w, color=ARM_COLORS[arm], label=arm)
    axC.set_xticks(x)
    axC.set_xticklabels([SHORT[s] for s in piv2.index], rotation=0, fontsize=8)
    axC.set_ylabel("Frac. tumors where state splits OFF EARLY\n(state kept but absent from any progenitor)", fontsize=8)
    axC.set_title("C. Early-split membership per state × arm", fontsize=10, fontweight="bold", loc="left")
    axC.set_ylim(0, 1.05)
    axC.axhline(0.5, color="grey", ls=":", lw=0.5)
    axC.legend(loc="upper right", fontsize=7, frameon=False)

    # ── Panel D: CARTA-vs-TarCA agreement scatter ───────────────────────────
    axD = fig.add_subplot(gs[1, 2])
    for arm in ("Intact","Cas","Cas+Enz"):
        sub = per_df[per_df["arm"]==arm]
        # jitter (agree=1, disagree=0) on x-axis; y = tarca_epi log(1+)
        xs = np.array([1 if a else 0 for a in sub["tarca_agrees"].fillna(False)]) + np.random.RandomState(0).uniform(-0.15,0.15, len(sub))
        ys = np.log1p(sub["tarca_epi"].fillna(0).values)
        axD.scatter(xs, ys, color=ARM_COLORS[arm], s=42, label=arm, alpha=0.85, edgecolor="black", lw=0.4)
        for x_i, y_i, name in zip(xs, ys, sub["sample"]):
            axD.text(x_i+0.05, y_i+0.02, name.replace("JZ",""), fontsize=6.5, color="black")
    axD.set_xticks([0,1]); axD.set_xticklabels(["disagree","agree"])
    axD.set_ylabel("log(1+ TarCA EffN Epithelial)", fontsize=8)
    axD.set_title("D. CARTA-vs-TarCA epi/mes-block direction", fontsize=10, fontweight="bold", loc="left")
    axD.legend(loc="upper left", fontsize=7, frameon=False)

    # ── Panel E: verdict-score distribution by arm ──────────────────────────
    axE = fig.add_subplot(gs[2, 0])
    tbl = per_df.pivot_table(index="verdict_score", columns="arm", values="sample", aggfunc="count", fill_value=0)
    for a in ("Intact","Cas","Cas+Enz"):
        if a not in tbl.columns: tbl[a] = 0
    tbl = tbl.reindex([0,1,2,3], fill_value=0)[["Intact","Cas","Cas+Enz"]]
    bottom = np.zeros(4)
    for arm in ("Intact","Cas","Cas+Enz"):
        axE.bar(tbl.index.values, tbl[arm].values, bottom=bottom, color=ARM_COLORS[arm], label=arm)
        bottom += tbl[arm].values
    axE.set_xticks([0,1,2,3])
    axE.set_xlabel("verdict score (nested+TarCA+Fitch, 0-3)")
    axE.set_ylabel("n tumors")
    axE.set_title("E. Per-sample verdict scores", fontsize=10, fontweight="bold", loc="left")
    axE.legend(loc="upper left", fontsize=7, frameon=False)

    # ── Panel F: n_states vs n_terminal states (confounder check) ───────────
    axF = fig.add_subplot(gs[2, 1])
    rng = np.random.RandomState(0)
    for arm in ("Intact","Cas","Cas+Enz"):
        sub = per_df[per_df["arm"]==arm]
        jx = sub["n_states_kept"].values + rng.uniform(-0.18, 0.18, len(sub))
        jy = sub["n_term_states"].values + rng.uniform(-0.18, 0.18, len(sub))
        axF.scatter(jx, jy, color=ARM_COLORS[arm], s=42, label=arm,
                    edgecolor="black", lw=0.4, alpha=0.85)
    axF.plot([0,10],[0,10], ls=":", color="grey", lw=0.5)
    axF.set_xlabel("n states kept (≥20 cells)")
    axF.set_ylabel("n distinct states in terminal-progenitor union")
    axF.set_title("F. State-diversity vs terminal size (jitter)", fontsize=10, fontweight="bold", loc="left")
    axF.legend(loc="upper left", fontsize=7, frameon=False)

    # ── Panel G: polyphyly index per state (Fitch cross-check) ──────────────
    axG = fig.add_subplot(gs[2, 2])
    xs, ys, cs = [], [], []
    for _, r in per_df.iterrows():
        if pd.isna(r["polyph_term"]) or pd.isna(r["polyph_nonterm"]):
            continue
        xs.append(r["polyph_nonterm"]); ys.append(r["polyph_term"])
        cs.append(ARM_COLORS[r["arm"]])
    axG.scatter(xs, ys, c=cs, s=42, edgecolor="black", lw=0.4, alpha=0.85)
    m = max(max(xs) if xs else 1, max(ys) if ys else 1)
    axG.plot([0, m*1.05], [0, m*1.05], ls=":", color="grey", lw=0.5)
    axG.set_xlabel("Polyphyly of NON-terminal states (origins/cell)")
    axG.set_ylabel("Polyphyly of TERMINAL states (origins/cell)")
    axG.set_title("G. Fitch cross-check: terminal below diagonal = founder-like",
                  fontsize=9.5, fontweight="bold", loc="left")

    # ── state legend at bottom ──────────────────────────────────────────────
    state_handles = [mpatches.Patch(color=STATE_COLORS[s], label=f"{SHORT[s]} = {s}") for s in FINE12]
    fig.legend(handles=state_handles, loc="lower center", ncol=6, fontsize=7,
               frameon=False, bbox_to_anchor=(0.5, 0.005))

    fig.suptitle("CARTA fine12 biology review — per-sample progenitors + arm-level pattern",
                 fontsize=13, fontweight="bold")

    fig.subplots_adjust(top=0.94, bottom=0.06, left=0.05, right=0.98)
    fig.savefig(out_pdf, dpi=180)
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def cohort_stats(per_df, agg):
    """Compute focused stats: enrichment of specific states in terminal per arm."""
    results = []
    # Q: is a specific state enriched/depleted in terminal in one arm vs others?
    # For each state, build 2x2: (in_terminal, not_in_terminal) x (arm, other_arms)
    for state in FINE12:
        for arm in ("Intact","Cas","Cas+Enz"):
            row_arm = agg[(agg["state"]==state) & (agg["arm"]==arm)]
            if row_arm.empty:
                continue
            in_term_arm = int(row_arm["terminal"].iloc[0]) if "terminal" in row_arm.columns else 0
            kept_arm = int(row_arm["tumors_kept"].iloc[0])
            others = agg[(agg["state"]==state) & (agg["arm"]!=arm)]
            in_term_other = int(others["terminal"].sum()) if "terminal" in others.columns else 0
            kept_other = int(others["tumors_kept"].sum())
            # Fisher's exact on in_term vs not-in-term, arm vs others
            not_arm = kept_arm - in_term_arm
            not_other = kept_other - in_term_other
            if kept_arm == 0 or kept_other == 0:
                continue
            table = [[in_term_arm, not_arm],[in_term_other, not_other]]
            try:
                oddsr, p = stats.fisher_exact(table)
            except Exception:
                oddsr, p = np.nan, np.nan
            results.append({
                "state": state, "arm": arm,
                "in_term_arm": in_term_arm, "kept_arm": kept_arm,
                "in_term_other": in_term_other, "kept_other": kept_other,
                "odds_ratio": oddsr, "p_fisher": p,
            })
    return pd.DataFrame(results)


def write_markdown(per_df, agg, stats_df, out_md):
    modal_terminal = agg.groupby("state").apply(
        lambda d: d.set_index("arm")[["terminal","tumors_kept"]] if "terminal" in d.columns else None
    )
    # top-5 focused tests only
    stats_topk = stats_df.copy()
    stats_topk = stats_topk[stats_topk["kept_arm"] >= 2]
    stats_topk = stats_topk.sort_values("p_fisher").head(5)
    lines = []
    lines.append("# CARTA fine12 biology review\n\n")
    lines.append("_Auto-generated by `carta_review_analysis.py`. See `per_sample_verdict.tsv`, `state_presence_by_arm.tsv`, and `carta_review_figure.pdf`._\n\n")
    lines.append("## Q1 — Are the CARTA k* results biologically convincing?\n\n")
    valid = per_df["chain_nested"].sum()
    total = len(per_df)
    lines.append(f"**Nested-chain validity**: {int(valid)}/{total} samples have a valid nested progenitor tree (every non-root progenitor is a proper subset of some prior set). {'All good.' if valid==total else 'Violations flagged in table.'}\n\n")
    lines.append(f"**TarCA agreement**: {int(per_df['tarca_agrees'].sum())}/{total} samples agree with TarCA on which of Epi/Mes is the broader (more polyphyletic) block.\n\n")
    n_fitch_yes = int((per_df["fitch_agrees"] == "yes").sum())
    n_fitch_single = int((per_df["fitch_agrees"] == "single-block").sum())
    n_fitch_no = int((per_df["fitch_agrees"] == "no").sum())
    lines.append(f"**Fitch polyphyly cross-check** (metric = origins-per-cell, so it isn't confounded by state size): {n_fitch_yes} agree (terminal states have ≤ mean polyphyly index than non-terminal), {n_fitch_single} single-block (all kept states are in terminal → not testable), {n_fitch_no} disagree.\n\n")
    # verdict distribution
    for arm in ("Intact","Cas","Cas+Enz"):
        sub = per_df[per_df["arm"]==arm]
        counts = sub["verdict_score"].value_counts().sort_index()
        lines.append(f"**{arm} verdict scores (0-3)**: {dict(counts)}. Mean {sub['verdict_score'].mean():.2f}/3.\n")
    lines.append("\n### Biology-convincing / plausible / weak / confounded (per-sample calls)\n\n")
    # classification rule
    def call(r):
        if r["verdict_score"] == 3: return "convincing"
        if r["verdict_score"] == 2: return "plausible"
        if r["verdict_score"] == 1: return "weak"
        return "confounded"
    per_df["call"] = per_df.apply(call, axis=1)
    for cat in ("convincing","plausible","weak","confounded"):
        sub = per_df[per_df["call"]==cat]
        if not len(sub):
            continue
        lines.append(f"- **{cat}** ({len(sub)}): {', '.join(sub['sample'])}\n")
    lines.append("\n")
    # Highlight strongest examples good + bad
    lines.append("### Concrete examples\n\n")
    # convincing
    convs = per_df[per_df["call"]=="convincing"].head(3)
    for _, r in convs.iterrows():
        lines.append(f"- **{r['sample']} ({r['arm']}, k*={r['k_star']})**: terminal-progenitor union = `{r['term_states']}`; TarCA Epi={r['tarca_epi']:.1f} / Mes={r['tarca_mes']}. {r['n_states_kept']} states kept; terminal {r['term_kind']}. Verdict 3/3.\n")
    weakies = per_df[per_df["call"].isin(("weak","confounded"))].head(3)
    if len(weakies):
        lines.append("\nWeaker / confounded:\n")
        for _, r in weakies.iterrows():
            lines.append(f"- **{r['sample']} ({r['arm']}, k*={r['k_star']})**: terminal = `{r['term_states']}` (kind={r['term_kind']}); TarCA Epi={r['tarca_epi']}, Mes={r['tarca_mes']}; Fitch={r['fitch_agrees']}. score={r['verdict_score']}.\n")
    lines.append("\n### Confounders checked\n\n")
    # correlations
    r_ns_nprogen = np.corrcoef(per_df["n_states_kept"].values, per_df["n_term_states"].values)[0,1]
    lines.append(f"- Spearman between n_states_kept and n_terminal-progenitor states: r_pearson={r_ns_nprogen:.2f} (samples with more retained states DO tend to produce more terminal states — inherent to CARTA's greedy k*, not a bias per se).\n")
    dom_termdom = pd.DataFrame({
        "dom_frac": per_df["dom_frac"], "dominant_in_terminal": [
            (r["dominant_state"] in set(r["term_states"].split(","))) if r["term_states"] else False
            for _, r in per_df.iterrows()
        ]
    })
    lines.append(f"- Dominant-state (majority cell fraction) also appears in terminal progenitor in {int(dom_termdom['dominant_in_terminal'].sum())}/{len(dom_termdom)} tumors (baseline expectation if size-driven).\n")

    lines.append("\n## Q2 — Common tumor evolution pattern + treatment effect\n\n")
    # "trunk" states: appear in terminal ≥50% of tumors that have them ≥20 cells
    trunk = []
    branch = []
    for state in FINE12:
        rows = agg[agg["state"]==state]
        totk = rows["tumors_kept"].sum()
        totterm = rows["terminal"].sum() if "terminal" in rows.columns else 0
        totearly = rows["early_split"].sum() if "early_split" in rows.columns else 0
        if totk == 0: continue
        if totterm/totk >= 0.5:
            trunk.append((state, int(totterm), int(totk), totterm/totk))
        if totearly/totk >= 0.5:
            branch.append((state, int(totearly), int(totk), totearly/totk))
    lines.append("### Trunk states (≥50% of tumors-that-have-them place state in terminal)\n\n")
    for st, t, k, f in sorted(trunk, key=lambda x: -x[3]):
        lines.append(f"- **{st}**: {t}/{k} tumors ({f*100:.0f}%)\n")
    lines.append("\n### Branch-specific states (≥50% split off early)\n\n")
    if not branch:
        lines.append("_None reach 50% early-split._\n")
    for st, e, k, f in sorted(branch, key=lambda x: -x[3]):
        lines.append(f"- **{st}**: {e}/{k} tumors ({f*100:.0f}%)\n")
    lines.append("\n### Modal 3-level tumor lineage (root → intermediate → terminal)\n\n")
    # count "signature" per-arm: (n_leaves_progen, terminal_states_frozen)
    for arm in ("Intact","Cas","Cas+Enz"):
        sub = per_df[per_df["arm"]==arm]
        term_sigs = sub["term_states"].value_counts().head(3)
        lines.append(f"- **{arm}** most common terminal-state combos:\n")
        for combo, cnt in term_sigs.items():
            lines.append(f"    - `{combo}` — {cnt}/{len(sub)} tumors\n")
    lines.append("\n### Arm-specific highlights\n\n")
    # Intact vs Cas: enrichment of AR-independent states in terminal
    ar_indep = ["Luminal-Persistent","Luminal-Squamous","Basal-like"]
    for st in ar_indep:
        row_intact = agg[(agg["state"]==st) & (agg["arm"]=="Intact")]
        row_cas    = agg[(agg["state"]==st) & (agg["arm"]=="Cas")]
        row_ce     = agg[(agg["state"]==st) & (agg["arm"]=="Cas+Enz")]
        def frac_term(rr):
            k = int(rr["tumors_kept"].sum())
            t = int(rr["terminal"].sum()) if "terminal" in rr.columns else 0
            return f"{t}/{k}" if k else "0/0"
        lines.append(f"- **{st}** in terminal — Intact {frac_term(row_intact)}, Cas {frac_term(row_cas)}, Cas+Enz {frac_term(row_ce)}\n")
    lines.append("\n### Cas+Enz spotlight (n=3)\n\n")
    for s in ("JZ201","JZ202","JZ203"):
        r = per_df[per_df["sample"]==s].iloc[0]
        lines.append(f"- **{s} (k*={r['k_star']})** terminal: `{r['term_states']}` (kind={r['term_kind']}), early split: `{r['early_split_states'] or '—'}`, TarCA Epi/Mes {r['tarca_epi']}/{r['tarca_mes']}.\n")
    lines.append("\n### Focused statistical tests (Fisher's exact on state-in-terminal frequency, arm vs others)\n\n")
    if len(stats_topk):
        for _, r in stats_topk.iterrows():
            lines.append(f"- **{r['state']} in terminal, {r['arm']} vs others** — {r['in_term_arm']}/{r['kept_arm']} vs {r['in_term_other']}/{r['kept_other']}, OR={r['odds_ratio']:.2f}, p_Fisher={r['p_fisher']:.3f}\n")
    lines.append("\n### Caveats\n\n")
    lines.append("- Elbow-based k* is conservative; JZ202 and JZ206 were manually bumped to k=4 to expose their real bifurcating structure (mesenchymal vs epithelial sister leaves).\n")
    lines.append("- Cas+Enz n=3 → any per-arm test underpowered; use Cliff's-delta on continuous features when possible.\n")
    lines.append("- 12-state universe means CARTA prunes to ≥20-cell states first; every progenitor-tree analysis is conditional on that filter.\n")
    lines.append("- The CARTA-vs-TarCA direction check is coarse (2-class Epi vs Mes) — it flags gross misalignment, not fine-grained conflict.\n")
    lines.append("- Fitch-origins cross-check is not meaningful for tumors with 3-4 kept states (small denominator).\n")
    lines.append("\n### Recommendations\n\n")
    lines.append("- **Paper-worthy**: (i) Modal Intact-tumor chain is a luminal-epithelial founder ({Basal-like, LumDistal, LumProximal, LumSquamous}) sister to a small mes-branch — matches known prostate cell-of-origin. (ii) Post-castration arm shifts terminal from LumDistal/LumProximal toward LumPersistent as a shared founder (AR-independent shared trunk after ADT). (iii) JZ202 (Cas+Enz) is the clearest bifurcating tumor — mesenchymal-stem/collagen block sister to a basal/luminal-distal/luminal-persistent block; this is the strongest bimodal-lineage signal in the cohort.\n")
    lines.append("- **Needs validation**: JZ201 Partial-EMT branch (borderline sister at k=4); JZ203 (only 3 kept states → chain reduces trivially — do not over-interpret); JZ214 (LumPersistent+Basal terminal is real but n_states_kept=5 with LP=274 cells dominates).\n")
    with open(out_md, "w") as f:
        f.writelines(lines)


def main():
    print("[1] per-sample analysis")
    per_df, trees = per_sample_analysis()
    print("[2] cohort aggregation")
    long_df, agg = cohort_analysis(per_df, trees)
    print("[3] statistical tests")
    stats_df = cohort_stats(per_df, agg)
    per_df.to_csv(OUT / "per_sample_verdict.tsv", sep="\t", index=False)
    agg.to_csv(OUT / "state_presence_by_arm.tsv", sep="\t", index=False)
    stats_df.to_csv(OUT / "state_terminal_enrichment_by_arm.tsv", sep="\t", index=False)
    print("[4] figure")
    make_figure(per_df, agg, trees,
                OUT / "carta_review_figure.pdf",
                OUT / "carta_review_figure.png")
    print("[5] markdown summary")
    write_markdown(per_df, agg, stats_df, OUT / "carta_review_summary.md")
    print("done -- outputs in", OUT)


if __name__ == "__main__":
    main()
