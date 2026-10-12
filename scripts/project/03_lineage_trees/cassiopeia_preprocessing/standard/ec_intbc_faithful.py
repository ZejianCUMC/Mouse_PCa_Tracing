"""Faithful fast re-implementation of cassiopeia.preprocess.utilities.error_correct_intbc (v2.1.0).
Goal: byte-identical output to the original, only faster.
 - same grouping/aggregation/sort calls as the original
 - distance = the ORIGINAL ngs_tools.sequence.levenshtein_distance (IUPAC/N-aware NW alignment)
 - identical pair order; condition logic identical (distance evaluated only when UMI/prop pass: same AND result)
 - corrections replayed in the original order on an array -> last write wins exactly like sequential .loc writes
 - outer (cellBC, allele) groups processed in parallel (order-preserving joblib)"""
import numpy as np, pandas as pd, warnings

def _group_corrections(key, intbcs, umis, prop, umi_thresh, dist_thresh):
    import ngs_tools as ngs
    cellBC, allele = key
    out = []
    n = len(intbcs)
    for i1 in range(n):
        b1, u1 = intbcs[i1], umis[i1]
        for i2 in range(i1 + 1, n):
            b2, u2 = intbcs[i2], umis[i2]
            if not (u2 <= umi_thresh):
                continue
            if not ((u2 / (u1 + u2)) < prop):
                continue
            if ngs.sequence.levenshtein_distance(b1, b2) <= dist_thresh:
                out.append(((cellBC, b2, allele), b1))
    return out

def _batch(items, prop, umi_thresh, dist_thresh):
    return [_group_corrections(k, i, u, prop, umi_thresh, dist_thresh) for k, i, u in items]

def error_correct_intbc_faithful(molecule_table, prop=0.5, umi_count_thresh=10, dist_thresh=1, n_jobs=16, verbose=True):
    from joblib import Parallel, delayed
    if prop > 0.5:
        warnings.warn("No intBC correction was done because `prop` is greater than 0.5.")
        return molecule_table
    g3 = molecule_table.groupby(["cellBC", "intBC", "allele"], sort=False)
    idx = g3.groups
    agg = g3.agg({"UMI": "count", "readCount": "sum"}).sort_values("UMI", ascending=False).reset_index()
    items = [(k, t["intBC"].values, t["UMI"].values)
             for k, t in agg.groupby(["cellBC", "allele"], sort=False) if len(t) > 1]
    if verbose: print(f"[faithful EC] {len(items):,} multi-intBC groups, n_jobs={n_jobs}", flush=True)
    B = 2000
    chunks = [items[i:i + B] for i in range(0, len(items), B)]
    res = Parallel(n_jobs=n_jobs, backend="loky")(delayed(_batch)(c, prop, umi_count_thresh, dist_thresh) for c in chunks)
    corr = [x for chunk in res for grp in chunk for x in grp]
    if verbose: print(f"[faithful EC] {len(corr):,} corrections (replayed in original order)", flush=True)
    col = molecule_table["intBC"].to_numpy(dtype=object, copy=True)
    pos_of = molecule_table.index.get_indexer
    for key, new in corr:                       # sequential replay == original .loc writes, last write wins
        col[pos_of(idx[key])] = new
    molecule_table["intBC"] = col
    return molecule_table
