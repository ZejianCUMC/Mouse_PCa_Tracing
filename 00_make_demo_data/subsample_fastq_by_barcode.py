#!/usr/bin/env python
# =============================================================================
# Demo-data builder: subsample raw FASTQs of one 10x run down to a few mice
# =============================================================================
# Description : Picks a small set of 10x cell barcodes from run MJZ019 and keeps
#               ONLY the read pairs carrying those barcodes, in all three
#               libraries of the run (GEX, hashtag/HTO, lineage amplicon TSL).
#               The barcode set =
#                 (a) up to N_CELLS HTO-singlet cells of each chosen mouse,
#                 (b) N_OTHER random non-singlet cell barcodes (doublets /
#                     negatives) so demultiplexing has something to reject,
#                 (c) N_EMPTY low-UMI "empty" droplets so CellBender can learn
#                     the ambient-RNA profile.
#               Barcode match = exact match of R1[0:16]. The result is a real,
#               small dataset that runs through the whole tutorial in ~hours.
# Input       : original FASTQs + Cell Ranger outs + HTODemux singlet table
# Output      : OUTDIR/fastq/{GEX,HTO,TSL}/*.fastq.gz + OUTDIR/manifest.json
# Conda env   : any python3 with numpy + h5py (HPC env `scvi`); needs `pigz`
# Validated   : 2026-09-27
# =============================================================================

# ── IMPORTS ──────────────────────────────────────────────────────────────────
import argparse, json, os, subprocess, sys, time
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd
import h5py

# ── DEFAULTS (all tunable parameters) ────────────────────────────────────────
DEFAULTS = {
    "samples": "JZ202,JZ204",   # mice to keep (Cas+Enz JZ202 / Castration JZ204)
    "n_cells": 800,             # max HTO-singlet cells kept per mouse
    "n_other": 300,             # extra non-singlet cell barcodes (doublets/negatives)
    "n_empty": 15000,           # empty droplets for CellBender
    "empty_min_umi": 5,         # empty droplet = GEX UMI in [empty_min_umi, empty_max_umi]
    "empty_max_umi": 200,
    "seed": 42,                 # lab default seed
    "pigz_threads": 2,          # compression threads per output file
}

# ── HELPER FUNCTIONS ─────────────────────────────────────────────────────────
def log(m): print(f"[{time.strftime('%F %T')}] {m}", flush=True)

def gex_umi_per_barcode(raw_h5):
    """Total Gene-Expression UMIs per barcode from a Cell Ranger raw_feature_bc_matrix.h5."""
    with h5py.File(raw_h5, "r") as f:
        g = f["matrix"]
        bcs = np.array([b.decode() for b in g["barcodes"][:]])
        ftype = np.array([t.decode() for t in g["features/feature_type"][:]])
        data, indices, indptr = g["data"][:], g["indices"][:], g["indptr"][:]
    is_gex = (ftype == "Gene Expression")
    col = np.repeat(np.arange(len(bcs)), np.diff(indptr))
    umi = np.bincount(col, weights=data * is_gex[indices], minlength=len(bcs))
    return pd.Series(umi, index=bcs)

def choose_barcodes(p):
    rng = np.random.default_rng(p["seed"])
    meta = pd.read_csv(p["hash_meta"])                     # HTO singlets: cell_barcode, sample_id, ...
    filt = pd.read_csv(os.path.join(p["cr_outs"], "filtered_feature_bc_matrix", "barcodes.tsv.gz"),
                       header=None)[0].tolist()
    keep, per_sample = [], {}
    for s in p["samples"].split(","):
        cells = meta.loc[meta["sample_id"] == s, "cell_barcode"].tolist()
        pick = list(rng.choice(cells, size=min(p["n_cells"], len(cells)), replace=False))
        keep += pick; per_sample[s] = len(pick)
    singlets = set(meta["cell_barcode"])
    non_singlet = sorted(set(filt) - singlets)
    other = list(rng.choice(non_singlet, size=min(p["n_other"], len(non_singlet)), replace=False))
    umi = gex_umi_per_barcode(os.path.join(p["cr_outs"], "raw_feature_bc_matrix.h5"))
    pool = umi[(umi >= p["empty_min_umi"]) & (umi <= p["empty_max_umi"])].index
    pool = sorted(set(pool) - set(filt))
    empty = list(rng.choice(pool, size=min(p["n_empty"], len(pool)), replace=False))
    all_bc = {b.split("-")[0] for b in keep + other + empty}          # 16-bp raw barcode
    summary = {"per_sample_cells": per_sample, "other_cells": len(other),
               "empty_droplets": len(empty), "total_barcodes": len(all_bc)}
    return all_bc, summary

def filter_pair(job):
    """Stream one R1/R2 pair; keep records whose R1[0:16] is in the barcode set."""
    r1, r2, o1, o2, bcfile, threads = job
    bcs = set(open(bcfile).read().split())
    p1 = subprocess.Popen(["pigz", "-dc", r1], stdout=subprocess.PIPE, bufsize=1 << 20)
    p2 = subprocess.Popen(["pigz", "-dc", r2], stdout=subprocess.PIPE, bufsize=1 << 20)
    w1 = subprocess.Popen(["pigz", "-p", str(threads), "-c"], stdin=subprocess.PIPE, stdout=open(o1, "wb"))
    w2 = subprocess.Popen(["pigz", "-p", str(threads), "-c"], stdin=subprocess.PIPE, stdout=open(o2, "wb"))
    f1, f2 = p1.stdout, p2.stdout
    n = k = 0
    while True:
        a = [f1.readline() for _ in range(4)]
        b = [f2.readline() for _ in range(4)]
        if not a[0]:
            break
        n += 1
        if a[1][:16].decode() in bcs:
            k += 1
            w1.stdin.write(b"".join(a)); w2.stdin.write(b"".join(b))
    for w in (w1, w2):
        w.stdin.close(); w.wait()
    p1.wait(); p2.wait()
    return os.path.basename(o1), n, k

# ── MAIN PIPELINE FUNCTION ───────────────────────────────────────────────────
def run_pipeline(input_path: str, output_dir: str, **kwargs) -> dict:
    """input_path = JSON listing the source FASTQs / Cell Ranger outs / singlet table."""
    params = {**DEFAULTS, **json.load(open(input_path)), **kwargs}
    os.makedirs(output_dir, exist_ok=True)
    bcs, summary = choose_barcodes(params)
    bcfile = os.path.join(output_dir, "selected_barcodes.txt")
    with open(bcfile, "w") as f:
        f.write("\n".join(sorted(bcs)) + "\n")
    log(f"barcodes selected: {summary}")
    jobs = []
    for lib, (r1, r2, o1, o2) in params["fastqs"].items():
        d = os.path.join(output_dir, "fastq", lib); os.makedirs(d, exist_ok=True)
        for o in (o1, o2):
            if os.path.exists(os.path.join(d, o)):
                sys.exit(f"refusing to overwrite existing {os.path.join(d, o)}")
        jobs.append((r1, r2, os.path.join(d, o1), os.path.join(d, o2), bcfile, params["pigz_threads"]))
    with ProcessPoolExecutor(max_workers=len(jobs)) as ex:
        for name, n, k in ex.map(filter_pair, jobs):
            log(f"{name}: kept {k:,} / {n:,} read pairs ({100 * k / max(n, 1):.2f}%)")
            summary[name] = {"reads_in": n, "reads_kept": k}
    params_used = {k: v for k, v in params.items()}
    json.dump({"summary": summary, "params_used": params_used},
              open(os.path.join(output_dir, "manifest.json"), "w"), indent=2)
    return {"output_path": output_dir, "params_used": params_used, "summary": summary}

# ── CLI ENTRY POINT ──────────────────────────────────────────────────────────
if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sources", required=True, help="JSON with hash_meta, cr_outs, fastqs")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    run_pipeline(a.sources, a.outdir)
