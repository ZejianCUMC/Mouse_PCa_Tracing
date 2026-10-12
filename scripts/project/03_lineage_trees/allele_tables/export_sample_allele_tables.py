# =============================================================================
# Per-sample allele table export (Npp53 + Np53Rb1) Pipeline
# =============================================================================
# Description : Write one allele table per mouse, keeping only rows whose cell is a
#               valid final scRNA-seq tumour cell of that mouse and dropping un-detected
#               rows (empty allele or any empty r1/r2/r3 site). Uncut sites ("[None]")
#               are kept. Name: {Mouse_model}_{treatment}_{sample_id}_allele_table.csv
# Input       : Npp53: per_sample JZ*_allele_table.csv + data/npp53_merge/obs_metadata.csv
#               Np53Rb1: frozen tree inputs (cassiopeia_*_seed42/inputs/<S>/) + final metadata
# Output      : per-project folder of CSVs + SUMMARY.csv (+ Dropbox copies)
# Conda env   : /home/jrchristin/zejian/miniforge3/envs/scvi/bin/python
# Validated   : 2026-09-25
# Parameters  : See DEFAULTS dict below
# =============================================================================

# ── IMPORTS ──────────────────────────────────────────────────────────────────
import json, re, shutil, hashlib
from pathlib import Path
import pandas as pd

# ── DEFAULTS (all tunable parameters) ────────────────────────────────────────
MB = Path("/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based")
DB = Path("/mnt/7efffd2d-39a7-4d6a-87c7-773685593413/dropbox/Dropbox/share/Tracing")
DEFAULTS = {
    "npp53_per_sample": MB / "Npp53_tracing/lineage_tracing/allele_tables/per_sample",
    "npp53_obs": MB / "Npp53_tracing/data/npp53_merge/obs_metadata.csv",            # 20,978 valid tumour cells
    "npp53_out": MB / "Npp53_tracing/lineage_tracing/allele_tables/per_sample_final_20260925",
    "npp53_dropbox": DB / "Npp53_tracing/allele_tables/per_sample_final_20260925",
    "rb1_cas_runs": [MB / "Npp53Rb1_Tracing/lineage_tracing/cassiopeia_eight_20260923_seed42",
                     MB / "Npp53Rb1_Tracing/lineage_tracing/cassiopeia_mjz026_028_20260925_seed42"],
    "rb1_meta": MB / "Npp53Rb1_Tracing/analysis/tumor_states_20260922/exclude_cl5_cl8_seed42/"
                     "exclude_current_c17_c19_seed42/metadata_PRIVATE.csv",              # 40,141 valid tumour cells
    "rb1_out": MB / "Npp53Rb1_Tracing/lineage_tracing/allele_tables_per_sample_20260925",
    "rb1_dropbox": DB / "Npp53Rb1_tracing/allele_tables/per_sample_20260925",
    "sites": ["allele", "r1", "r2", "r3"],      # empty string in any of these = un-detected -> drop row
    "cols": ["cell_id", "cellBC", "intBC", "allele", "r1", "r2", "r3", "lineageGrp", "UMI", "readCount",
             "sample_id", "run_id", "treatment", "model", "putative_state"],
}
TREAT = {"Intact": "Intact", "Castration": "Castration", "Castration+ENZA": "CastrationENZA",
         "Castration + Enzalutamide": "CastrationENZA"}

# ── HELPER FUNCTIONS ─────────────────────────────────────────────────────────
def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(8 << 20), b""): h.update(b)
    return h.hexdigest()

def finish(df, sample, valid_cells, p):
    """Keep valid tumour cells of this sample; drop un-detected rows; order columns/rows."""
    n0, c0 = len(df), df.cell_id.nunique()
    df = df[df.cell_id.isin(valid_cells)]
    n1, c1 = len(df), df.cell_id.nunique()
    und = (df[p["sites"]].astype(str) == "").any(axis=1)
    df = df[~und]
    assert df[["intBC", "allele", "r1", "r2", "r3"]].ne("").all().all()
    df = df[p["cols"]].sort_values(["cell_id", "intBC", "allele"], kind="stable").reset_index(drop=True)
    stats = dict(sample_id=sample, rows_input=n0, cells_input=c0, rows_valid_cells=n1, cells_valid=c1,
                 undetected_rows_dropped=int(und.sum()), rows_final=len(df), cells_final=df.cell_id.nunique(),
                 intBCs=df.intBC.nunique(), valid_rna_tumour_cells=len(valid_cells),
                 pct_rna_cells_with_alleles=round(100 * df.cell_id.nunique() / len(valid_cells), 2))
    return df, stats

def write(df, stats, out, dbx):
    model, tr = df.model.iloc[0], TREAT[df.treatment.iloc[0]]
    name = f"{model}_{tr}_{stats['sample_id']}_allele_table.csv"
    f = out / name
    assert not f.exists(), f"refusing to overwrite {f}"
    df.to_csv(f, index=False)
    stats.update(file=name, model=model, treatment=df.treatment.iloc[0], run_id=df.run_id.iloc[0], md5=md5(f))
    return stats

# ── MAIN PIPELINE FUNCTION ───────────────────────────────────────────────────
def run_pipeline(input_path: str = None, output_dir: str = None, **kwargs) -> dict:
    p = {**DEFAULTS, **kwargs}
    summary = {}
    # ---- Npp53 (25 mice) ----
    out = p["npp53_out"]; out.mkdir(parents=True, exist_ok=False)
    obs = pd.read_csv(p["npp53_obs"], index_col=0)
    rows = []
    for f in sorted(p["npp53_per_sample"].glob("JZ*_allele_table.csv")):
        s = f.name.split("_")[0]
        a = pd.read_csv(f, keep_default_na=False, dtype=str)
        o = obs[obs.sample_id == s]
        a = a.rename(columns={"cell_full": "cell_id", "run": "run_id"})
        a = a.drop(columns=["treatment", "putative_state"]).merge(
            o[["treatment", "model", "putative_state"]], left_on="cell_id", right_index=True, how="left")
        df, st = finish(a, s, set(o.index), p)
        rows.append(write(df, st, out, p["npp53_dropbox"])); print("Npp53", rows[-1]["file"], st["cells_final"], flush=True)
    summary["Npp53"] = pd.DataFrame(rows); summary["Npp53"].to_csv(out / "SUMMARY.csv", index=False)
    # ---- Np53Rb1 (22 mice; full JZ135 table) ----
    out = p["rb1_out"]; out.mkdir(parents=True, exist_ok=False)
    meta = pd.read_csv(p["rb1_meta"], index_col=0)
    rows = []
    for run in p["rb1_cas_runs"]:
        for s in json.loads((run / "PREPARED.json").read_text())["samples"]:
            a = pd.read_csv(run / "inputs" / s / "allele_table.csv.gz", keep_default_na=False, dtype=str)
            m = meta[meta.sample_id == s]
            a = a.rename(columns={"cellBC": "cell_id"})
            a = a.merge(m[["source_barcode", "run_id", "treatment", "model"]], left_on="cell_id", right_index=True, how="left")
            a["cellBC"] = a.source_barcode.str.replace(r"-\d+$", "", regex=True)
            a["sample_id"] = s
            ls = pd.read_csv(run / "inputs" / s / "cell_metadata.csv", index_col=0)
            a["putative_state"] = a.cell_id.map(ls.putative_state)
            df, st = finish(a, s, set(m.index), p)
            rows.append(write(df, st, out, p["rb1_dropbox"])); print("Np53Rb1", rows[-1]["file"], st["cells_final"], flush=True)
    summary["Np53Rb1"] = pd.DataFrame(rows); summary["Np53Rb1"].to_csv(out / "SUMMARY.csv", index=False)
    # ---- Dropbox copies (deliverable summaries of sample-level tables) ----
    for key, o, d in [("Npp53", p["npp53_out"], p["npp53_dropbox"]), ("Np53Rb1", p["rb1_out"], p["rb1_dropbox"])]:
        d.mkdir(parents=True, exist_ok=False)
        for f in o.glob("*.csv"): shutil.copy2(f, d / f.name)
        bad = [f.name for f in o.glob("*_allele_table.csv") if md5(f) != md5(d / f.name)]
        assert not bad, bad
    return {"output_path": [str(p["npp53_out"]), str(p["rb1_out"])], "params_used": {k: str(v) for k, v in p.items()},
            "summary": {k: dict(files=len(v), cells=int(v.cells_final.sum()), rows=int(v.rows_final.sum())) for k, v in summary.items()}}

# ── CLI ENTRY POINT ──────────────────────────────────────────────────────────
if __name__ == "__main__":
    print(json.dumps(run_pipeline()["summary"], indent=2))
