"""
build_qc.py — raw-sequencing QC for the new MJZ019/021/022 TSL batch, with
cross-batch comparison against MJZ005, MJZ008, MJZ013-15 (all Element AVITI
bases2fastq 'TSL' libraries).

All numbers read from the bases2fastq stats on HPC (IndexAssignment.csv,
Metrics.csv, RunParameters.json, *_stats.json) surveyed 2026-07-15 —
no FASTQ decompression needed. Self-contained HTML in the house style.

Run: /home/jrchristin/zejian/miniforge3/envs/scvi/bin/python build_qc.py
"""
from __future__ import annotations
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, "/media/jrchristin/Backup10TB1/agent/1-project/23-ancestorOT/report/lib")
import note_kit as nk  # noqa: E402
from note_kit import (C_INK, C_MUTE, C_GREEN, C_GOLD, C_RED, C_BLUE,  # noqa: E402
                      C_ORANGE, C_PURPLE, to_b64)
import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

# ── extracted data (bases2fastq stats, 2026-07-15) ──────────────────────────
BATCH = ["MJZ005", "MJZ008", "MJZ013-15", "MJZ019-22"]
BCOL = {"MJZ005": C_MUTE, "MJZ008": C_BLUE, "MJZ013-15": C_GOLD, "MJZ019-22": C_RED}
RUN = {
    "MJZ005": ("MJZ005TSL", "2025-11-16", "28 / 280"),
    "MJZ008": ("TSLJZ008", "2026-01-27", "28 / 280"),
    "MJZ013-15": ("MJZ0131415TSL", "2026-05-05", "29 / 289"),
    "MJZ019-22": ("MJZ019-21-22TSL", "2026-07-13", "29 / 289"),
}
# per-sample: (sample, reads_pairs, yield_gb, pct_q30, R1_GB, R2_GB, total_GB)
# pct_q30 = authoritative per-sample value from <sample>_stats.json
SAMPLES = {
    "MJZ005": [("MJZ005TSL", 139_430_305, 41.61, 98.35, 3.4, 11.0, 14.4)],
    "MJZ008": [("JZ008", 158_893_074, 46.89, 94.34, 4.7, 20.0, 24.7)],
    "MJZ013-15": [("JZ013", 91_468_884, 27.79, 98.05, 2.2, 7.0, 9.2),
                  ("JZ014", 151_450_429, 45.13, 98.04, 3.5, 12.0, 15.5),
                  ("JZ015", 89_609_324, 27.84, 97.73, 2.1, 7.0, 9.1)],
    "MJZ019-22": [("MJZ019", 115_690_399, 30.93, 97.44, 2.7, 8.9, 12.0),
                  ("MJZ021", 78_913_261, 24.40, 96.38, 1.9, 7.3, 9.2),
                  ("MJZ022", 75_073_734, 23.18, 96.42, 1.8, 7.1, 8.8)],
}
# per-sample extras: (pct_q40, mean_q, pct_mismatch, pct_of_run_assigned)
EXTRA = {
    "MJZ005TSL": (94.10, 43.54, 0.446, 28.4), "JZ008": (83.41, 41.41, 0.540, 23.2),
    "JZ013": (93.90, 43.37, 0.551, 24.9), "JZ014": (93.87, 43.38, 0.478, 41.2),
    "JZ015": (92.96, 43.18, 0.918, 24.4), "MJZ019": (92.66, 43.03, 0.401, 42.3),
    "MJZ021": (90.03, 42.47, 0.513, 28.9), "MJZ022": (89.97, 42.47, 1.033, 27.5),
}
# run-level: (pct_assigned, assigned_yield_gb, total_yield_gb, q30)
RUNLVL = {
    "MJZ005": (28.4, 41.61, 146.44, 97.76),
    "MJZ008": (23.2, 46.89, 202.09, 94.91),
    "MJZ013-15": (90.5, 100.77, 111.37, 97.86),
    "MJZ019-22": (98.7, 78.52, 79.57, 96.77),
}


def _samps_flat():
    out = []
    for b in BATCH:
        for s in SAMPLES[b]:
            out.append((b,) + s)
    return out


# =========================================================================
# Fig 1 — reads per sample, grouped by batch
# =========================================================================
def fig_reads() -> str:
    fig, ax = plt.subplots(figsize=(12.5, 5.2))
    flat = _samps_flat()
    x = np.arange(len(flat))
    vals = [r[2] / 1e6 for r in flat]
    cols = [BCOL[r[0]] for r in flat]
    labs = [r[1] for r in flat]
    ax.bar(x, vals, color=cols, edgecolor="white")
    for i, v in enumerate(vals):
        ax.text(i, v + 2, f"{v:.0f}M", ha="center", fontsize=8.8, fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(labs, rotation=25, ha="right", fontsize=9)
    ax.set_ylabel("read pairs (millions)", fontsize=10.5)
    ax.set_ylim(0, 175)
    ax.axhspan(0, 50, color=C_RED, alpha=0.05)
    ax.text(0.2, 46, "≈50M — practical floor for deep intBC recovery", fontsize=8, color=C_RED, style="italic")
    ax.set_title("Reads per TSL sample (Element AVITI 'polonies' = read pairs)",
                 fontsize=12.5, loc="left", fontweight="bold")
    import matplotlib.patches as _mp
    handles = [_mp.Patch(facecolor=BCOL[b], label=b) for b in BATCH]
    ax.legend(handles=handles, title="batch", fontsize=8.6, frameon=False,
              loc="upper right", ncol=2)
    plt.tight_layout()
    return to_b64(fig)


# =========================================================================
# Fig 2 — % assigned trend (the headline)
# =========================================================================
def fig_assigned() -> str:
    fig, ax = plt.subplots(figsize=(11, 5.2))
    x = np.arange(len(BATCH))
    pa = [RUNLVL[b][0] for b in BATCH]
    n_s = [len(SAMPLES[b]) for b in BATCH]
    ax.plot(x, pa, "-o", color=C_RED, lw=2.4, markersize=10, zorder=3)
    for i, (v, n) in enumerate(zip(pa, n_s)):
        ax.text(i, v + 4, f"{v:.0f}%", ha="center", fontsize=11, fontweight="bold", color=C_RED)
        ax.text(i, 8, f"{n} sample{'s' if n>1 else ''}\nmultiplexed", ha="center",
                fontsize=8.4, color=C_MUTE)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{b}\n{RUN[b][1]}" for b in BATCH], fontsize=9)
    ax.set_ylabel("% of flowcell reads assigned to target index", fontsize=10.5)
    ax.set_ylim(0, 108)
    ax.set_title("Index-assignment efficiency improved sharply over the four runs",
                 fontsize=12.5, loc="left", fontweight="bold")
    ax.text(0.02, 0.44, "single-sample runs on shared flowcells\n→ most polonies were other indices",
            transform=ax.transAxes, fontsize=8.4, color=C_MUTE, style="italic")
    ax.annotate("dedicated multiplexed runs\n→ nearly all reads are yours",
                xy=(3, 98.7), xytext=(2.2, 66), fontsize=8.6, color=C_GREEN, style="italic",
                arrowprops=dict(arrowstyle="->", color=C_GREEN, lw=1.3))
    plt.tight_layout()
    return to_b64(fig)


# =========================================================================
# Fig 3 — total vs assigned yield per batch
# =========================================================================
def fig_yield() -> str:
    fig, ax = plt.subplots(figsize=(11, 5.0))
    x = np.arange(len(BATCH)); w = 0.38
    tot = [RUNLVL[b][2] for b in BATCH]
    asg = [RUNLVL[b][1] for b in BATCH]
    ax.bar(x - w/2, tot, w, label="total flowcell yield", color="#cfc8b8", edgecolor="white")
    ax.bar(x + w/2, asg, w, label="assigned to these samples", color=C_GREEN, edgecolor="white")
    for i in range(len(BATCH)):
        ax.text(i - w/2, tot[i] + 3, f"{tot[i]:.0f}", ha="center", fontsize=8.4, color=C_MUTE)
        ax.text(i + w/2, asg[i] + 3, f"{asg[i]:.0f}", ha="center", fontsize=8.4, fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(BATCH, fontsize=9.5)
    ax.set_ylabel("yield (Gb)", fontsize=10.5); ax.set_ylim(0, 225)
    ax.set_title("Total vs assigned yield — the old runs bought yield they didn't keep",
                 fontsize=12, loc="left", fontweight="bold")
    ax.legend(fontsize=9, frameon=False)
    plt.tight_layout()
    return to_b64(fig)


# =========================================================================
# Fig 4 — per-sample file size (R1 + R2 stacked)
# =========================================================================
def fig_size() -> str:
    fig, ax = plt.subplots(figsize=(12.5, 4.8))
    flat = _samps_flat()
    x = np.arange(len(flat))
    r1 = [r[5] for r in flat]; r2 = [r[6] for r in flat]
    ax.bar(x, r1, color="#9ecae1", edgecolor="white", label="R1 (barcode+UMI, 28 bp)")
    ax.bar(x, r2, bottom=r1, color=C_ORANGE, edgecolor="white", label="R2 (insert, 280–289 bp)")
    for i, r in enumerate(flat):
        ax.text(i, r[7] + 0.4, f"{r[7]:.0f}G", ha="center", fontsize=8.4, fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels([r[1] for r in flat], rotation=25, ha="right", fontsize=9)
    ax.set_ylabel("gzipped FASTQ size (GB)", fontsize=10.5); ax.set_ylim(0, 28)
    ax.set_title("Per-sample FASTQ footprint — R2 (long insert) dominates",
                 fontsize=12, loc="left", fontweight="bold")
    ax.legend(fontsize=9, frameon=False, loc="upper right")
    plt.tight_layout()
    return to_b64(fig)


# =========================================================================
# tables
# =========================================================================
def tbl_new() -> str:
    body = ""
    for s, reads, y, q30, r1, r2, tot in SAMPLES["MJZ019-22"]:
        body += (f"<tr><td><b>{s}</b></td><td>{reads:,}</td><td>{y:.2f} Gb</td>"
                 f"<td>{q30:.1f}%</td><td>{r1:.1f}G + {r2:.1f}G = <b>{tot:.0f}G</b></td></tr>")
    body += ("<tr class='tot'><td><b>Σ run</b></td><td><b>269,677,394</b></td>"
             "<td><b>78.52 Gb</b> (98.7% assigned)</td><td>96.8%</td><td><b>≈30G</b> fastq</td></tr>")
    return ("<table><thead><tr><th>Sample</th><th>Read pairs</th><th>Yield</th>"
            f"<th>%Q30</th><th>FASTQ size (R1 / R2)</th></tr></thead><tbody>{body}</tbody></table>")


def tbl_combined() -> str:
    body = ""
    for b in BATCH:
        cyc = RUN[b][2]; date = RUN[b][1]
        for i, (s, reads, y, q30, r1, r2, tot) in enumerate(SAMPLES[b]):
            q40, mq, mm, pr = EXTRA.get(s, (0, 0, 0, 0))
            ref = " ref" if b in ("MJZ005", "MJZ008") else ""
            bcell = (f"<td rowspan='{len(SAMPLES[b])}'><b>{b}</b><br>"
                     f"<span class='mono2'>{cyc} · {date}</span></td>") if i == 0 else ""
            body += (f"<tr class='{('refrow' if ref else '')}'>{bcell}"
                     f"<td><b>{s}</b></td><td>{reads:,}</td><td>{y:.2f}</td>"
                     f"<td>{q30:.2f}</td><td>{q40:.2f}</td><td>{mq:.2f}</td>"
                     f"<td>{mm:.2f}</td><td>{pr:.0f}%</td>"
                     f"<td>{tot:.0f}G <span class='mono2'>({r1:.1f}+{r2:.1f})</span></td></tr>")
    return ("<table class='comb'><thead><tr><th>Batch (R1/R2 · date)</th><th>Sample</th>"
            "<th>Read pairs</th><th>Yield Gb</th><th>%Q30</th><th>%Q40</th><th>mean Q</th>"
            "<th>%mis</th><th>%run</th><th>FASTQ (R1+R2)</th></tr></thead>"
            f"<tbody>{body}</tbody></table>")


def tbl_compare() -> str:
    body = ""
    for b in BATCH:
        run, date, cyc = RUN[b]
        pa, ay, ty, q30 = RUNLVL[b]
        ns = len(SAMPLES[b])
        rd = sum(s[1] for s in SAMPLES[b]) / 1e6
        body += (f"<tr><td><b>{b}</b></td><td class='mono'>{run}</td><td>{date}</td>"
                 f"<td>{cyc}</td><td>{ns}</td><td>{rd:.0f}M</td>"
                 f"<td>{q30:.1f}%</td><td>{pa:.0f}%</td><td>{ay:.0f}/{ty:.0f} Gb</td></tr>")
    return ("<table><thead><tr><th>Batch</th><th>Run</th><th>Date</th>"
            "<th>R1/R2 cyc</th><th>n</th><th>Σ reads</th><th>Q30</th>"
            f"<th>%assigned</th><th>assigned/total</th></tr></thead><tbody>{body}</tbody></table>")


BODY = r"""
<div class="crumb">12-Jingbo / 2-Tracing / 0-raw &nbsp;·&nbsp; raw-sequencing QC</div>

<h1>MJZ019 / 021 / 022 — raw-sequencing QC &amp; cross-batch comparison</h1>
<div class="subtitle">
New Element-AVITI <b>bases2fastq</b> batch <code>MJZ019-21-22TSL</code>
(2026-07-13), compared against the three prior tracing-library runs
(<b>MJZ005</b>, <b>MJZ008</b>, <b>MJZ013–15</b>). All read counts, yields and
quality come straight from the bases2fastq stats files — no FASTQ decompression,
so <b>no compute job was needed</b>.
</div>

<div class="callout tip">
<span class="callout-label">Headline</span>
The new batch is <b>clean and deep</b>: <b>269.7M read pairs</b> across 3 samples
(MJZ019 116M, MJZ021 79M, MJZ022 75M), <b>98.7% assigned</b> to the correct index,
<b>Q30 96.8%</b>, ~30 GB of FASTQ. Same "TSL" library design as every prior run
(10x barcode R1 = 28 bp, long insert R2 ≈ 289 bp, dual 8 bp index). Index-assignment
efficiency has climbed run-over-run from <b>23% → 99%</b>.
</div>

<h2>1 &middot; The new batch (MJZ019 / 021 / 022)</h2>
__NEW__
<p class="small">Run <code>MJZ019-21-22TSL</code> · instrument AV100007 · flowcell
2551406424 · 2 lanes · 300-cycle kit · dual 8 bp index. Reads are balanced across
the three samples (42 / 29 / 27% of the run). <b>MJZ020 is absent</b> — only 019,
021, 022 were sequenced.</p>
<figure><img src="data:image/png;base64,{fig_reads}" alt="reads per sample"></figure>

<h2>2 &middot; Read structure &amp; chemistry (identical across all four batches)</h2>
<table><thead><tr><th>Feature</th><th>Value</th><th>Meaning</th></tr></thead><tbody>
<tr><td>Platform</td><td>Element AVITI (AV100007), Cloudbreak-FS, 300-cycle kit</td><td>bases2fastq demux</td></tr>
<tr><td>R1</td><td>28–29 cycles, mask <code>Y28N*</code></td><td>16 bp cell barcode + 12 bp UMI (10x)</td></tr>
<tr><td>R2</td><td>280–289 cycles (adapter-trimmed)</td><td><b>long insert</b> — reads the lineage/target-seq recorder amplicon</td></tr>
<tr><td>Index</td><td>dual 8 bp (i7 + i5)</td><td>sample demultiplex</td></tr>
<tr><td>Adapters</td><td>R1 <code>AGATCGG…ACTCCAGTCA</code>, R2 <code>AGATCGG…AGAGTGT</code></td><td>standard TruSeq; trimmed + N-masked</td></tr>
</tbody></table>
<div class="callout bio">
<span class="callout-label">Library type</span>
Every run is a <b>"TSL"</b> library with a short 10x barcode R1 and a very long
(~289 bp) R2 that is heavily adapter-trimmed — the signature of a
<b>single-cell lineage / target-seq (recorder-amplicon) library</b>, not a plain
short-R2 GEX library. The design is unchanged from MJZ005 onward, so the new batch
is directly comparable to the earlier tracing runs.
</div>

<h2>3 &middot; Cross-batch comparison</h2>
<p>Full per-sample table across all four batches (per-sample values from each
<code>&lt;sample&gt;_stats.json</code>). <b>MJZ005 &amp; MJZ008 (shaded) are the
established "good" reference samples</b> — the deep, successful tracing runs that
produced JZ135–138.</p>
__COMBINED__
<p class="small">MJZ005TSL is the highest-quality library (Q30 98.4%, Q40 94.1%);
MJZ008/JZ008 is the deepest (159M pairs) but has the weakest base quality of the
set (Q40 83.4% — worth noting for its target-seq base-calling). <code>%run</code>
= share of the flowcell's assigned reads that went to that sample.</p>
<h3>Run-level summary</h3>
__COMPARE__
<figure><img src="data:image/png;base64,{fig_assigned}" alt="assigned trend"></figure>
<p>The single most useful comparison: <b>index-assignment efficiency</b>. The two
earliest runs (MJZ005, MJZ008) were single-sample and recovered only <b>23–28%</b>
of flowcell polonies to the target index — the rest were other indices (shared
flowcell / low index diversity), so the big <b>146–202 Gb total yields were mostly
not kept</b>. The multiplexed runs (MJZ013–15, MJZ019–22) assign <b>90–99%</b> —
far more efficient use of the flowcell.</p>
<figure><img src="data:image/png;base64,{fig_yield}" alt="yield"></figure>
<figure><img src="data:image/png;base64,{fig_size}" alt="file size"></figure>
<p>Per sample, R2 (the long insert) dominates the file size. The new batch's
7–9 GB R2 files are in line with MJZ013–15; MJZ008's 20 GB R2 reflects its very
deep single-sample run.</p>

<h2>4 &middot; Disk footprint &amp; a flag</h2>
<table><thead><tr><th>Artifact</th><th>Size</th><th>Note</th></tr></thead><tbody>
<tr><td><code>MJZ019_21_22_bases2fastq/</code> (FASTQ + stats)</td><td>≈30 GB</td><td>the deliverable</td></tr>
<tr><td><code>20260713_AV100007_MJZ019-21-22TSL.tar</code></td><td>50 GB</td><td>raw run-folder archive</td></tr>
<tr><td><code>20260713_AV100007_MJZ019-21-22TSL/</code> (extracted)</td><td>50 GB</td><td>extracted run folder</td></tr>
</tbody></table>
<div class="callout warn">
<span class="callout-label">Flag (no action taken)</span>
The raw run folder exists <b>both</b> as a 50 GB tar <b>and</b> as a 50 GB extracted
directory — ~50 GB is redundant. If the tar is the keep-copy, the extracted folder
(or vice-versa) could be removed to reclaim ~50 GB. <b>I did not delete anything</b>
— HPC deletions need your go-ahead.
</div>

<h2>5 &middot; Open items</h2>
<ul>
<li><b>Biological identity of MJZ019/021/022</b> — the sequencing manifests carry
only index sequences. What each sample is (genotype / model, as in the
<code>MJZ16to18.csv</code> sheet) needs your wet-lab sample sheet. I could not find
a <code>MJZ19–22</code> description on HPC.</li>
<li><b>GEX vs lineage</b> — these are the long-R2 "TSL" libraries. If a matched
short-R2 <b>GEX</b> run exists for 019/021/022 (as for 013–18), point me to it and
I'll characterise it too.</li>
<li><b>Deeper FASTQ-level QC</b> (per-base quality re-derivation, adapter content,
gzip integrity, exact read-length distribution) would require reading the actual
reads — <b>that</b> is an sbatch job (FastQC/seqkit), which I can submit on request.
The counts/features above did not need it.</li>
</ul>

<div class="footer">
Compiled 2026-07-15 from bases2fastq stats on
<code>hpc.c2b2.columbia.edu:/groups/ms3625_gp/zw2994/…/0-raw/</code>
(IndexAssignment.csv · Metrics.csv · RunParameters.json · *_stats.json).
Per-sample table: <code>raw_seq_qc/raw_seq_qc_metrics.csv</code>. Rebuild:
<code>raw_seq_qc/build_qc.py</code>.
</div>
"""

EXTRA_CSS = """
.small { font-size: 0.86rem; color: var(--mute); }
.mono { font-family: "JetBrains Mono","SF Mono",Menlo,monospace; font-size: 0.82rem; }
tr.tot td { border-top: 2px solid var(--line); background: #f9f7f0; }
.callout.warn ul, .callout.bio ul { margin: 0.4em 0 0; padding-left: 1.2em; }
table.comb { font-size: 0.85rem; }
table.comb td, table.comb th { padding: 5px 8px; }
.mono2 { font-family: "JetBrains Mono","SF Mono",Menlo,monospace; font-size: 0.76rem; color: var(--mute); }
tr.refrow td { background: #fbf3e8; }
"""


def main():
    figs = {"fig_reads": fig_reads(), "fig_assigned": fig_assigned(),
            "fig_yield": fig_yield(), "fig_size": fig_size()}
    body = (BODY.replace("__NEW__", tbl_new())
                .replace("__COMBINED__", tbl_combined())
                .replace("__COMPARE__", tbl_compare()))
    head = nk.HEAD.replace("</style>", EXTRA_CSS + "\n</style>")
    html = head.replace("{title}", "MJZ019-22 raw-seq QC") + body + nk.FOOT
    for k, v in figs.items():
        html = html.replace("{" + k + "}", v)
    out = HERE / "MJZ019_22_raw_seq_QC.html"
    out.write_text(html, encoding="utf-8")
    print(f"Wrote {out}  ({out.stat().st_size/1024:.0f} KB)  imgs={len(figs)}")


if __name__ == "__main__":
    main()
