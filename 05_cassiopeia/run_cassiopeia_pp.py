#!/usr/bin/env python
# =============================================================================
# 05 · Cassiopeia preprocessing: lineage amplicon FASTQ -> allele table
# =============================================================================
# Description : Target-site (TSL) reads -> per-cell, per-intBC indel alleles at
#               the 3 Cas9 cut sites (r1/r2/r3). Recipe used for every Npp53
#               library (MJZ005/008/013-15/019-025). Resumable: each step
#               writes a checkpoint CSV and is skipped on restart.
#   R1 = 16 bp cell barcode + 12 bp UMI (10x 3' v3);  R2 = amplicon
#   1 FASTQ -> unmapped BAM      2 filter_bam (Q>=10)   3 cellBC -> 10x whitelist
#   4 collapse_umis              5 resolve_umi_sequence 6 align to refOct1.fa
#   7 call_alleles (cut sites 112/166/220; intBC = ref interval 20-34)
#   8 error_correct_umis         9 filter_molecule_table
#  10 call_lineage_groups -> <LIB>_allele_table.csv
# Input       : TSL R1/R2 fastq.gz, refOct1.fa, 10x v3 whitelist
# Output      : OUTDIR/<LIB>_allele_table.csv
#               (cellBC, intBC, allele, r1, r2, r3, lineageGrp, UMI, readCount)
# Conda env   : HPC ${ENVS}/cassiopeia_env  (Cassiopeia 2.1)
# Validated   : 2026-07-17 on MJZ019/021/022 (full libraries)
# Usage       : run_cassiopeia_pp.py LIB R1 R2 OUTDIR NTHREADS REF WHITELIST
# =============================================================================
import os, sys, time
import pandas as pd
import cassiopeia as cas

# ── DEFAULTS (all tunable parameters) ────────────────────────────────────────
DEFAULTS = {
    "chemistry": "10xv3",
    "quality_threshold": 10,          # filter_bam: min mean base quality
    "max_hq_mismatches": 3,           # collapse_umis
    "max_indels": 2,                  # collapse_umis
    "min_umi_per_cell": 10,           # cells need >= 10 UMIs (resolve / filter / lineage groups)
    "min_avg_reads_per_umi": 2.0,
    "gap_open_penalty": 20, "gap_extend_penalty": 1,   # align_sequences (Smith-Waterman)
    "barcode_interval": (20, 34),     # intBC position in refOct1.fa
    "cutsite_locations": [112, 166, 220],               # r1 / r2 / r3 cut sites in refOct1.fa
    "cutsite_width": 12,
    "max_umi_distance": 2,            # error_correct_umis
    "intbc_prop_thresh": 0.5, "intbc_umi_thresh": 10, "intbc_dist_thresh": 1,
    "doublet_threshold": 0.35,        # filter_molecule_table allele-conflict doublets
    "min_cluster_prop": 0.005,        # call_lineage_groups
    "min_intbc_thresh": 0.05, "inter_doublet_threshold": 0.35, "kinship_thresh": 0.25,
}

def run_pipeline(LIB, R1, R2, OUTDIR, NTHREADS, REF, WHITELIST, **kw):
    p = {**DEFAULTS, **kw}
    os.makedirs(OUTDIR, exist_ok=True)
    ck = lambda n: os.path.join(OUTDIR, f"{LIB}_{n}.csv")
    log = lambda m: print(f"[{time.strftime('%F %T')}] {LIB}: {m}", flush=True)
    log(f"START n_threads={NTHREADS}")
    cas.pp.setup(OUTDIR, verbose=True)

    if os.path.exists(ck("collapsed")):
        log("resume: collapsed"); umi = pd.read_csv(ck("collapsed"))
    else:
        bam = cas.pp.convert_fastqs_to_unmapped_bam([R1, R2], chemistry=p["chemistry"],
                  output_directory=OUTDIR, name=LIB, n_threads=NTHREADS); log("1 convert done")
        bam = cas.pp.filter_bam(bam, output_directory=OUTDIR,
                  quality_threshold=p["quality_threshold"], n_threads=NTHREADS); log("2 filter done")
        bam = cas.pp.error_correct_cellbcs_to_whitelist(bam, whitelist=WHITELIST,
                  output_directory=OUTDIR, n_threads=NTHREADS); log("3 cellBC-correct done")
        umi = cas.pp.collapse_umis(bam, output_directory=OUTDIR, max_hq_mismatches=p["max_hq_mismatches"],
                  max_indels=p["max_indels"], method="likelihood", n_threads=NTHREADS)
        umi.to_csv(ck("collapsed"), index=False); log(f"4 collapse done {umi.shape}")

    if os.path.exists(ck("resolved")):
        log("resume: resolved"); umi = pd.read_csv(ck("resolved"))
    else:
        umi = cas.pp.resolve_umi_sequence(umi, output_directory=OUTDIR,
                  min_umi_per_cell=p["min_umi_per_cell"], min_avg_reads_per_umi=p["min_avg_reads_per_umi"], plot=False)
        umi.to_csv(ck("resolved"), index=False); log(f"5 resolve done {umi.shape}")

    if os.path.exists(ck("aligned")):
        log("resume: aligned"); umi = pd.read_csv(ck("aligned"))
    else:
        umi = cas.pp.align_sequences(umi, ref_filepath=REF, gap_open_penalty=p["gap_open_penalty"],
                  gap_extend_penalty=p["gap_extend_penalty"], n_threads=NTHREADS)
        umi.to_csv(ck("aligned"), index=False); log(f"6 align done {umi.shape}")

    if os.path.exists(ck("call_alleles")):
        log("resume: call_alleles"); umi = pd.read_csv(ck("call_alleles"))
    else:
        umi = cas.pp.call_alleles(umi, ref_filepath=REF, barcode_interval=p["barcode_interval"],
                  cutsite_locations=p["cutsite_locations"], cutsite_width=p["cutsite_width"],
                  context=True, context_size=5)
        umi.to_csv(ck("call_alleles"), index=False); log(f"7 call_alleles done {umi.shape}")

    if os.path.exists(ck("error_correct_umis")):
        log("resume: error_correct_umis"); umi = pd.read_csv(ck("error_correct_umis"))
    else:
        umi = cas.pp.error_correct_umis(umi, max_umi_distance=p["max_umi_distance"],
                  allow_allele_conflicts=False, n_threads=NTHREADS)
        umi.to_csv(ck("error_correct_umis"), index=False); log(f"8 error_correct_umis done {umi.shape}")

    if os.path.exists(ck("filtered_molecule_table")):
        log("resume: filtered_molecule_table"); umi = pd.read_csv(ck("filtered_molecule_table"))
    else:
        umi = cas.pp.filter_molecule_table(umi, output_directory=OUTDIR,
                  min_umi_per_cell=p["min_umi_per_cell"], min_avg_reads_per_umi=p["min_avg_reads_per_umi"],
                  min_reads_per_umi=-1, intbc_prop_thresh=p["intbc_prop_thresh"],
                  intbc_umi_thresh=p["intbc_umi_thresh"], intbc_dist_thresh=p["intbc_dist_thresh"],
                  doublet_threshold=p["doublet_threshold"], allow_allele_conflicts=False, plot=False)
        umi.to_csv(ck("filtered_molecule_table"), index=False); log(f"9 filter_molecule_table done {umi.shape}")

    allele = cas.pp.call_lineage_groups(umi, output_directory=OUTDIR,
                 min_umi_per_cell=p["min_umi_per_cell"], min_avg_reads_per_umi=p["min_avg_reads_per_umi"],
                 min_cluster_prop=p["min_cluster_prop"], min_intbc_thresh=p["min_intbc_thresh"],
                 inter_doublet_threshold=p["inter_doublet_threshold"], kinship_thresh=p["kinship_thresh"], plot=False)
    allele.to_csv(ck("allele_table"), index=False)
    summary = {"rows": len(allele), "cells": allele["cellBC"].nunique(), "intBCs": allele["intBC"].nunique()}
    log(f"10 DONE allele_table {summary}")
    return {"output_path": ck("allele_table"), "params_used": p, "summary": summary}

if __name__ == "__main__":
    if len(sys.argv) != 8:
        sys.exit(__doc__ if __doc__ else "usage: LIB R1 R2 OUTDIR NTHREADS REF WHITELIST")
    LIB, R1, R2, OUTDIR, NT, REF, WL = sys.argv[1:]
    run_pipeline(LIB, R1, R2, OUTDIR, int(NT), REF, WL)
