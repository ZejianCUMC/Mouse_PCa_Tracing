"""[STANDARD for all future TSL libraries — Zejian 2026-10-06] Cassiopeia preprocessing recipe of cassiopeia_MJZ19_22/scripts/run_cassiopeia_pp.py
with a FIXED filter_molecule_table min_reads_per_umi=1 (identical to the MJZ019 re-run). Do not use -1 (dynamic: gave 6 for MJZ019, 0-2 elsewhere).
Run uninterrupted; if a resume is unavoidable use the validated faithful wrapper (keep_default_na=False + ec_intbc_faithful).
Cassiopeia target-seq preprocessing for one 10x TSL library.
Replicates the established MJZ005/008/013-15 recipe (Cassiopeia_MJZ008.ipynb).
Resumable: each step writes a checkpoint CSV and is skipped if present, so a
restart never repeats the multi-hour cellBC-correction / collapse steps.

Usage: run_cassiopeia_pp.py LIB R1.fastq.gz R2.fastq.gz OUTDIR NTHREADS
Output: OUTDIR/<LIB>_allele_table.csv  (cellBC,intBC,r1,r2,r3,allele,lineageGrp,UMI,readCount,...)
"""
import os, sys, time
import pandas as pd
import cassiopeia as cas

LIB, R1, R2, OUTDIR, NTHREADS = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5])
REF = "/groups/ms3625_gp/zw2994/1-project/12-Jingbo/Nov13/refOct1.fa"
WHITELIST = "/groups/ms3625_gp/zw2994/1-project/12-Jingbo/Nov13/3M-3pgex-may-2023.txt"
os.makedirs(OUTDIR, exist_ok=True)
def ck(n): return os.path.join(OUTDIR, f"{LIB}_{n}.csv")
def log(m): print(f"[{time.strftime('%F %T')}] {LIB}: {m}", flush=True)

log(f"START n_threads={NTHREADS}  R1={R1}  R2={R2}")
cas.pp.setup(OUTDIR, verbose=True)

# 1-4: FASTQ -> unmapped BAM -> filter -> cellBC->whitelist -> collapse UMIs
if os.path.exists(ck("collapsed")):
    log("resume: load collapsed"); umi = pd.read_csv(ck("collapsed"))
else:
    bam = cas.pp.convert_fastqs_to_unmapped_bam([R1, R2], chemistry='10xv3',
              output_directory=OUTDIR, name=LIB, n_threads=NTHREADS); log("convert done")
    bam = cas.pp.filter_bam(bam, output_directory=OUTDIR, quality_threshold=10, n_threads=NTHREADS); log("filter done")
    bam = cas.pp.error_correct_cellbcs_to_whitelist(bam, whitelist=WHITELIST,
              output_directory=OUTDIR, n_threads=NTHREADS); log("cellBC-correct done")
    umi = cas.pp.collapse_umis(bam, output_directory=OUTDIR, max_hq_mismatches=3,
              max_indels=2, method='likelihood', n_threads=NTHREADS)
    umi.to_csv(ck("collapsed"), index=False); log(f"collapse done  {umi.shape}")

# 5: resolve UMI sequence
if os.path.exists(ck("resolved")):
    log("resume: load resolved"); umi = pd.read_csv(ck("resolved"))
else:
    umi = cas.pp.resolve_umi_sequence(umi, output_directory=OUTDIR,
              min_umi_per_cell=10, min_avg_reads_per_umi=2.0, plot=False)
    umi.to_csv(ck("resolved"), index=False); log(f"resolve done  {umi.shape}")

# 6: align to target-site reference
if os.path.exists(ck("aligned")):
    log("resume: load aligned"); umi = pd.read_csv(ck("aligned"))
else:
    umi = cas.pp.align_sequences(umi, ref_filepath=REF, gap_open_penalty=20,
              gap_extend_penalty=1, n_threads=NTHREADS)
    umi.to_csv(ck("aligned"), index=False); log(f"align done  {umi.shape}")

# 7: call alleles at the 3 cut sites
if os.path.exists(ck("call_alleles")):
    log("resume: load call_alleles"); umi = pd.read_csv(ck("call_alleles"))
else:
    umi = cas.pp.call_alleles(umi, ref_filepath=REF, barcode_interval=(20, 34),
              cutsite_locations=[112, 166, 220], cutsite_width=12, context=True, context_size=5)
    umi.to_csv(ck("call_alleles"), index=False); log(f"call_alleles done  {umi.shape}")

# 8: error-correct UMIs
if os.path.exists(ck("error_correct_umis")):
    log("resume: load error_correct_umis"); umi = pd.read_csv(ck("error_correct_umis"))
else:
    umi = cas.pp.error_correct_umis(umi, max_umi_distance=2,
              allow_allele_conflicts=False, n_threads=NTHREADS)
    umi.to_csv(ck("error_correct_umis"), index=False); log(f"error_correct_umis done  {umi.shape}")

# 9: filter molecule table
if os.path.exists(ck("filtered_molecule_table")):
    log("resume: load filtered_molecule_table"); umi = pd.read_csv(ck("filtered_molecule_table"))
else:
    umi = cas.pp.filter_molecule_table(umi, output_directory=OUTDIR,
              min_umi_per_cell=10, min_avg_reads_per_umi=2.0, min_reads_per_umi=1,   # RERUN 2026-10-04: fixed 1 (was -1 = dynamic 99th pct//10 -> 6 for MJZ019)
              intbc_prop_thresh=0.5, intbc_umi_thresh=10, intbc_dist_thresh=1,
              doublet_threshold=0.35, allow_allele_conflicts=False, plot=False)
    umi.to_csv(ck("filtered_molecule_table"), index=False); log(f"filter_molecule_table done  {umi.shape}")

# 10: call lineage groups -> FINAL allele table
allele = cas.pp.call_lineage_groups(umi, output_directory=OUTDIR,
             min_umi_per_cell=10, min_avg_reads_per_umi=2.0, min_cluster_prop=0.005,
             min_intbc_thresh=0.05, inter_doublet_threshold=0.35, kinship_thresh=0.25, plot=False)
allele.to_csv(ck("allele_table"), index=False)
log(f"DONE  allele_table {allele.shape}  cells={allele['cellBC'].nunique()}  intBCs={allele['intBC'].nunique()}")
