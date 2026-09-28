# Mouse_PCa_Tracing
This repository saves crucial analysis workflow to handle single cell level lineage tracing works in mouse PCa.

It goes from raw 10x scRNA-seq + hashtag + Cas9 lineage-recorder FASTQs to **tumour cell states** and a
**Cassiopeia greedy lineage tree per mouse**.

**Start here → [HANDSOFF.md](HANDSOFF.md)**, the step-by-step tutorial. Either start from FASTQs on the c2b2 HPC, or start at step 03 on your own machine with the data folder Zejian provides.

```
01 Cell Ranger → 02 CellBender → 03 HTODemux + QC + DoubletFinder → 04 tumour cells + cell states
05 Cassiopeia (lineage FASTQ → allele table) ───────────────────────────────┐
                                                                            └→ 06 greedy tree per mouse
```

This repository contains **code only**. No sequencing data, metadata, intermediate files, trees or
cell-state results are stored here (enforced by `.gitignore`).
