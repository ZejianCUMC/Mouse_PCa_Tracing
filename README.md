# Mouse_PCa_Tracing
This repository saves crucial analysis workflow to handle single cell level lineage tracing works in mouse PCa.

All scripts and commands are in `scripts/`:

- **`scripts/tutorial/`**: a teaching workflow from raw 10x scRNA-seq + hashtag + Cas9 lineage-recorder FASTQs
  to **tumour cell states** and a **Cassiopeia greedy lineage tree per mouse**. **Start here →
  [HANDSOFF.md](HANDSOFF.md)**, the step-by-step tutorial. Either start from FASTQs on the c2b2 HPC, or start at
  step 03 on your own machine with the data folder Zejian provides. Run every command from `scripts/tutorial/`.

  ```
  01 Cell Ranger → 02 CellBender → 03 HTODemux + QC + DoubletFinder → 04 tumour cells + cell states
  05 Cassiopeia (lineage FASTQ → allele table) ───────────────────────────────┐
                                                                              └→ 06 greedy tree per mouse
  ```

- **`scripts/project/`**: the current analysis scripts of the whole tracing project (Npp53, Npp53Rb1 and the
  merged cohort), organised by stage: demultiplexing, cell states, Cassiopeia trees, barcode QC, fitness
  signatures, survival, evolution features and CARTA fate maps. Start with
  [scripts/project/README.md](scripts/project/README.md).
