"""Current FSU / FSD gene sets per tumor model for the progression GSEA (Zejian 2026-10-12: re-produce the Npp53 part-2 Fig 9.1 GSEA "based on
current FSU/FSD for each tumor model"): the normalized-count meta-analysis lists FSU_top1000 / FSD_top1000 (arms metaNorm_20261010_<model>_bcr),
mapped to human with the pipeline's orthologs.py (MGI homology classes, all human members, HGNC rescue). Writes <model>_{FSU,FSD}_human.txt. Env: moscot."""
import sys, json
from pathlib import Path
O = Path(__file__).resolve().parent; O5 = O.parents[1]; sys.path.insert(0, str(O5)); import orthologs as ORT
for g in ["Npp53", "Npp53Rb1"]:
    L = json.load(open(O5 / f"metaNorm_20261010_{g}_bcr" / "fsu_lists.json"))
    for d in ["FSU", "FSD"]:
        m = L[f"META_{g}_{d}_top1000_n3000"]["genes"]; h = ORT.to_human(m)
        (O / f"{g}_{d}_human.txt").write_text("\n".join(h) + "\n"); print(g, d, "mouse", len(m), "human", len(h))
