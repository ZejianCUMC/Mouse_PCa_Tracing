# QC: MJZ019 re-run (min_reads_per_umi=1, fast faithful branch) vs original (dynamic threshold 6), per mouse JZ201-205,
# restricted to the same valid tumour cells used for the v8 trees (obs_metadata, Ionocyte/Cycling removed).
import pandas as pd, numpy as np
MB="/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing"
obs=pd.read_csv(f"{MB}/data/npp53_merge/obs_metadata.csv",index_col=0)
obs=obs[(obs.run=="MJZ019")&~obs.putative_state.isin(["Ionocyte","Cycling"])].copy()
obs["bc"]=obs.index.str.replace(r"^MJZ019_","",regex=True).str.replace(r"-\d+$","",regex=True)
new=pd.read_csv("fast_faithful/MJZ019/MJZ019_allele_table.csv",dtype=str,keep_default_na=False)
rows=[]; F=None
for s in ["JZ201","JZ202","JZ203","JZ204","JZ205"]:
    cells=set(obs[obs.sample_id==s].bc)
    old=pd.read_csv(f"{MB}/lineage_tracing/allele_tables/per_sample/{s}_allele_table.csv",dtype=str,keep_default_na=False)
    old=old[old.cellBC.isin(cells)]; nw=new[new.cellBC.isin(cells)]
    for tag,a in [("original",old),("rerun",nw)]:
        und=(a[["allele","r1","r2","r3"]]=="").any(axis=1)
        det=a.groupby("intBC").cellBC.nunique().sort_values(ascending=False)
        top=det.index[:min(25,len(det))]
        aa=a[a.intBC.isin(top)].drop_duplicates(["cellBC","intBC"])
        obs_sites=sum((aa[c]!="").sum() for c in ["r1","r2","r3"])
        top10=det.index[:10]; a10=a[a.intBC.isin(top10)].drop_duplicates(["cellBC","intBC"])
        o10=sum((a10[c]!="").sum() for c in ["r1","r2","r3"])
        rows.append(dict(sample=s,run=tag,rna_cells=len(cells),cells_with_alleles=a.cellBC.nunique(),
            pct_cells=round(100*a.cellBC.nunique()/len(cells),1),rows=len(a),undetected_rows=int(und.sum()),
            intBCs=len(det),intBCs_per_cell_median=float(a.groupby("cellBC").intBC.nunique().median()),
            UMIs_per_cell_median=float(a.assign(U=a.UMI.astype(int)).groupby("cellBC").U.sum().median()),
            missing_top25=round(100*(1-obs_sites/(len(cells)*len(top)*3)),1),missing_top10=round(100*(1-o10/(len(cells)*10*3)),1)))
d=pd.DataFrame(rows); d.to_csv("qc/qc_rerun_vs_original.csv",index=False)
pd.set_option("display.width",220); print(d.to_string(index=False))
