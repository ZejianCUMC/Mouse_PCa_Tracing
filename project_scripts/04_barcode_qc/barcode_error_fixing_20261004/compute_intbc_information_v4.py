import sys, json, glob, os, pandas as pd, numpy as np
sys.path.insert(0,"/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/scripts/barcode_heatmaps_20260925"); import gen_barcode_heatmaps as v1
sys.path.insert(0,"/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/scripts/barcode_heatmaps_20261004"); import gen_barcode_heatmaps_v3 as v3
B="/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based"
rows=[]
for proj,d in [("Npp53",f"{B}/Npp53_tracing/figures/barcode_heatmaps_v4_20261005/_display_matrices"),("Np53Rb1",f"{B}/Npp53Rb1_Tracing/figures/barcode_heatmaps_v4_20261005/_display_matrices")]:
  for f in sorted(glob.glob(f"{d}/*_display_matrix.csv")):
    s=os.path.basename(f).replace("_display_matrix.csv","")
    cm,_=v1.to_numeric_cm(pd.read_csv(f,index_col=0)); keys=json.load(open(f.replace("_display_matrix.csv","_columns.json")))
    err=v3.error_intbcs(cm,keys,v3.DEFAULTS); keep=[j for j,k in enumerate(keys) if k[:-2] not in err]
    cm=cm.iloc[:,keep]; keys=[keys[j] for j in keep]
    cm2,sites,bounds=v1.group_columns(cm,keys); N=len(cm2)
    for i,(a,b) in enumerate(bounds):
        X=cm2.iloc[:,a:b]
        det=(X>=0).any(axis=1)                                      # detected (any site observed)
        # informative site-state: an indel allele carried by >=2 cells and < 95% of detected cells (splits cells)
        inf_cell=pd.Series(False,index=X.index)
        for c in X.columns:
            v=X[c]; vc=v[v>0].value_counts(); nd=(v>=0).sum()
            good=vc[(vc>=2)&(vc<0.95*max(nd,1))].index
            inf_cell|=v.isin(good)
        # Shannon entropy of allele combination among detected cells (bits)
        combo=X[det].astype(str).agg("|".join,axis=1); p=combo.value_counts(normalize=True); H=float(-(p*np.log2(p)).sum())
        rows.append(dict(project=proj,sample=s,panel_intBC=i+1,cells=N,detected=round(det.mean(),3),edited=round((X>0).any(axis=1).mean(),3),
             informative=round(inf_cell.mean(),3),entropy_bits=round(H,2),n_states=int(p.size)))
d=pd.DataFrame(rows); d.to_csv("/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/scripts/barcode_error_report_20261004/data/intbc_information_v4.csv",index=False)
pd.set_option("display.width",200)
print(d[d["sample"]=="JZ201"].drop(columns=["project","sample","cells"]).to_string(index=False))
