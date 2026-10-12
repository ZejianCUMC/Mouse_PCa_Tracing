# Scan every intBC used in each heatmap/tree (47 samples) for near-duplicates of a more abundant intBC.
import sys, json, glob, os, pandas as pd, numpy as np
sys.path.insert(0,"/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/scripts/barcode_heatmaps_20260925")
import gen_barcode_heatmaps as v1
B="/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based"
def lev(a,b):
    d=list(range(len(b)+1))
    for i,x in enumerate(a,1):
        p,d[0]=d[0],i
        for j,y in enumerate(b,1): p,d[j]=d[j],min(d[j]+1,d[j-1]+1,p+(x!=y))
    return d[-1]
ham=lambda a,b: sum(x!=y for x,y in zip(a,b))+abs(len(a)-len(b))
rows=[]; sets={}
for proj,d in [("Npp53",f"{B}/Npp53_tracing/figures/barcode_heatmaps_v2_20261004"),("Np53Rb1",f"{B}/Npp53Rb1_Tracing/figures/barcode_heatmaps_v2_20261004")]:
    for f in sorted(glob.glob(f"{d}/_display_matrices/*_display_matrix.csv")):
        s=os.path.basename(f).replace("_display_matrix.csv","")
        cm=pd.read_csv(f,index_col=0,dtype=str); keys=json.load(open(f.replace("_display_matrix.csv","_columns.json")))
        grp={}
        for j,k in enumerate(keys): grp.setdefault(k[:-2],[]).append(j)
        det={g:int((cm.iloc[:,js]!="-1").any(axis=1).sum()) for g,js in grp.items()}
        sets[(proj,s)]=det
        gs=sorted(det,key=lambda g:-det[g])
        for g in gs:
            par=[(h,lev(g,h)) for h in gs if det[h]>det[g] and h!=g]
            par=[(h,dd) for h,dd in par if dd<=2]
            best=min(par,key=lambda x:(x[1],-det[x[0]])) if par else (None,None)
            rows.append(dict(project=proj,sample=s,cells=len(cm),intBC=g,n_cells=det[g],pct=round(100*det[g]/len(cm),1),
                             parent=best[0],dist=best[1],parent_cells=det.get(best[0]) if best[0] else None,
                             ratio=round(det[g]/det[best[0]],3) if best[0] else None))
d=pd.DataFrame(rows); d.to_csv("intbc_scan.csv",index=False)
pd.set_option("display.width",250)
fl=d[d.parent.notna()]
print("total intBC entries:",len(d)," flagged (Levenshtein<=2 to a more-detected intBC in same mouse):",len(fl))
print(fl.groupby(["project","sample","dist"]).size().unstack(fill_value=0).to_string())
print("\nflagged ratio (child/parent cells) summary:"); print(fl.ratio.describe().round(3).to_string())
print("\nall flagged entries:"); print(fl[["project","sample","intBC","n_cells","pct","parent","dist","parent_cells","ratio"]].to_string(index=False))
# clean (unflagged) intBC sets per project/library
print("\nclean intBC count per sample:"); print(d[d.parent.isna()].groupby(["project","sample"]).size().to_string())
